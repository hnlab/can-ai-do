"""output distribution of properties and fingerprint.
"""
import gzip
import json
import pickle
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import scipy.sparse as sp
from scipy.spatial import distance

from tqdm import tqdm
# from multiprocessing.dummy import Pool
import multiprocessing as mp

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.DataStructs import BulkTanimotoSimilarity

from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
PandasTools.RenderImagesInAllDataFrames(images=True)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f',
                    '--fold_list',
                    required=True,
                    help="k-fold config in json, a dict or a list of list")
parser.add_argument('-d',
                    '--datadir',
                    default='./all',
                    help="datadir, default is ./all")
parser.add_argument('--use_dude_ism', action='store_true')
parser.add_argument('--use_dude_sdf', action='store_true')
parser.add_argument('--MW500',
                    action='store_true',
                    help="remove actives with HeavyAtomMolWt > 500.")
parser.add_argument('-o',
                    '--output',
                    default='result.jpg',
                    help="distribution figures")
parser.add_argument('--single', action='store_true')
parser.add_argument(
    '--feat_imports',
    nargs='+',
    help="random forest results in .json format having feature_importances")
parser.add_argument(
    '--zinc',
    required=True,
    help='ZINC .smi from http://zinc12.docking.org/db/bysubset/6/6_p0.smi.gz')
args = parser.parse_args()

nBits = 2048


def mfp2(m):
    # radius 2 MorganFingerprint equal ECFP4
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=nBits)
    return fp


def getProp(mol):
    # mw = Descriptors.ExactMolWt(mol)
    mwha = Descriptors.HeavyAtomMolWt(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    q = Chem.GetFormalCharge(mol)
    return tuple([mwha, mw, logp, rotb, hbd, hba, q])


def ForwardMol2MolSupplier(file_obj, sanitize=True):
    lines = []
    for line in file_obj:
        if line.startswith(b"@<TRIPOS>MOLECULE"):
            if lines:
                block = b''.join(lines)
                yield Chem.MolFromMol2Block(block, sanitize=sanitize)
                lines = []
        lines.append(line)
    else:
        if lines:
            block = b''.join(lines)
            yield Chem.MolFromMol2Block(block, sanitize=sanitize)
    file_obj.close()


def load_smiles(names, MolWt=None, MW500=False, fpAsArray=False, bits=None):
    datadir = Path(args.datadir)
    all_fps = []
    all_props = []
    all_labels = []
    for name in names:
        tdir = datadir / name

        activeFile = tdir / 'actives_final.smi'
        if activeFile.exists():
            # generate in this work
            active_supp = Chem.SmilesMolSupplier(str(activeFile),
                                                 titleLine=False)
        else:
            # from DUD-E
            if args.use_dude_ism:  # no charge and hydrogens information.
                activeFile = tdir / 'actives_final.ism'
                active_supp = Chem.SmilesMolSupplier(str(activeFile),
                                                     titleLine=False)
            elif args.use_dude_sdf:  # duplicate
                activeFile = tdir / 'actives_final.sdf.gz'
                active_supp = Chem.ForwardSDMolSupplier(gzip.open(activeFile))
            else:
                activeFile = tdir / 'actives_final.mol2.gz'
                active_supp = ForwardMol2MolSupplier(gzip.open(activeFile))

        decoyFile = tdir / 'decoys_final.smi'
        if decoyFile.exists():
            # generate in this work
            decoy_supp = Chem.SmilesMolSupplier(str(decoyFile),
                                                titleLine=False)
        else:
            # from DUD-E
            if args.use_dude_ism:
                decoyFile = tdir / 'decoys_final.ism'
                decoy_supp = Chem.SmilesMolSupplier(str(decoyFile),
                                                    titleLine=False)
            elif args.use_dude_sdf:
                decoyFile = tdir / 'decoys_final.sdf.gz'
                decoy_supp = Chem.ForwardSDMolSupplier(gzip.open(decoyFile))
            else:
                decoyFile = tdir / 'decoys_final.mol2.gz'
                decoy_supp = ForwardMol2MolSupplier(gzip.open(decoyFile))

        fpf = activeFile.with_name(activeFile.name + '.fp.pkl')
        propf = activeFile.with_name(activeFile.name + '.prop.pkl')
        labelf = activeFile.with_name(activeFile.name + '.label.pkl')
        fpf_mw500 = fpf.with_suffix('.MW500.pkl')
        propf_mw500 = propf.with_suffix('.MW500.pkl')
        labelf_mw500 = labelf.with_suffix('.MW500.pkl')
        if not all([f.exists() for f in (fpf, propf, labelf)]):
            fps = []
            props = []
            labels = []
            fps_mw500 = []
            props_mw500 = []
            labels_mw500 = []
            canonical_smiles_set = set()
            for m in active_supp:
                if m is None: continue
                smiles = Chem.MolToSmiles(m)
                if smiles in canonical_smiles_set:
                    continue
                else:
                    canonical_smiles_set.add(smiles)
                fp = mfp2(m)
                fps.append(fp)
                p = getProp(m)
                props.append(p)
                labels.append(1)
                # p:[mwha, mw, logp, rotb, hbd, hba, q]
                if p[0] > 500:
                    continue
                fps_mw500.append(fp)
                props_mw500.append(p)
                labels_mw500.append(1)
            frac = len(fps_mw500) / len(fps)
            decoy_mols = [m for m in decoy_supp if m is not None]
            select_num = int(frac * len(decoy_mols))
            np.random.seed(123)
            inds = np.random.choice(len(decoy_mols), select_num, replace=False)
            for i, m in enumerate(decoy_mols):
                fp = mfp2(m)
                fps.append(fp)
                p = getProp(m)
                props.append(p)
                labels.append(0)
                if i in inds:
                    fps_mw500.append(fp)
                    props_mw500.append(p)
                    labels_mw500.append(0)

            with open(fpf, 'wb') as f:
                pickle.dump(fps, f)
            with open(propf, 'wb') as f:
                pickle.dump(props, f)
            with open(labelf, 'wb') as f:
                pickle.dump(labels, f)

            with open(fpf_mw500, 'wb') as f:
                pickle.dump(fps_mw500, f)
            with open(propf_mw500, 'wb') as f:
                pickle.dump(props_mw500, f)
            with open(labelf_mw500, 'wb') as f:
                pickle.dump(labels_mw500, f)

            for file_name, fps_list in ((fpf, fps), (fpf_mw500, fps_mw500)):
                fpf_np = file_name.with_suffix('.np.pkl')
                fps_np = []
                for fp in fps_list:
                    fp_np = np.zeros((0, ), dtype=np.int8)
                    # faster, https://github.com/rdkit/rdkit/pull/2557
                    DataStructs.ConvertToNumpyArray(fp, fp_np)
                    fps_np.append(fp_np)
                fps_np = np.array(fps_np, dtype=np.int8)
                with open(fpf_np, 'wb') as f:
                    pickle.dump(fps_np, f)

        if MW500:
            fpf = fpf_mw500
            propf = propf_mw500
            labelf = labelf_mw500

        if bits is not None:
            fpAsArrays = True

        if fpAsArray:
            fpf_np = fpf.with_suffix('.np.pkl')
            with open(fpf_np, 'rb') as f:
                fps = pickle.load(f)
        else:
            with open(fpf, 'rb') as f:
                fps = pickle.load(f)

        if bits is not None:
            fps = fps[:, bits]

        with open(propf, 'rb') as f:
            props = pickle.load(f)
        with open(labelf, 'rb') as f:
            labels = pickle.load(f)
        all_fps.extend(fps)
        all_props.extend(props)
        all_labels.extend(labels)

    # prop: [mwha, mw, logp, rotb, hbd, hba, q]
    all_props = np.array(all_props)
    if MolWt == 'HeavyAtomMolWt':
        all_props = all_props[:, (0, 2, 3, 4, 5, 6)]
    if MolWt == 'MolWt':
        all_props = all_props[:, (1, 2, 3, 4, 5, 6)]
    return all_fps, all_props, all_labels


with open(args.fold_list) as f:
    folds = json.load(f)
    if type(folds) is list:
        folds = {'{}'.format(fold): fold for fold in folds}
    targets = [i for fold in folds.values() for i in fold]

iter_targets = [[i] for i in targets]
p = mp.Pool()
for _ in tqdm(p.imap_unordered(load_smiles, iter_targets),
              desc='Converting smiles into fingerprints and properties',
              total=len(targets)):
    pass
p.close()

output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)
print(f"loading fps and properties with MW500={args.MW500}")
fps, props, labels = load_smiles(targets, MW500=args.MW500)
props = np.array(props)
labels = np.array(labels)
active_mask = labels == 1
decoy_mask = labels == 0
# mw, logp, rotb, hbd, hba, q
prop_keys = ['mwha', 'mw', 'logp', 'rotb', 'hbd', 'hba', 'q']
prop_names = {
    'mwha': 'Molecular Weight of Heavy Atoms (Ignoring Hydrogens)',
    'mw': 'Molecular Weight',
    'logp': 'Calculated LogP',
    'rotb': 'Number of Rotatable Bonds',
    'hbd': 'Number of Hydrogen Bond Donors',
    'hba': 'Number of Hydrogen Bond Acceptors',
    'q': 'Net Charge'
}
prop_bins = {
    'mwha': np.linspace(100, 700, 61),
    'mw': np.linspace(100, 700, 61),
    'logp': np.linspace(-8, 10, 46),
    'rotb': np.linspace(0, 20, 21) + 0.5,
    'hbd': np.linspace(0, 20, 21) + 0.5,
    'hba': np.linspace(0, 20, 21) + 0.5,
    'q': np.linspace(-4, 3, 8) + 0.5
}
fig, axes = plt.subplots(nrows=7, ncols=1, figsize=(8, 9))
fig.subplots_adjust(hspace=.8, top=0.95, bottom=0.05)
for p_key, ps, ax in zip(prop_keys, props.T, axes):
    hist_kws = None
    if p_key in ['logp']:
        ax.set_xticks(np.linspace(-8, 10, 10))
    if p_key in ['rotb', 'hbd', 'hba']:
        hist_kws = {'rwidth': 0.6}
        ax.set_xticks(np.linspace(0, 20, 11))
    if p_key in ['q']:
        hist_kws = {'rwidth': 0.2}
        ax.set_xticks(np.linspace(-4, 3, 8))
    bins = prop_bins[p_key]
    # ax.set_title(targets)
    ax.set_xlabel(prop_names[p_key])
    print(f"{p_key}: min {min(ps)} max {max(ps)}")
    decoy = ps[decoy_mask]
    sns.distplot(decoy,
                 label='Decoys',
                 bins=bins,
                 kde=False,
                 norm_hist=True,
                 color='blue',
                 hist_kws=hist_kws,
                 ax=ax)
    active = ps[active_mask]
    sns.distplot(active,
                 label='Actives',
                 bins=bins,
                 kde=False,
                 norm_hist=True,
                 color='red',
                 hist_kws=hist_kws,
                 ax=ax)
    ax.legend()

fig.savefig(output, dpi=300)
print(f"result figure saved at {args.output}")

if args.single:
    for target in targets:
        fps, props, labels = load_smiles([target], MW500=args.MW500)
        props = np.array(props)
        labels = np.array(labels)
        active_mask = labels == 1
        decoy_mask = labels == 0
        # mw, logp, rotb, hbd, hba, q
        prop_keys = ['mwha', 'mw', 'logp', 'rotb', 'hbd', 'hba', 'q']
        prop_names = {
            'mwha': 'Molecular Weight of Heavy Atoms (Ignoring Hydrogens)',
            'mw': 'Molecular Weight',
            'logp': 'Calculated LogP',
            'rotb': 'Number of Rotatable Bonds',
            'hbd': 'Number of Hydrogen Bond Donors',
            'hba': 'Number of Hydrogen Bond Acceptors',
            'q': 'Net Charge'
        }
        prop_bins = {
            'mwha': np.linspace(100, 700, 61),
            'mw': np.linspace(100, 700, 61),
            'logp': np.linspace(-8, 10, 46),
            'rotb': np.linspace(0, 20, 21) + 0.5,
            'hbd': np.linspace(0, 20, 21) + 0.5,
            'hba': np.linspace(0, 20, 21) + 0.5,
            'q': np.linspace(-4, 3, 8) + 0.5
        }
        fig, axes = plt.subplots(nrows=7, ncols=1, figsize=(8, 9))
        fig.subplots_adjust(hspace=.8, top=0.95, bottom=0.05)
        for p_key, ps, ax in zip(prop_keys, props.T, axes):
            hist_kws = None
            if p_key in ['logp']:
                ax.set_xticks(np.linspace(-8, 10, 10))
            if p_key in ['rotb', 'hbd', 'hba']:
                hist_kws = {'rwidth': 0.6}
                ax.set_xticks(np.linspace(0, 20, 11))
            if p_key in ['q']:
                hist_kws = {'rwidth': 0.2}
                ax.set_xticks(np.linspace(-4, 3, 8))
            bins = prop_bins[p_key]
            # ax.set_title(targets)
            ax.set_xlabel(prop_names[p_key])
            # print(f"{p_key}: min {min(ps)} max {max(ps)}")
            decoy = ps[decoy_mask]
            sns.distplot(decoy,
                         label='Decoys',
                         bins=bins,
                         kde=False,
                         norm_hist=True,
                         color='blue',
                         hist_kws=hist_kws,
                         ax=ax)
            active = ps[active_mask]
            sns.distplot(active,
                         label='Actives',
                         bins=bins,
                         kde=False,
                         norm_hist=True,
                         color='red',
                         hist_kws=hist_kws,
                         ax=ax)
            ax.legend()
        target_output = output.with_suffix(f'.{target}.jpg')
        fig.savefig(target_output, dpi=300)
        plt.close(fig)  # remove warning about too many opened figures
        print(f"result figure saved at {target_output}")


def smiles2fp(line):
    smiles, name, *_ = line.split()
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    fp = mfp2(m)
    return fp


zinc = Path(args.zinc)
zinc_freq_file = zinc.with_suffix('.bits_freq.json')
if zinc_freq_file.exists():
    with open(zinc_freq_file) as f:
        zinc_freq = json.load(f)
    print(f'load zinc frequency from {zinc_freq_file}')
else:
    # http://rdkit.blogspot.com/2016/02/morgan-fingerprint-bit-statistics.html
    freq = [0 for i in range(nBits)]
    vaild_count = 0
    with open(args.zinc) as f:
        total = sum(1 for line in f)
    pbar = tqdm(desc='count ZINC bits freq', total=total, unit=' mol')
    with open(args.zinc) as f:
        # 5000000 chars about 100000 lines
        nchar = 5000000
        lines = f.readlines(nchar)
        while lines:
            with mp.Pool() as p:
                fps = p.map(smiles2fp, lines)
            for fp in fps:
                if fp is None:
                    continue
                vaild_count += 1
                for i in fp.GetOnBits():
                    freq[i] += 1
            pbar.update(len(lines))
            lines = f.readlines(nchar)
    pbar.close()

    zinc_freq = np.array(freq) / vaild_count
    with open(zinc_freq_file, 'w') as f:
        json.dump(zinc_freq.tolist(), f)
    print(f'zinc freq saved to {zinc_freq_file}')

print('counting DUD-E fingerprint bits ...')
fps, probe, labels = load_smiles(targets, MW500=args.MW500, fpAsArray=True)
fps = np.asarray(fps)
labels = np.asarray(labels)
active_fps = fps[labels == 1]
decoy_fps = fps[labels == 0]
active_bits_freq = active_fps.sum(axis=0) / len(active_fps)
decoy_bits_freq = decoy_fps.sum(axis=0) / len(decoy_fps)
bits_factor = np.divide(active_bits_freq,
                        decoy_bits_freq,
                        out=np.zeros(nBits),
                        where=decoy_bits_freq != 0)
bits_factor = np.log2(bits_factor, out=np.zeros(nBits), where=bits_factor != 0)
mean_freq = (active_bits_freq + decoy_bits_freq) / 2

freq_cutoff = 0.03
factor_cutoff = 1
freq_mask = mean_freq >= freq_cutoff
pos_mask = bits_factor >= factor_cutoff
neg_mask = bits_factor <= -factor_cutoff
factor_mask = pos_mask | neg_mask
sign_mask = freq_mask & factor_mask

sizes = np.array([4 for i in bits_factor])
sizes[sign_mask] = 16
colors = np.array(['black' for i in bits_factor])
colors[freq_mask & pos_mask] = 'red'
colors[freq_mask & neg_mask] = 'blue'
fig, ax = plt.subplots()
ax.vlines(factor_cutoff, -1, 1, colors='gray', zorder=1)
ax.vlines(-factor_cutoff, -1, 1, colors='gray', zorder=1)
ax.hlines(freq_cutoff, -6, 6, colors='gray', zorder=1)
ax.scatter(bits_factor, mean_freq, s=sizes, c=colors, zorder=2)
ax.set_xlim(-5, 5)
ax.set_xticks(np.linspace(-4, 4, 9))
ax.set_ylim(-0.02, 1)
ax.set_yticks(np.append(np.linspace(0, 1, 6), freq_cutoff))
ax.set_xlabel('log2(active_bit_freq / decoy_bit_freq)')
ax.set_ylabel('(active_bit_freq + decoy_bit_freq) / 2')
jpg = output.with_suffix(f'.bits_freq_vs_factor.jpg')
fig.savefig(jpg, dpi=300)
ax.set_ylim(-0.02, 0.3)
ax.set_yticks(np.append(np.linspace(0, 0.3, 4), freq_cutoff))
print(f"bits freq vs factor saved at {jpg}")
jpg = output.with_suffix(f'.bits_freq0.4_vs_factor.jpg')
fig.savefig(jpg, dpi=300)
print(f"bits freq vs factor saved at {jpg}")
plt.close(fig)

data = list(
    zip(range(nBits), zinc_freq, bits_factor, mean_freq, active_bits_freq,
        decoy_bits_freq))
cols = [
    'bit', 'zinc_freq', 'log2_FC', 'mean_freq', 'active_freq', 'decoy_freq'
]
df = pd.DataFrame(data, columns=cols)
df = df.sort_values(by='zinc_freq', ascending=False)
df = df.reset_index(drop=True)
csv = output.with_suffix(f'.bits_freq.csv')
with open(csv, 'w') as f:
    df.to_csv(f, index=False)
print(f"bits data saved at {csv}")

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(df.index, df.active_freq, color='red', label='Actives', marker='+')
ax.scatter(df.index, df.decoy_freq, color='blue', label='Decoys', marker='+')
ax.scatter(df.index, df.zinc_freq, color='black', label='ZINC', marker='+')
# ax.legend()
ax.set_ylim([-0.05, 1.05])
ax.set_yticks(np.linspace(0, 1, 11))
ax.set_xlabel('2048 bits of Morgan Fingerprint sorted by frequency in ZINC')
ax.set_ylabel('frequency')
ax2 = fig.add_axes([0.29, 0.365, 0.6, 0.5])
ax2.scatter(df.index, df.active_freq, color='red', label='Actives', marker='+')
ax2.scatter(df.index, df.decoy_freq, color='blue', label='Decoys', marker='+')
ax2.scatter(df.index, df.zinc_freq, color='black', label='ZINC', marker='+')
ax2.legend()
ax2.set_ylim([-0.005, 0.1])
jpg = output.with_suffix('.all_bits_vs_zinc.jpg')
fig.savefig(jpg, dpi=300)
print(f'all bits saved to {jpg}')

df = df.loc[(df.mean_freq >= 0.03) & (np.abs(df.log2_FC) >= 1)]
df = df.reset_index(drop=True)
fig, ax = plt.subplots(figsize=(8, 6))
for i, (lower, upper) in enumerate(zip(df.active_freq, df.decoy_freq)):
    if lower > upper:
        lower, upper = upper, lower
    ax.vlines(i, lower, upper, colors='silver', zorder=0)
ax.scatter(df.index, df.active_freq, color='red', label='Actives', marker='+')
ax.scatter(df.index, df.decoy_freq, color='blue', label='Decoys', marker='+')
ax.scatter(df.index, df.zinc_freq, color='black', label='ZINC', marker='+')
ax.set_xlabel(f'{len(df)} high frequency and change significant bits in DUD-E')
ax.set_ylabel('frequency')
ax.set_xticks([i * 5 for i in range(int(len(df) / 5) + 1)])
ax.legend()
jpg = output.with_suffix('.significant_bits_vs_zinc.jpg')
fig.savefig(jpg, dpi=300)
print(f'significant bits plot saved to {jpg}')

datadir = Path(args.datadir)
mol2_files = []
for t in targets:
    mol2_files.append(datadir / t / 'actives_final.mol2.gz')
    mol2_files.append(datadir / t / 'decoys_final.mol2.gz')

bits_data = []
np.random.seed(123)
n_example = 4
bits_examples = [[] for i in range(n_example)]
for bit in tqdm(df.bit, desc='bits example'):
    bit_examples = []
    example_count = 0
    for file_name in np.random.permutation(mol2_files):
        if example_count >= n_example:
            break
        with gzip.open(file_name) as f:
            for m in ForwardMol2MolSupplier(f):  # not Kekulize
                if m is None:
                    continue
                # remove 3D info for Draw Bits, may fail for Kekulize
                m = Chem.MolFromSmiles(Chem.MolToSmiles(m))
                if m is None:
                    continue
                info = {}
                fp = AllChem.GetMorganFingerprintAsBitVect(m,
                                                           2,
                                                           nBits=nBits,
                                                           bitInfo=info)
                if bit in set(fp.GetOnBits()):
                    bit_examples.append((m, bit, info))
                    example_count += 1
                    break
    # http://rdkit.blogspot.com/2018/10/using-new-fingerprint-bit-rendering-code.html
    examples_svg = [Draw.DrawMorganBit(*i) for i in bit_examples]
    # remove two '\n' at the start and end of svg
    examples_svg = [i.strip() for i in examples_svg]
    examples_svg = [i.replace('\n<svg', '<svg') for i in examples_svg]
    for list_, svg in zip(bits_examples, examples_svg):
        list_.append(svg)

svg_cols = [f'example_{i+1}' for i in range(4)]
for col, list_ in zip(svg_cols, bits_examples):
    df[col] = list_
df = df.sort_values(by='zinc_freq', ascending=False)
df = df.reset_index(drop=True)
df['bit'] = [f'bit:<br>{i}' for i in df['bit']]
df['log2_FC'] = [f'log2(FC):<br>{i:.2f}' for i in df['log2_FC']]
for col in ['active_freq', 'decoy_freq', 'mean_freq', 'zinc_freq']:
    df[col] = [f'{col}:<br>{100*i:.1f}%' for i in df[col]]

html_file = output.with_suffix(".bits_example.html")
# https://stackoverflow.com/questions/50807744/apply-css-class-to-pandas-dataframe-using-to-html/50939211
style = """
/* includes alternating gray and white with on-hover color */
svg { 
    width: 150px;
    height: 150px;
    display: block;
}
.dataframe {
    font-size: 11pt; 
    font-family: Arial;
    /* border-collapse: collapse; */
    border: 1px solid silver;

}
.dataframe td, th {
    text-align: right;
    padding: 2px;
}
.dataframe tr:nth-child(even) {
    background: #E0E0E0;
}
.dataframe tr:hover {
    background: silver;
    cursor: pointer;
}
"""
table = df.to_html()
with open(html_file, 'w') as f:
    f.write(f'<html>\n'
            f'<head><title>HTML Pandas Dataframe with CSS</title>\n'
            f'<style>{style}</style>\n'
            f'</head>\n'
            f'<body>\n{table}\n</body>\n'
            f'</html>')
print(f'bits example saved at {html_file}')

if args.feat_imports:
    prop_data = []
    prop_cols = ('file', 'fold_name', 'mwha', 'logp', 'rotb', 'hbd', 'hba',
                 'q')
    fp_data = []
    top_n = 5
    fp_cols = ('file', 'fold_name', *np.arange(top_n))
    import_bits = set()
    for feat_import_file in args.feat_imports:
        feat_import_file = Path(feat_import_file)
        results = json.load(open(feat_import_file))
        repeat_results, repeat_means = results
        repeat_result = repeat_results[0]
        feature_sets = ('fp', 'prop')
        folds = repeat_result['folds']
        feat_set_imports = {
            'fp': repeat_result['fp']['feature_importances'],
            'prop': repeat_result['prop']['feature_importances']
        }
        feat_imports = feat_set_imports['prop']
        fold_weights = np.array([len(i) for i in folds.values()])
        # targets in folds is in test set, need number of targets in training set.
        fold_weights = sum(fold_weights) - fold_weights
        fold_weights = np.array(fold_weights) / sum(fold_weights)
        fold_weights = fold_weights.reshape(-1, 1)
        for fold_name, fold_feat_imports in zip(folds, feat_imports):
            feat_percents = [f'{100*i:.2f}%' for i in fold_feat_imports]
            prop_data.append(
                [feat_import_file.stem, fold_name, *feat_percents])
        mean_feat_imports = np.sum(feat_imports * fold_weights, axis=0)
        feat_percents = [f'{100*i:.2f}%' for i in mean_feat_imports]
        prop_data.append([feat_import_file.stem, 'mean', *feat_percents])
        feat_imports = feat_set_imports['fp']
        for fold_name, fold_feat_imports in zip(folds, feat_imports):
            fold_feat_imports = np.array(fold_feat_imports)
            sort_idx = np.argsort(-fold_feat_imports)[:top_n]
            import_bits.update(sort_idx)
            feat_percents = [
                f'{i}({100*fold_feat_imports[i]:.2f}%)' for i in sort_idx
            ]
            fp_data.append([feat_import_file.stem, fold_name, *feat_percents])
        mean_feat_imports = np.sum(feat_imports * fold_weights, axis=0)
        sort_idx = np.argsort(-mean_feat_imports)[:top_n]
        import_bits.update(sort_idx)
        feat_percents = [
            f'{i}({100*mean_feat_imports[i]:.2f}%)' for i in sort_idx
        ]
        fp_data.append([feat_import_file.stem, 'mean', *feat_percents])

    df = pd.DataFrame(prop_data, columns=prop_cols)
    csv = output.with_suffix('.prop_imports.csv')
    df.to_csv(csv, index=False)
    print(f"feature_importances saved at {csv}")
    df = pd.DataFrame(fp_data, columns=fp_cols)
    csv = output.with_suffix('.fp_imports.csv')
    df.to_csv(csv, index=False)
    print(f"feature_importances saved at {csv}")
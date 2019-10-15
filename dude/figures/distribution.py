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


def count_bits(targets):
    fps, _, labels = load_smiles(targets, MW500=args.MW500)
    active_count = 0
    decoy_count = 0
    active_bits_count = np.zeros(nBits, dtype=int)
    decoy_bits_count = np.zeros(nBits, dtype=int)
    for fp, label in zip(fps, labels):
        if label == 1:
            active_count += 1
            for bit in fp.GetOnBits():
                active_bits_count[bit] += 1
        else:
            decoy_count += 1
            for bit in fp.GetOnBits():
                decoy_bits_count[bit] += 1
    return (active_count, decoy_count, active_bits_count, decoy_bits_count)


print('counting fingerprint bits ...')
with mp.Pool() as p:
    iter_targets = [[i] for i in targets]
    counters = p.map(count_bits, iter_targets)

active_num = 0
decoy_num = 0
active_bits_count = np.zeros(nBits, dtype=int)
decoy_bits_count = np.zeros(nBits, dtype=int)
for c in counters:
    active_num += c[0]
    decoy_num += c[1]
    active_bits_count += c[2]
    decoy_bits_count += c[3]
active_bits_freq = active_bits_count / active_num
decoy_bits_freq = decoy_bits_count / decoy_num

high_freq_idx = np.flatnonzero(decoy_bits_freq >= 1 / nBits)
mean_freq = (active_bits_freq + decoy_bits_freq) / 2
mean_high_freq = mean_freq[high_freq_idx]
bits_factor = active_bits_freq[high_freq_idx] / decoy_bits_freq[high_freq_idx]
bits_factor = np.log2(bits_factor)

fig, ax = plt.subplots()
sns.distplot(
    bits_factor,
    # bins=bins,
    kde=False,
    norm_hist=True,
    # color='blue',
    # hist_kws=hist_kws,
    ax=ax)
jpg = output.with_suffix(f'.bits_factor.jpg')
fig.savefig(jpg, dpi=300)
plt.close(fig)
print(f"bits factor distribution saved at {jpg}")

fig, ax = plt.subplots()
ax.scatter(bits_factor, mean_high_freq)
jpg = output.with_suffix(f'.bits_freq_vs_factor.jpg')
fig.savefig(jpg, dpi=300)
plt.close(fig)
print(f"bits freq vs factor saved at {jpg}")

significant_bits = high_freq_idx[(np.abs(bits_factor) >= 1)
                                 & (mean_high_freq >= 0.01)]
jsonf = output.with_suffix(f'.significant_bits.json')
with open(jsonf, 'w') as f:
    json.dump(significant_bits.tolist(), f)
print(f"significant_bits saved at {jsonf}")

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

    bit_data = []
    import_bits = sorted(import_bits)
    for bit in import_bits:
        freq_in_active = active_bits_freq[bit]
        freq_in_decoy = decoy_bits_freq[bit]
        bit_data.append((bit, freq_in_active, freq_in_decoy))
    df = pd.DataFrame(bit_data,
                      columns=['bit', 'freq_in_active', 'freq_in_decoy'])
    csv = output.with_suffix('.import_bits.csv')
    df.to_csv(csv, index=False)
    print(f"important bits saved at {csv}")

    from rdkit.Chem import Draw

    datadir = Path(args.datadir)
    files = []
    for t in targets:
        files.append(datadir / t / 'actives_final.mol2.gz')
        files.append(datadir / t / 'decoys_final.mol2.gz')

    bit_example_path = output.with_suffix('.bit_examples')
    bit_example_path.mkdir(exist_ok=True)
    for bit in import_bits:
        bit_examples = []
        example_count = 0
        for file_name in np.random.permutation(files):
            if example_count >= 9:
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
        svg_text = Draw.DrawMorganBits(bit_examples,
                                       molsPerRow=3)  # legends = []
        svg_file = bit_example_path / f"bit{bit:04d}_example.svg"
        svg_file.write_text(svg_text)

"""output distribution of properties and fingerprint.
"""
import gzip
import json
import pickle
import argparse
import numpy as np
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
args = parser.parse_args()


def mfp2(m):
    # radius 2 MorganFingerprint equal ECFP4
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=512)
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


def load_smiles(names, MolWt=None, MW500=False):
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

        if MW500:
            fpf = fpf_mw500
            propf = propf_mw500
            labelf = labelf_mw500

        with open(fpf, 'rb') as f:
            fps = pickle.load(f)
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
#%%
if args.single:
    for target in targets:
        fps, props, labels = load_smiles([target])
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
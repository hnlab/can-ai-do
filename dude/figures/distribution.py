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
parser.add_argument('--removeHeavyMW500',
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


def load_smiles(names,
                HeavyAtomMolWt=True,
                removeHeavyMW500=False,
                subsample=False):
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
        fpf = activeFile.with_name(activeFile.name + '.fp.MWHA.pkl')
        propf = activeFile.with_name(activeFile.name + '.prop.MWHA.pkl')
        labelf = activeFile.with_name(activeFile.name + '.labelf.MWHA.pkl')
        if fpf.exists() and propf.exists() and labelf.exists():
            with open(fpf, 'rb') as f:
                fps = pickle.load(f)
            with open(propf, 'rb') as f:
                props = pickle.load(f)
            with open(labelf, 'rb') as f:
                labels = pickle.load(f)
        else:
            fps = []
            props = []
            labels = []
            canonical_smiles_set = set()
            for m in active_supp:
                if m is None: continue
                smiles = Chem.MolToSmiles(m)
                if smiles in canonical_smiles_set:
                    continue
                else:
                    canonical_smiles_set.add(smiles)
                fps.append(mfp2(m))
                props.append(getProp(m))
                labels.append(1)
            for m in decoy_supp:
                if m is None: continue
                fps.append(mfp2(m))
                props.append(getProp(m))
                labels.append(0)
            with open(fpf, 'wb') as f:
                pickle.dump(fps, f)
            with open(propf, 'wb') as f:
                pickle.dump(props, f)
            with open(labelf, 'wb') as f:
                pickle.dump(labels, f)
        if removeHeavyMW500:
            props = np.array(props)
            labels = np.array(labels)
            fps = np.array(fps)
            actives_props = props[labels == 1]
            decoys_props = props[labels == 0]
            actives_fps = fps[labels == 1]
            decoys_fps = fps[labels == 0]
            big_active_mask = actives_props[:, 0] > 500
            nbig = sum(big_active_mask)
            fraction = nbig / len(actives_props)
            # print(f"left {len(actives) - nbig:4d}/{len(actives):4d} from target {name}")
            if subsample:
                perm = np.random.permutation(len(actives_props))
                actives_props = actives_props[perm[nbig:]]
                actives_fps = actives_fps[perm[nbig:]]
            else:
                actives_props = actives_props[~big_active_mask]
                actives_fps = actives_fps[~big_active_mask]
            perm = np.random.permutation(len(decoys_props))
            rmN = int(fraction * len(decoys_props))
            decoys_props = decoys_props[perm[rmN:]]
            decoys_fps = decoys_fps[perm[rmN:]]
            props = np.vstack((actives_props, decoys_props))
            fps = np.vstack((actives_fps, decoys_fps))
            labels = np.hstack(
                (np.ones(len(actives_props)), np.zeros(len(decoys_props))))
        all_fps.extend(fps)
        all_props.extend(props)
        all_labels.extend(labels)
    # [mwha, mw, logp, rotb, hbd, hba, q]
    all_props = np.array(all_props)
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
print(f"loading fps and properties with removeHeavyMW500={args.removeHeavyMW500}")
fps, props, labels = load_smiles(targets, removeHeavyMW500=args.removeHeavyMW500)
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
        print(f"result figure saved at {target_output}")
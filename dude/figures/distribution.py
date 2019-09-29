#%%
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

#%%
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f',
                    '--fold_list',
                    required=True,
                    help="k-fold config in json, a dict or a list of list")
parser.add_argument('-d',
                    '--datadir',
                    required=True,
                    default='./all',
                    help="datadir, default is ./all")
parser.add_argument('-o',
                    '--output',
                    default='result.jpg',
                    help="distribution figures")
parser.add_argument('--use_dude_ism', action='store_true')
args = parser.parse_args()
# args = parser.parse_args([
#     "-f",
#     "/pubhome/jcyang/git/can-ai-do/dude/2split/diverse1.json",
#     # "/pubhome/jcyang/git/can-ai-do/dude/2split/crossFamilySplit/family3fold.json",
#     "-d",
#     # "/pubhome/jcyang/git/can-ai-do/dude/figures/diverse",
#     "/pubhome/jcyang/git/can-ai-do/dude/figures/diverse.heavyMW500",
#     # "/pubhome/jcyang/tmp/dude/all",
#     "-o",
#     "/pubhome/jcyang/git/can-ai-do/dude/figures/result"
# ])

#%%


def mfp2(m):
    # radius 2 MorganFingerprint equal ECFP4
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
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


def load_smiles(names):
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
            if args.use_dude_ism:
                activeFile = tdir / 'actives_final.ism'
                active_supp = Chem.SmilesMolSupplier(str(activeFile),titleLine=False)
            else:
                activeFile = tdir / 'actives_final.sdf.gz'
                active_supp = Chem.ForwardSDMolSupplier(gzip.open(activeFile))

        decoyFile = tdir / 'decoys_final.smi'
        if decoyFile.exists():
            # generate in this work
            decoy_supp = Chem.SmilesMolSupplier(str(decoyFile),
                                                titleLine=False)
        else:
            # from DUD-E
            if args.use_dude_ism:
                decoyFile = tdir / 'decoys_final.ism'
                decoy_supp = Chem.SmilesMolSupplier(str(decoyFile),titleLine=False)
            else:
                decoyFile = tdir / 'decoys_final.sdf.gz'
                decoy_supp = Chem.ForwardSDMolSupplier(gzip.open(decoyFile))

        propf = activeFile.with_name(activeFile.name + '.prop.MWHA.pkl')
        labelf = activeFile.with_name(activeFile.name + '.labelf.MWHA.pkl')
        if propf.exists() and labelf.exists():
            with open(propf, 'rb') as f:
                props = pickle.load(f)
            with open(labelf, 'rb') as f:
                labels = pickle.load(f)
        else:
            props = []
            labels = []
            for m in active_supp:
                if m is None: continue
                props.append(getProp(m))
                labels.append(1)
            for m in decoy_supp:
                if m is None: continue
                props.append(getProp(m))
                labels.append(0)
            with open(propf, 'wb') as f:
                pickle.dump(props, f)
            with open(labelf, 'wb') as f:
                pickle.dump(labels, f)
        all_props.extend(props)
        all_labels.extend(labels)
    return all_props, all_labels


with open(args.fold_list) as f:
    folds = json.load(f)
    if type(folds) is list:
        folds = {'{}'.format(fold): fold for fold in folds}
    targets = [i for fold in folds.values() for i in fold]

p = mp.Pool()
iter_targets = [[i] for i in targets]
for _ in tqdm(p.imap_unordered(load_smiles, iter_targets),
              desc='Converting smiles into fingerprints and properties',
              total=len(targets)):
    pass
p.close()

#%%
props, labels = load_smiles(targets)
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

fig.savefig(args.output, dpi=300)
print(f"result figure saved at {args.output}")
#%%

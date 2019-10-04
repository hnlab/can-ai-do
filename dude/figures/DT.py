"""A fingerprint + random forest model.
Try to generate independent and identically distributed figerprint as decoy.
"""
#%%
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
import multiprocessing as mp

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.DataStructs import BulkTanimotoSimilarity

from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier

from sklearn import tree

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
parser.add_argument(
    '--use_MW',
    action='store_false',
    help="use MolWt for random forset, default is HeavyAtomMolWt.")
parser.add_argument('--removeHeavyMW500',
                    action='store_true',
                    help="remove actives with HeavyAtomMolWt > 500.")
parser.add_argument(
    '--subsample',
    action='store_true',
    help=
    "randomly remove same number of actives and decoys as --removeHeavyMW500")
parser.add_argument(
    '-o',
    '--output',
    default='result',
    help=("prefix of output. default is 'result'." +
          "will output 2 files, one is result.performance.json," +
          "other is result.importance_features.json."))
args = parser.parse_args()


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
            if args.use_dude_ism:
                activeFile = tdir / 'actives_final.ism'
                active_supp = Chem.SmilesMolSupplier(str(activeFile),
                                                     titleLine=False)
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
                decoy_supp = Chem.SmilesMolSupplier(str(decoyFile),
                                                    titleLine=False)
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
        if removeHeavyMW500:
            props = np.array(props)
            labels = np.array(labels)
            actives = props[labels == 1]
            decoys = props[labels == 0]
            big_active_mask = actives[:, 0] > 500
            nbig = sum(big_active_mask)
            fraction = nbig / len(actives)
            # print(f"left {len(actives) - nbig:4d}/{len(actives):4d} from target {name}")
            if subsample:
                perm = np.random.permutation(len(actives))
                actives = actives[perm[nbig:]]
            else:
                actives = actives[~big_active_mask]
            perm = np.random.permutation(len(decoys))
            rmN = int(fraction * len(decoys))
            decoys = decoys[perm[rmN:]]
            props = np.vstack((actives, decoys))
            labels = np.hstack((np.ones(len(actives)), np.zeros(len(decoys))))
        all_props.extend(props)
        all_labels.extend(labels)
    # [mwha, mw, logp, rotb, hbd, hba, q]
    all_props = np.array(all_props)
    if HeavyAtomMolWt:
        all_props = all_props[:, (0, 2, 3, 4, 5, 6)]
    else:
        all_props = all_props[:, (1, 2, 3, 4, 5, 6)]
    return all_props, all_labels


def enrichment_factor(y_true, y_pred, first=0.01):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    n = len(y_pred)
    first_n = max(int(first * n), 1)
    indices = np.argsort(y_pred)[-first_n:]
    first_active_percent = np.sum(y_true[indices] == 1,
                                  dtype=np.float) / first_n
    active_percent = np.sum(y_true == 1, dtype=np.float) / n
    return first_active_percent / active_percent


def random_forest(fold):
    train_props, train_labels = load_smiles(
        fold,
        HeavyAtomMolWt=args.use_MW,
        removeHeavyMW500=args.removeHeavyMW500,
        subsample=args.subsample,
    )
    # XY = {'fp': (train_fps, train_labels), 'prop': (train_props, train_labels)}
    XY = {'prop': (train_props, train_labels)}
    results = {}
    for key, (X, Y) in XY.items():
        # clf = RandomForestClassifier(
        # n_estimators=1,
        # max_depth=10,
        # min_samples_split=10,
        # min_samples_split=5,
        # min_samples_leaf=2,
        # class_weight='balanced',
        # class_weight = {0: 1, 1: 50},
        # random_state=0,
        # n_jobs=8,
        # )
        clf = tree.DecisionTreeClassifier(
            max_depth=6,
            min_samples_split=5,
            #class_weight = {0: 1, 1: 50},
        )
        clf = clf.fit(X, Y)
        pred = clf.predict_proba(X)
        ROC = metrics.roc_auc_score(Y, pred[:, 1])
        EF1 = enrichment_factor(Y, pred[:, 1], first=0.01)
        print(fold, ROC, EF1)
        import graphviz 
        dot_data = tree.export_graphviz(
            clf,
            out_file=None,
            feature_names=['mwha', 'logp', 'rotb', 'hbd', 'hba', 'q'],
            class_names=['decoy', 'active'],
            # proportion=True,
            filled=True,
            rounded=True,
            special_characters=True,
        )
        graph = graphviz.Source(dot_data, format='png')
        graph.render(f"/tmp/tree_{fold[0]}")
    return results


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

p = mp.Pool()
result = p.map(random_forest, folds.values())
p.close()
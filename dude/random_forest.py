"""A fingerprint + random forest model.
Try to generate independent and identically distributed figerprint as decoy.
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
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    q = Chem.GetFormalCharge(mol)
    return tuple([mw, logp, rotb, hbd, hba, q])


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

        fpf = activeFile.with_name(activeFile.name + '.fp.pkl')
        propf = activeFile.with_name(activeFile.name + '.prop.pkl')
        labelf = activeFile.with_name(activeFile.name + '.labelf.pkl')
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
            for m in active_supp:
                if m is None: continue
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
        all_fps.extend(fps)
        all_props.extend(props)
        all_labels.extend(labels)
    return all_fps, all_props, all_labels


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


def random_forest(train_test):
    train_names, test_names = train_test
    train_fps, train_props, train_labels = load_smiles(train_names)
    XY = {'fp': (train_fps, train_labels), 'prop': (train_props, train_labels)}
    result = {}
    for key, (X, Y) in XY.items():
        clf = RandomForestClassifier(
            n_estimators=32,
            max_depth=10,
            # min_samples_split=10,
            min_samples_split=5,
            min_samples_leaf=2,
            # class_weight='balanced',
            class_weight={
                0: 1,
                1: 50
            },
            random_state=0,
            # n_jobs=8,
        )
        clf.fit(X, Y)
        result[key] = {'ROC': 0, 'EF1': 0}
        for test_name in test_names:
            test_fps, test_props, test_labels = load_smiles([test_name])
            if key == 'fp':
                test_X = test_fps
            if key == 'prop':
                test_X = test_props
            pred = clf.predict_proba(test_X)
            ROC = metrics.roc_auc_score(test_labels, pred[:, 1])
            # auc_prc = metrics.average_precision_score(test_labels, pred[:, 1])
            EF1 = enrichment_factor(test_labels, pred[:, 1], first=0.01)
            result[key]['ROC'] += ROC / len(test_names)
            result[key]['EF1'] += EF1 / len(test_names)
    return result


def most_simi(train_test):
    train_names, test_names = train_test
    train_fps, train_props, train_labels = load_smiles(train_names)
    test_fps, test_props, test_labels = load_smiles(test_names)
    train_actives = [fp for fp, y in zip(train_fps, train_labels) if y == 1]
    train_decoys = [fp for fp, y in zip(train_fps, train_labels) if y == 0]
    test_actives = [fp for fp, y in zip(test_fps, test_labels) if y == 1]
    test_decoys = [fp for fp, y in zip(test_fps, test_labels) if y == 0]
    fps_pairs = {
        'test_actives vs test_decoys': (test_actives, test_decoys),
        'test_actives vs train_actives': (test_actives, train_actives),
        'test_decoys vs train_actives': (test_decoys, train_actives),
        'test_actives vs train_decoys': (test_actives, train_decoys),
        'test_decoys vs train_decoys': (test_decoys, train_decoys),
    }
    most_simi = {}
    for k, (fps_a, fps_b) in fps_pairs.items():
        most_simi[k] = [max(BulkTanimotoSimilarity(i, fps_b)) for i in fps_a]
    return most_simi


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

train_test_pairs = []
fold_names = []
for k, fold in folds.items():
    fold_names.append(k)
    test_names = fold
    train_names = [name for ki, vi in folds.items() for name in vi if ki != k]
    train_test_pairs.append((train_names, test_names))

nfold = len(train_test_pairs)

p = mp.Pool(min(nfold, mp.cpu_count()))
iter_result = tqdm(p.imap(random_forest, train_test_pairs),
                   desc='Benchmarking random forest model on each fold',
                   total=nfold)
performance_on_fold = [i for i in iter_result]
p.close()

output = Path(args.output)
# add .suffix in with_suffix() for output with dot '.'
with open(output.with_suffix(output.suffix + '.json'), 'w') as f:
    result = []
    mean = {}
    for name, performance in zip(fold_names, performance_on_fold):
        for key in performance:
            if key not in mean:
                mean[key] = {}
            for k, v in performance[key].items():
                if k not in mean[key]:
                    mean[key][k] = 0
                mean[key][k] += performance[key][k] / nfold
        performance = performance.copy()
        performance['fold'] = name
        result.append(performance)
    result.append(mean)
    json.dump(result, f, sort_keys=True, indent=4)
    print('save performance at {}'.format(f.name))

"""A fingerprint + random forest model.
Try to generate independent and identically distributed figerprint as decoy.
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
    '--use_MW',
    action='store_false',
    help="use MolWt for random forset, default is HeavyAtomMolWt.")
parser.add_argument('--random_fold',
                    action='store_true',
                    help="use random folds")
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
    if HeavyAtomMolWt:
        all_props = all_props[:, (0, 2, 3, 4, 5, 6)]
    else:
        all_props = all_props[:, (1, 2, 3, 4, 5, 6)]
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
    train_fps, train_props, train_labels = load_smiles(
        train_names,
        HeavyAtomMolWt=args.use_MW,
        removeHeavyMW500=args.removeHeavyMW500,
        subsample=args.subsample,
    )
    XY = {'fp': (train_fps, train_labels), 'prop': (train_props, train_labels)}
    # XY = {'prop': (train_props, train_labels)}
    results = {}
    for key, (X, Y) in XY.items():
        result = {'ROC': {}, 'EF1': {}}
        clf = RandomForestClassifier(
            n_estimators=100,
            # max_depth=10,
            # min_samples_split=10,
            # min_samples_split=5,
            # min_samples_leaf=2,
            # class_weight='balanced',
            # class_weight={
            #     0: 1,
            #     1: 50
            # },
            random_state=0,
            # n_jobs=8,
        )
        clf.fit(X, Y)
        for test_name in test_names:
            test_fps, test_props, test_labels = load_smiles(
                [test_name],
                HeavyAtomMolWt=args.use_MW,
                removeHeavyMW500=args.removeHeavyMW500,
                subsample=args.subsample,
            )
            if key == 'fp':
                test_X = test_fps
            if key == 'prop':
                test_X = test_props
            pred = clf.predict_proba(test_X)
            ROC = metrics.roc_auc_score(test_labels, pred[:, 1])
            # auc_prc = metrics.average_precision_score(test_labels, pred[:, 1])
            EF1 = enrichment_factor(test_labels, pred[:, 1], first=0.01)
            result['ROC'][test_name] = ROC
            result['EF1'][test_name] = EF1
        result['feature_importances'] = list(clf.feature_importances_)
        results[key] = result
    return results


with open(args.fold_list) as f:
    folds = json.load(f)
    if type(folds) is list:
        folds = {'{}'.format(fold): fold for fold in folds}
    targets = [i for fold in folds.values() for i in fold]

iter_targets = [[i] for i in targets]
active_num = 0
decoy_num = 0
p = mp.Pool()
for _, _, labels in tqdm(
        p.imap_unordered(load_smiles, iter_targets),
        desc='Converting smiles into fingerprints and properties',
        total=len(targets)):
    labels = np.asarray(labels)
    active_num += sum(labels == 1)
    decoy_num += sum(labels == 0)
p.close()
print(f"total actives:{active_num}, total decoy:{decoy_num}")

feature_sets = ('prop', 'fp')
metric_names = ('ROC', 'EF1')
np.random.seed(123)
repeat = 1
repeat_results = []
repeat_means = []
for r in range(repeat):
    tmp_folds = folds
    if args.random_fold:
        perm = np.random.permutation(len(targets))
        targets = np.array(targets)
        tmp_folds = {}
        start = 0
        for k, fold in folds.items():
            end = start + len(fold)
            tmp_folds[k] = list(targets[perm[start:end]])
            start = end

    train_test_pairs = []
    fold_names = []
    for k, fold in tmp_folds.items():
        fold_names.append(k)
        test_names = fold
        train_names = [
            name for ki, vi in tmp_folds.items() for name in vi if ki != k
        ]
        train_test_pairs.append((train_names, test_names))

    nfold = len(train_test_pairs)

    p = mp.Pool(min(nfold, mp.cpu_count()))
    iter_result = tqdm(p.imap(random_forest, train_test_pairs),
                       desc='Benchmarking random forest model on each fold',
                       total=nfold)
    performance_on_fold = [i for i in iter_result]
    p.close()

    result = {}
    for name, performance in zip(fold_names, performance_on_fold):
        for feat in feature_sets:
            if feat not in result:
                result[feat] = {}
            print(performance.keys())
            print(performance[feat].keys())
            result[feat]['feature_importances'] = performance[feat]['feature_importances']
            for metric in metric_names:
                if metric not in result[feat]:
                    result[feat][metric] = {}
                result[feat][metric].update(performance[feat][metric])
    mean = {}
    for feat in feature_sets:
        mean[feat] = {}
        for metric in metric_names:
            mean[feat][metric] = np.mean(list(result[feat][metric].values()))
    result['folds'] = tmp_folds
    result['mean'] = mean
    repeat_results.append(result)
    repeat_means.append(mean)
    print(mean)

target_performances = []
for result in repeat_results:
    for feat in feature_sets:
        if feat not in result:
            continue
        for metric in metric_names:
            for target, value in result[feat][metric].items():
                target_performances.append((target, feat, metric, value))
df = pd.DataFrame(data=target_performances,
                  columns=['target', 'feat', 'metric', 'value'])
output = Path(args.output)
# add .suffix in with_suffix() for output with dot '.'
with open(output.with_suffix(output.suffix + '.json'), 'w') as f:
    final_result = [repeat_results, repeat_means]
    json.dump(final_result, f, sort_keys=True, indent=4)
    print(f'save performance at {f.name}')
csv = output.with_suffix(output.suffix + '.csv')
df.to_csv(csv, index=False)
print(f'save target performance at {csv}')
for feat in feature_sets:
    EF1 = df[(df.feat == feat) & (df.metric == 'EF1')]
    if EF1.empty:
        continue
    grouped = EF1.groupby(['target', 'feat', 'metric']).mean().reset_index()
    sorted_ = grouped.sort_values(by=['value'], ascending=False)
    sorted_EF1 = EF1.set_index('target').loc[sorted_.target]
    sorted_csv = output.with_suffix(output.suffix + f'.{feat}.sorted.csv')
    sorted_EF1.to_csv(sorted_csv)
    print(f'save sorted target performance at {sorted_csv}')

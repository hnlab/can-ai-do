"""Random forest on DUD-E
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
parser.add_argument('--use_dude_sdf', action='store_true')
parser.add_argument(
    '--use_MW',
    action='store_false',
    help="use MolWt for random forset, default is HeavyAtomMolWt.")
parser.add_argument('--random_fold',
                    action='store_true',
                    help="use random folds")
parser.add_argument('--MW500',
                    action='store_true',
                    help="remove actives with HeavyAtomMolWt > 500.")
parser.add_argument('--bits', help='only using FP bits in the json file.')
parser.add_argument(
    '-o',
    '--output',
    default='result',
    help=("prefix of output. default is 'result'." +
          "will output 2 files, one is result.performance.json," +
          "other is result.importance_features.json."))
args = parser.parse_args()

if args.use_MW:
    args.MolWt = 'MolWt'
else:
    args.MolWt = 'HeavyAtomMolWt'

if args.bits:
    with open(args.bits) as f:
        args.bits = json.load(f)
    print(f"only using {len(args.bits)} FP bits {args.bits}")

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


def load_dude(names, MolWt=None, MW500=False, fpAsArray=False, bits=None):
    if bits is not None:
        fpAsArray = True
    datadir = Path(args.datadir)
    all_ids = []
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

        idf = activeFile.with_name(activeFile.name + '.id.pkl')
        fpf = activeFile.with_name(activeFile.name + '.fp.pkl')
        propf = activeFile.with_name(activeFile.name + '.prop.pkl')
        labelf = activeFile.with_name(activeFile.name + '.label.pkl')
        idf_mw500 = idf.with_suffix('.MW500.pkl')
        fpf_mw500 = fpf.with_suffix('.MW500.pkl')
        propf_mw500 = propf.with_suffix('.MW500.pkl')
        labelf_mw500 = labelf.with_suffix('.MW500.pkl')

        if not all([f.exists() for f in (idf, fpf, propf, labelf)]):
            ids = []
            fps = []
            props = []
            labels = []
            ids_mw500 = []
            fps_mw500 = []
            props_mw500 = []
            labels_mw500 = []
            for m in active_supp:
                if m is None: continue
                mol_id = f'{name}_{m.GetProp("_Name")}'
                ids.append(mol_id)
                fp = mfp2(m)
                fps.append(fp)
                p = getProp(m)
                props.append(p)
                labels.append(1)
                # p:[mwha, mw, logp, rotb, hbd, hba, q]
                if p[0] > 500:
                    continue
                ids_mw500.append(mol_id)
                fps_mw500.append(fp)
                props_mw500.append(p)
                labels_mw500.append(1)
            frac = len(fps_mw500) / len(fps)
            decoy_mols = [m for m in decoy_supp if m is not None]
            select_num = int(frac * len(decoy_mols))
            np.random.seed(123)
            inds = np.random.choice(len(decoy_mols), select_num, replace=False)
            for i, m in enumerate(decoy_mols):
                mol_id = f'{name}_{m.GetProp("_Name")}'
                ids.append(mol_id)
                fp = mfp2(m)
                fps.append(fp)
                p = getProp(m)
                props.append(p)
                labels.append(0)
                if i in inds:
                    ids_mw500.append(mol_id)
                    fps_mw500.append(fp)
                    props_mw500.append(p)
                    labels_mw500.append(0)

            with open(idf, 'wb') as f:
                pickle.dump(ids, f)
            with open(fpf, 'wb') as f:
                pickle.dump(fps, f)
            with open(propf, 'wb') as f:
                pickle.dump(props, f)
            with open(labelf, 'wb') as f:
                pickle.dump(labels, f)

            with open(idf_mw500, 'wb') as f:
                pickle.dump(ids_mw500, f)
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
            idf = idf_mw500
            fpf = fpf_mw500
            propf = propf_mw500
            labelf = labelf_mw500

        if fpAsArray:
            fpf_np = fpf.with_suffix('.np.pkl')
            with open(fpf_np, 'rb') as f:
                fps = pickle.load(f)
        else:
            with open(fpf, 'rb') as f:
                fps = pickle.load(f)

        if bits is not None:
            fps = fps[:, bits]

        with open(idf, 'rb') as f:
            ids = pickle.load(f)
        with open(propf, 'rb') as f:
            props = pickle.load(f)
        with open(labelf, 'rb') as f:
            labels = pickle.load(f)

        all_ids.extend(ids)
        all_props.extend(props)
        all_labels.extend(labels)

        all_fps.append(fps)

    if fpAsArray:
        all_fps = np.vstack(all_fps)
    else:
        all_fps = sum(all_fps, [])  # flatten list of list

    # prop: [mwha, mw, logp, rotb, hbd, hba, q]
    all_props = np.array(all_props)
    if MolWt == 'HeavyAtomMolWt':
        all_props = all_props[:, (0, 2, 3, 4, 5, 6)]
    if MolWt == 'MolWt':
        all_props = all_props[:, (1, 2, 3, 4, 5, 6)]
    return all_ids, all_fps, all_props, all_labels


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
    train_ids, train_fps, train_props, train_labels = load_dude(
        train_names,
        MolWt=args.MolWt,
        fpAsArray=True,
        bits=args.bits,
        MW500=args.MW500)
    XY = {'fp': (train_fps, train_labels), 'prop': (train_props, train_labels)}
    # XY = {'prop': (train_props, train_labels)}
    results = {}
    for key, (X, Y) in XY.items():
        result = {'ROC': {}, 'EF1': {}}
        clf = RandomForestClassifier(n_estimators=100, random_state=0)
        clf.fit(X, Y)
        for test_name in test_names:
            test_ids, test_fps, test_props, test_labels = load_dude(
                [test_name],
                MolWt=args.MolWt,
                fpAsArray=True,
                bits=args.bits,
                MW500=args.MW500)
            test_ids = np.asarray(test_ids)
            test_labels = np.asarray(test_labels)

            if key == 'fp':
                test_X = test_fps
            if key == 'prop':
                test_X = test_props

            pred = clf.predict_proba(test_X)
            y_pred = pred[:, 1]

            sort_indices = np.argsort(-y_pred)
            test_ids = test_ids[sort_indices]
            uniq_ids, uniq_indices = np.unique(test_ids, return_index=True)
            y_pred = y_pred[sort_indices][uniq_indices]
            test_labels = test_labels[sort_indices][uniq_indices]

            ROC = metrics.roc_auc_score(test_labels, y_pred)
            EF1 = enrichment_factor(test_labels, y_pred, first=0.01)
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
p = mp.Pool()
for _ in tqdm(p.imap_unordered(load_dude, iter_targets),
              desc='Converting smiles into fingerprints and properties',
              total=len(targets)):
    pass
p.close()

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
            feat_imports = performance[feat]['feature_importances']
            if 'feature_importances' in result[feat]:
                result[feat]['feature_importances'].append(feat_imports)
            else:
                result[feat]['feature_importances'] = [feat_imports]
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

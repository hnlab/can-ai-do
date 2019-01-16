"""A fingerprint + random forest model.
Show feature importance
"""
import os
import sys
import json
import argparse
import numpy as np
import scipy.sparse as sp
from multiprocessing import Pool

from rdkit import Chem
from rdkit.Chem import AllChem

from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f', '--fold_list', required=True,
    help="k-fold config in json, a dict or a list of list")
parser.add_argument('-d', '--datadir', default='./all',
    help="datadir, default is ./all")
parser.add_argument('-m', '--mask_important_features', type=int, default=0,
    help="mask first n important features")
parser.add_argument('-o', '--output', default='result',
    help=("prefix of output. default is 'result'."
         +"will output 2 files, one is result.performance.json,"
         +"other is result.importance_features.json."))
args = parser.parse_args()

def mol_to_fp(m):
  # radius 2 MorganFingerprint equal ECFP4
  fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
  return fp

def load_smiles(names, datadir='./'):
  all_fps = []
  all_labels = []
  ws = []
  print("Loading data set ...")
  for name in names:
    file_name = os.path.join(datadir, name)
    activeFile = os.path.join(file_name, 'actives_final.ism')
    decoyFile = os.path.join(file_name, 'decoys_final.ism')
    fpf = os.path.join(file_name, 'fp.npz')
    labelf = os.path.join(file_name, 'label.npz')
    if os.path.exists(fpf):
      # print("Reloading data set {} ...".format(fpf))
      fps = sp.load_npz(fpf)
      labels = sp.load_npz(labelf)
    else:
      # print("Loading smiles {} ...".format(activeFile))
      active = [m for m in Chem.SmilesMolSupplier(activeFile, titleLine=False) if m is not None]
      # print("Loading smiles {} ...".format(decoyFile))
      decoy = [m for m in Chem.SmilesMolSupplier(decoyFile, titleLine=False) if m is not None]
      labels = sp.csr_matrix(np.hstack((np.ones(len(active)), np.zeros(len(decoy)))))
      ms = np.hstack((active, decoy))
      # print("Computing figerprints ...")
      with Pool() as p:
        fps = sp.csr_matrix(p.map(mol_to_fp, ms))
      #names = np.array([m.GetProp('_Name') for m in ms])
      sp.save_npz(fpf, fps)
      sp.save_npz(labelf, labels)
    all_fps.append(fps)
    all_labels.append(labels)
  return sp.vstack(all_fps), sp.hstack(all_labels).toarray().flatten()

with open(args.fold_list) as f:
  folds = json.load(f)
  if type(folds) is list:
    folds = {i:fold for i,fold in enumerate(folds)}
def enrichment_factor(y_true, y_pred, first=0.01):
  y_true = np.asarray(y_true)
  y_pred = np.asarray(y_pred)
  n = len(y_pred)
  first_n = first*n
  indices = np.argsort(y_pred)[-int(first_n):]
  first_active_percent = np.sum(y_true[indices]==1, dtype=np.float)/first_n
  active_percent = np.sum(y_true==1, dtype=np.float)/n
  return first_active_percent/active_percent

result = []
result_features = []
for k, fold in folds.items():
  test_names = fold
  train_names = [name for i, fold_ in folds.items() for name in fold_ if i != k]
  train_fps, train_labels = load_smiles(train_names, datadir=args.datadir)

  clf = RandomForestClassifier(
      n_estimators=32,
      max_depth=10,
      # min_samples_split=10,
      min_samples_split=5,
      min_samples_leaf=2,
      # class_weight='balanced',
      class_weight={0:1,1:50},
      random_state=0,
      n_jobs=8,
      )
  print("fiting model for {}".format(fold))
  clf = clf.fit(train_fps, train_labels)
  if args.mask_important_features:
    mask_indices = np.argsort(clf.feature_importances_)[-args.mask_important_features:]
    train_fps[:, mask_indices] = 0
    # retrain model
    clf = clf.fit(train_fps, train_labels)
  important_indices = np.argsort(clf.feature_importances_)[-10:]
  important_w = clf.feature_importances_[important_indices]
  important_fp = train_fps[:,important_indices]
  active_important_fp = important_fp[train_labels==1]
  decoy_important_fp = important_fp[train_labels==0]
  active_important_occ = active_important_fp.mean(axis=0)
  active_important_occ = tuple(np.squeeze(np.asarray(active_important_occ)))
  decoy_important_occ = decoy_important_fp.mean(axis=0).flatten()
  decoy_important_occ = tuple(np.squeeze(np.asarray(decoy_important_occ)))
  fold_result_features = [
    important_indices,
    important_w,
    active_important_occ,
    decoy_important_occ,
    ]
  result_features.append(fold_result_features)
  print("first 10 important features explain {}%.".format(100*sum(important_w)))
  print("important_fp, importantance, active_occ, decoy_occ")
  for idx, w, a_occ, d_occ in zip(*fold_result_features):
    print("{:4d}, {:.3f}, {:.3f}, {:.3f}".format(idx, w, a_occ, d_occ))
  for test_name in test_names:
    test_fps, test_labels = load_smiles([test_name], datadir=args.datadir)
    if args.mask_important_features:
      test_fps[:, mask_indices] = 0
    pred = clf.predict_proba(test_fps)
    auc_roc = metrics.roc_auc_score(test_labels, pred[:,1])
    print('{} AUC_ROC: {}'.format(test_name, auc_roc))
    auc_prc = metrics.average_precision_score(test_labels, pred[:,1])
    print('{} AUC_PRC: {}'.format(test_name, auc_prc))
    enrich_factor0_5 = enrichment_factor(test_labels, pred[:,1], first=0.005)
    print('{} enrich_factor_0.5%: {}'.format(test_name, enrich_factor0_5))
    enrich_factor1 = enrichment_factor(test_labels, pred[:,1], first=0.01)
    print('{} enrich_factor_1%: {}'.format(test_name, enrich_factor1))
    result.append([test_name, auc_roc, auc_prc, enrich_factor0_5, enrich_factor1])

avg_roc = sum([i[1] for i in result])/len(result)
avg_prc = sum([i[2] for i in result])/len(result)
avg_ef0_5 = sum([i[3] for i in result])/len(result)
avg_ef1 = sum([i[4] for i in result])/len(result)
result.append(['average', avg_roc, avg_prc, avg_ef0_5, avg_ef1])
print("average:\n"
  + "roc: {}\n".format(avg_roc)
  + "prc: {}\n".format(avg_prc)
  + "EF0.5%: {}\n".format(avg_ef0_5)
  + "EF1%: {}\n".format(avg_ef1)
  )
with open(args.output+'.performance.json', 'w') as f:
  json.dump(result, f, indent=2)

"""A fingerprint + random forest model.
Try to generate independent and identically distributed figerprint as decoy.
"""
import os
import sys
import json
import argparse
import numpy as np
import scipy.sparse as sp
from scipy.spatial import distance
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
parser.add_argument('-I', '--indepen_indent_distr', action="store_true",
    help="use independent and identically distributed fingerprint as decoy.")
parser.add_argument('-R', '--random_zinc_decoy', action="store_true",
    help="use random pick drug like mol. from zinc as decoys.")
parser.add_argument('-X', '--dudx', action="store_true",
    help="use reducing similar DUDX.")
parser.add_argument('-o', '--output', default='result',
    help=("prefix of output. default is 'result'."
         +"will output 2 files, one is result.performance.json,"
         +"other is result.importance_features.json."))
args = parser.parse_args()

def mol_to_fp(m):
  # radius 2 MorganFingerprint equal ECFP4
  fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
  return fp

# import numba
# @numba.njit()
# def jaccard(a,b):
#   ab = a.dot(b.T)[0,0]
#   return  1.0 - ab / (a.sum()+b.sum()-ab)
#
# def indepen_indent_distr(fps):
#   # fps = active_fps.toarray()
#   fps = fps[:10]
#   n_samples, n_bits = fps.shape
#   p = fps.mean(axis=0)
#   # matrix to 1d array
#   p = np.squeeze(np.asarray(p))
#   p[p>1] = 1.
#   print("generate decoy fps")
#   iid_fps = np.random.binomial(1, p, (n_samples*100, n_bits))
#   print("compute dist")
#   # dist = distance.cdist(iid_fps, fps, 'jaccard')
#   iid_fps = sp.csr_matrix(iid_fps)
#   dist = metrics.pairwise.pairwise_distances(
#       X=fps, Y=iid_fps, metric=jaccard, n_jobs=-1)
#   print(np.sum(dist<0.8)/(100*n_samples**2))
#   return None

def indepen_indent_distr(fps, w=50):
  fps = fps.toarray()
  fps = fps.astype(np.int8)
  n_samples, n_bits = fps.shape
  p = fps.mean(axis=0)
  # matrix to 1d array
  p = np.squeeze(np.asarray(p))
  p[p>1] = 1.
  print("generate decoy fps")
  iid_fps = np.random.binomial(1, p, (n_samples*w, n_bits))
  iid_fps = iid_fps.astype(np.int8)
  print("compute dist")
  # dist = distance.cdist(fps, iid_fps, 'jaccard') # slower
  dist = metrics.pairwise.pairwise_distances(
    X=fps, Y=iid_fps, metric='jaccard', n_jobs=4)
  print('{}% decoys have similarity > 0.5'.format(np.sum(dist<0.5)/(w*n_samples**2)))
  X = sp.vstack((
    sp.csr_matrix(fps),
    sp.csr_matrix(iid_fps)
    ))
  y = np.hstack((
    np.ones(n_samples, dtype=np.int8),
    np.zeros(n_samples*w, dtype=np.int8)
    ))
  return X, y

def load_smiles(names, datadir='./', iid=False, random=False, dudx=False):
  all_fps = []
  all_labels = []
  ws = []
  print("Loading data set ...")
  for name in names:
    file_name = os.path.join(datadir, name)
    activeFile = os.path.join(file_name, 'actives_final.ism')
    if random:
      fpf = os.path.join(file_name, 'fpR.npz')
      labelf = os.path.join(file_name, 'labelR.npz')
      decoyFile = os.path.join(file_name, 'decoys_random.smi')
    elif iid:
      fpf = os.path.join(file_name, 'fpI.npz')
      labelf = os.path.join(file_name, 'labelI.npz')
    elif dudx:
      fpf = os.path.join(file_name, 'fpX.npz')
      labelf = os.path.join(file_name, 'labelX.npz')
    else:
      fpf = os.path.join(file_name, 'fp.npz')
      labelf = os.path.join(file_name, 'label.npz')
      decoyFile = os.path.join(file_name, 'decoys_final.ism')
    if os.path.exists(fpf) and os.path.exists(labelf):
      # print("Reloading data set {} ...".format(fpf))
      fps = sp.load_npz(fpf)
      labels = sp.load_npz(labelf)
    else:
      # print("Loading smiles {} ...".format(activeFile))
      active = [m for m in Chem.SmilesMolSupplier(activeFile, titleLine=False) if m is not None]
      # print("Loading smiles {} ...".format(decoyFile))
      with Pool() as p:
        active_fps = sp.csr_matrix(p.map(mol_to_fp, active))
      if iid:
        fps, labels = indepen_indent_distr(active_fps)
        labels = sp.csr_matrix(labels)
      else:
        decoy = [m for m in Chem.SmilesMolSupplier(decoyFile, titleLine=False) if m is not None]
        with Pool() as p:
          decoy_fps = sp.csr_matrix(p.map(mol_to_fp, decoy))
        fps = sp.vstack((active_fps, decoy_fps))
        labels = sp.csr_matrix(np.hstack((np.ones(len(active)), np.zeros(len(decoy)))))
      # print("Computing figerprints ...")
      #names = np.array([m.GetProp('_Name') for m in ms])
      sp.save_npz(fpf, fps)
      sp.save_npz(labelf, labels)
    labels = np.squeeze(labels.toarray())
    all_fps.append(fps)
    all_labels.append(labels)
    # X = sp.vstack(all_fps)
    # y =  sp.hstack(all_labels).toarray().flatten()
    X = sp.vstack(all_fps)
    y = np.hstack(all_labels)
  return X, y

with open(args.fold_list) as f:
  folds = json.load(f)
  if type(folds) is list:
    folds = {i:fold for i,fold in enumerate(folds)}
def enrichment_factor(y_true, y_pred, first=0.01):
  y_true = np.asarray(y_true)
  y_pred = np.asarray(y_pred)
  n = len(y_pred)
  first_n = max(int(first*n), 1)
  indices = np.argsort(y_pred)[-first_n:]
  first_active_percent = np.sum(y_true[indices]==1, dtype=np.float)/first_n
  active_percent = np.sum(y_true==1, dtype=np.float)/n
  return first_active_percent/active_percent

result = []
result_features = []
for k, fold in folds.items():
  test_names = fold
  train_names = [name for i, fold_ in folds.items() for name in fold_ if i != k]
  train_fps, train_labels = load_smiles(
      train_names, datadir=args.datadir,
      iid=args.indepen_indent_distr,
      random=args.random_zinc_decoy,
      dudx=args.dudx,
      )
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
  result_features.append(tuple(clf.feature_importances_))
  # result_features.append(fold_result_features)
  print("first 10 important features explain {}%.".format(100*sum(important_w)))
  print("important_fp, importantance, active_occ, decoy_occ")
  for idx, w, a_occ, d_occ in zip(*fold_result_features):
    print("{:4d}, {:.3f}, {:.3f}, {:.3f}".format(idx, w, a_occ, d_occ))
  for test_name in test_names:
    test_fps, test_labels = load_smiles(
        [test_name], datadir=args.datadir,
        iid=args.indepen_indent_distr,
        random=args.random_zinc_decoy,
        dudx=args.dudx,
        )
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
with open(args.output+'.feature_importances.json', 'w') as f:
  json.dump(result_features, f, indent=2)

"""A fingerprint + random forest model.
Try to generate independent and identically distributed figerprint as decoy.
"""
import os
import sys
import json
import argparse
import numpy as np
from pathlib import Path
from tqdm import tqdm

import scipy.sparse as sp
from scipy.spatial import distance
from multiprocessing import Pool

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Descriptors

from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', '--index', required=True)
parser.add_argument(
    '-d', '--datadir', required=True, help="pdbbind datadir, like v2018")
parser.add_argument(
        '-u', '--uclust', help="uclust output, format: https://www.drive5.com/usearch/manual/opt_uc.html")
args = parser.parse_args()

DATADIR = Path(args.datadir)

def read_index(index_file):
    codes = []
    pKs = []
    with open(index_file) as f:
        for i in f:
            if i[0] == '#': continue
            code, reso, year, pK, *others = i.split()
            codes.append(code)
            pKs.append(float(pK))
    return codes, pKs

def getProp(mol):
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    q = Chem.GetFormalCharge(mol)
    return tuple([mw, logp, rotb, hbd, hba, q])

def load_fps(codes):
    print("Loading ligand fingerprint")
    fps = []
    for i, code in tqdm(enumerate(codes), total=len(codes)):
        # already converted ligand.mol2 to ligand.pdb by babel
        path = DATADIR / code / (code + '_ligand.pdb')
        if not path.exists():
            fps.append(None)
            continue
        mol = Chem.MolFromPDBFile(str(path))
        if mol is None:
            fps.append(None)
            continue
        # fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)
        fp = getProp(mol)
        fps.append(fp)
    notNone = sum([1 for i in fps if i is not None])
    print('succeed loaded {}/{}'.format(notNone, len(codes)))
    return fps


def load_clust(uclust_file, codes):
    clust_nums = [None for i in codes]
    all_clust_nums = []
    labels = []
    with open(uclust_file) as f:
        for line in f:
            fields = line.split()
            all_clust_nums.append( int(fields[1]))
            labels.append(fields[8])
    for i, code in enumerate(codes):
        try:
            idx = labels.index(code)
            clust_nums[i] = all_clust_nums[idx]
        except ValueError:
            continue
    return clust_nums

codes, pKs = read_index(args.index)
fps = load_fps(codes)
Nones = [i for i in range(len(codes)) if fps[i] is None]
fps = [j for i,j in enumerate(fps) if i not in Nones]
pKs = [j for i,j in enumerate(pKs) if i not in Nones]
codes = [j for i,j in enumerate(codes) if i not in Nones]
X = np.array(fps)

if args.uclust:
    clust_nums = load_clust(args.uclust, codes)
    Nones.extend([i for i in range(len(codes)) if clust_nums[i] is None])
    Nones = set(Nones)
    fps = [j for i,j in enumerate(fps) if i not in Nones]
    pKs = [j for i,j in enumerate(pKs) if i not in Nones]
    codes = [j for i,j in enumerate(codes) if i not in Nones]
    clust_nums = [j for i,j in enumerate(clust_nums) if i not in Nones]
    clust_nums = np.array(clust_nums, dtype=int)
    join_clust = np.zeros_like(clust_nums)
    for i, num in enumerate(set(clust_nums)):
        mask = clust_nums == num
        # all cluster smaller than 5 will set as cluster 0
        if sum(mask) >= 10:
            join_clust[mask] = i+1
    nb_clust = max(join_clust) + 1
    print(join_clust)
    one_hot = np.eye(nb_clust, dtype=int)[join_clust]
    X = np.hstack((one_hot, fps))
    X = one_hot
    print(X.shape)


pKs = np.array(pKs)
# filter None
for seed in (111, 222, 333):
    np.random.seed(seed)
    N = len(codes)
    perm = np.random.permutation(N)
    train_idx = perm[:int(N*0.8)]
    valid_idx = perm[int(N*0.8):int(N*0.9)]
    test_idx = perm[int(N*0.9):]
    train_X = X[train_idx]
    test_X = X[test_idx]
    train_pKs = pKs[train_idx]
    test_pKs = pKs[test_idx]

    clf = RandomForestRegressor(
        n_estimators=10,
        max_depth=15,
        # min_samples_split=10,
        min_samples_split=5,
        min_samples_leaf=1,
        random_state=0,
        n_jobs=8,
    )
    clf.fit(train_X, train_pKs)
    pred_pKs = clf.predict(test_X)
    r2 = np.corrcoef(test_pKs, pred_pKs)[0,1] ** 2
    print('seed {} r2: {}'.format(seed, r2))

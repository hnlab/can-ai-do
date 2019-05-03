"""Load mols from smiles into database.
"""
import argparse
from pathlib import Path
from datetime import datetime as dt

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import dendrogram, linkage

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-s", "--smiles", required=True)
parser.add_argument("-o", "--output", default='fig.png')
args = parser.parse_args()

start = dt.now()

mols = [
    m for m in Chem.SmilesMolSupplier(args.smiles, titleLine=False)
    if m is not None
]

fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in mols]

dists = [] # condensed distance matrix
N = len(fps)
for i in range(N):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:])
    dists.extend([1 - x for x in sims])
# print(sorted(dists)[:3])
Z = linkage(dists, method='single')
# Z = linkage(dists, method='ward')
# Z = linkage(dists, method='average')
# Z = linkage(dists, method='centroid')
# print(Z)
fig = plt.figure(figsize=(25, 10))
# dn = dendrogram(Z, distance_sort=True, color_threshold=0.2)
# dn = dendrogram(Z, count_sort=True, color_threshold=0.2)
dn = dendrogram(Z, count_sort=True, color_threshold=0.4)
plt.ylim(0,1)
fig.savefig(args.output)
print('Fig. saved at {}'.format(args.output))
# print("Total elapsed time: {}".format(dt.now() - start))

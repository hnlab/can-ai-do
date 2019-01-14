"""load pdbbind into ani dataset in hdf5.
"""
import sys
import json
import argparse
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import dendrogram, linkage  

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', '--input', required=True,
    help="input blast result in json format")
parser.add_argument('-o', '--output', default='data',
    help="prefix of output. default is data")
args = parser.parse_args()

names = []
dist_matrix = []
with open(args.input) as f:
  blastResult = json.load(f)

for i, out in enumerate(blastResult["BlastOutput2"]):
  for j, result in enumerate(out['report']['results']['bl2seq']):
    # 'query_title': 'kinase27_abl1 abl1 270 bp; generated with OpenBabel 2.4.1'
    namei = result['query_title'].split()[0]
    query_len = result['query_len']
    if j <= i:
      continue
    # print(json.dumps(result, indent=2))
    if len(result['hits']) > 1: print("#"*100)
    if result['hits']:
      namej = result['hits'][0]["description"][0]["title"].split()[0]
      identityNum = result['hits'][0]['hsps'][0]["identity"]
      subject_len = result['hits'][0]['len']
      identity = float(identityNum)/min(query_len, subject_len)
      d = 1.0 - identity
      if i == j:
        assert namei == namej
    else:
      d = 0.85
    dist_matrix.append(d)
  names.append(namei)

# dist_matrix is condensed distance matrix
linked = linkage(dist_matrix, method='single')

plt.figure(figsize=(10, 7))
dendrogram(linked,  
  orientation='top',
  labels=names,
  distance_sort='descending',
  show_leaf_counts=True)
plt.savefig("blast_linkage.svg")

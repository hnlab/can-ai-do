import json
import argparse
import numpy as np
from tqdm import tqdm
from pathlib import Path
from datetime import datetime as dt
start = dt.now()

from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from Bio import PDB
from Bio.Seq import Seq
from Bio import pairwise2 as pw2

from multiprocessing import Pool

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', '--index', required=True)
parser.add_argument('-d',
                    '--datadir',
                    default='./v2018',
                    help="datadir, default is ./v2018")
args = parser.parse_args()

DATADIR = Path(args.datadir)


def SeqFromPDBCode(code):
    protein_pdb = DATADIR / code / (code + '_protein.pdb')
    pocket_pdb = DATADIR / code / (code + '_pocket.pdb')
    parser = PDB.PDBParser(QUIET=True)
    chain_id = None
    try:
        pocket = parser.get_structure(code, pocket_pdb)
        protein = parser.get_structure(code, protein_pdb)
    except:
        return None
    longest_chain = None
    for chain in pocket.get_chains():
        if chain.id == ' ': continue
        if longest_chain is None or len(chain) > len(longest_chain):
            longest_chain = chain
    if longest_chain is None:
        return None
    ppb = PDB.PPBuilder()
    for chain in protein.get_chains():
        if chain.id == longest_chain.id:
            seqs = [i.get_sequence() for i in ppb.build_peptides(chain)]
            seq_str = ''.join([str(i) for i in seqs])
            a = seqs[0].alphabet
            return Seq(seq_str, a)


def identity_percent(seqA, seqB):
    # Question: python: find percentage of match between two sequences
    # https://www.biostars.org/p/208540/
    # seq_length = min(len(seqA), len(seqB))
    # matches = pw2.align.globalxx(seqA, seqB, one_alignment_only=True, score_only=True)
    # global_align = pw2.align.globalxx(seqA, seqB)
    global_align = pw2.align.globalxx(seqA, seqB, one_alignment_only=True)
    seq_length = global_align[0][4]
    matches = global_align[0][2]
    return (matches / seq_length) * 100


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


codes, pKs = read_index(args.index)
seqs = []
Nones = []
p = Pool()
# iter_seqs = p.imap_unordered(SeqFromPDBCode, codes)
# iter_seqs = map(SeqFromPDBCode, codes)
iter_seqs = p.imap(SeqFromPDBCode, codes)
for i, seq in tqdm(enumerate(iter_seqs), total=len(codes)):
    if seq is None:
        Nones.append(i)
    else:
        seqs.append(seq)
p.close()
print('succeeded {}/{}\n'.format(len(seqs), len(codes)))
codes = [j for i, j in enumerate(codes) if i not in Nones]
pKs = [j for i, j in enumerate(pKs) if i not in Nones]


def identity_dist(i, j):
    seqA = seqs[int(i[0])]
    seqB = seqs[int(j[0])]
    return 100 - identity_percent(seqA, seqB)


seqs_2d = [[i] for i in range(len(seqs))]

# dm = pdist(seqs_2d, identity_dist)

print("Calculate distance matrix and clustering ...\n")
D = pairwise_distances(seqs_2d, metric=identity_dist, n_jobs=-1)
# print("distance matrix:\n{}\n".format(D))

dm = squareform(D)

Z = linkage(dm)
# print("linkage matrix:\n{}\n".format(Z))

T = fcluster(Z, t=20, criterion='distance')
print("Found {} clusters with max 20% different.\n".format(max(T)+1))

# print("flat cluster:\n{}\n".format(T))
cluster_file = args.index + '.T.json'
with open(cluster_file, 'w') as f:
    json.dump(T.tolist(), f, indent=4)
print('Flat cluster result save at {}\n'.format(cluster_file))

print('Elapsed time {}.'.format(dt.now() - start))

import json
import argparse
import numpy as np
from tqdm import tqdm
from pathlib import Path
from datetime import datetime as dt
start = dt.now()

from Bio import PDB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

fasta_file = args.index + '.fasta'
with open(fasta_file, 'w') as f:
    for code, seq in zip(codes, seqs):
        SeqIO.write(SeqRecord(seq, id=code, description=''), f, 'fasta')
print('Sequences save at {}\n'.format(fasta_file))
print('Elapsed time {}.'.format(dt.now() - start))

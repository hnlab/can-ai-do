import json
import argparse
import numpy as np
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from Bio import PDB
from Bio.Seq import Seq
from Bio import pairwise2 as pw2

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', '--index', required=True)
parser.add_argument(
    '-d', '--datadir', default='./v2018', help="datadir, default is ./v2018")
args = parser.parse_args()

def SeqFromPDBFile(_id, pdb_file, pocket_file=None):
    parser = PDB.PDBParser()
    chain_id = None
    if pocket_file:
        pocket = parser.get_structure(_id, pocket_file)
        for chain in pocket.get_chains():
            if not chain.id: continue
            if len(chain) > 10:
                chain_id = chain.id
                break
    
    ppb = PDB.PPBuilder()
    protein = parser.get_structure(_id, pdb_file)
    if chain_id:
        for chain in protein.get_chains():
            if chain.id == chain_id:
                seqs = [i.get_sequence() for i in ppb.build_peptides(chain)]
                seq_str = ''.join([str(i) for i in seqs])
                a = seqs[0].alphabet
                return Seq(seq_str, a)
    else:
        longest_chain = None
        for chain in protein.get_chains():
            if longest_chain is None or len(chain) > len(longest_chain):
                longest_chain = chain
        seq_str = ''.join([str(i) for i in seqs])
        a = seqs[0].alphabet
        return Seq(seq_str, a)

def identity_percent(seqA, seqB):
    # Question: python: find percentage of match between two sequences
    # https://www.biostars.org/p/208540/
    global_align = pw2.align.globalxx(seqA, seqB)
    # seq_length = min(len(seqA), len(seqB))
    seq_length = global_align[0][4]
    matches = global_align[0][2]
    return (matches / seq_length) * 100

def read_index(index_file):
    with open(index_file) as f:
        for i in f:
            if i[0] == '#': continue
            code, reso, year, pK, *others = i.split()
            yield code, float(pK)

codes = []
pKs = []
ligands = []
fps = []
clusters = []
clusters_codes = []
datadir = Path(args.datadir)
for code, pK in tqdm(read_index(args.index)):
    protein_pdb = datadir / code / (code + '_protein.pdb')
    pocket_pdb = datadir / code / (code + '_pocket.pdb')
    ligand_pdb = Path('ligands') / (code + '_ligand.pdb')
    if not ligand_pdb.exists():
        continue
    if not pocket_pdb.exists():
        continue
    if not protein_pdb.exists():
        continue
    mol = Chem.MolFromPDBFile(str(ligand_pdb))
    if mol is None:
        continue
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)
    codes.append(code)
    pKs.append(pK)
    ligands.append(mol)
    fps.append(fp)
    seqA = SeqFromPDBFile(code, protein_pdb, pocket_pdb)
    for ci, clust in enumerate(clusters):
        is_break = False
        for seqB in clust:
            if identity_percent(seqA, seqB) > 80:
                clust.append(seqA)
                clusters_codes[ci].append(code)
                is_break = True
                break
        if is_break:
            break
    else:
        clusters.append([seqA])
        clusters_codes.append([code])

# inplace
clusters_codes.sort(key=len, reverse=True)
for i, clust in enumerate(clusters):
    print(i, len(clust))
with open(args.index + '.clusters.json', 'w') as f:
    json.dump(clusters_codes, f, indent=4)

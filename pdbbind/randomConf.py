"""generate random conformations from pdbbind ligands for test.
"""
import json
import argparse
import numpy as np
from tqdm import tqdm
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', '--index', required=True)
parser.add_argument('-d',
                    '--datadir',
                    required=True,
                    help="pdbbind datadir, like v2018")
parser.add_argument(
    '-o',
    '--output',
    default='output')
args = parser.parse_args()

DATADIR = Path(args.datadir)
OUTPUT = Path(args.output)
path_2d = OUTPUT / '2D'
path_2d.mkdir(parents=True, exist_ok=True)
path_random = OUTPUT / 'random'
path_random.mkdir(parents=True, exist_ok=True)

def read_index(index_file):
    codes = []
    pKs = []
    with open(index_file) as f:
        for i in f:
            if i[0] == '#': continue
            code, _, _, pK, *_ = i.split()
            codes.append(code)
            pKs.append(float(pK))
    return codes, pKs

codes, pKs = read_index(args.index)
for code in tqdm(codes, smoothing=0):
    path = DATADIR / code / (code + '_ligand.pdb')
    if path.exists():
        mol = Chem.MolFromPDBFile(str(path))
        if mol:
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
            pdb = path_2d / code / (code + '_ligand.pdb')
            sdf = path_2d / code / (code + '_ligand.sdf')
            pdb.parent.mkdir(exist_ok=True)
            Chem.MolToPDBFile(mol, str(pdb))
            if Chem.MolFromPDBFile(str(pdb)) is None:
                print('Can not load', pdb)
            writer = Chem.SDWriter(str(sdf))
            writer.write(mol)
            writer.close()
            if next(Chem.SDMolSupplier(str(sdf))) is None:
                print('Can not load ', sdf)

            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            pdb = path_random / code / (code + '_ligand.pdb')
            sdf = path_random / code / (code + '_ligand.sdf')
            pdb.parent.mkdir(exist_ok=True)
            Chem.MolToPDBFile(mol, str(pdb))
            if Chem.MolFromPDBFile(str(pdb)) is None:
                print('Can not load', pdb)
            writer = Chem.SDWriter(str(sdf))
            writer.write(mol)
            writer.close()
            if next(Chem.SDMolSupplier(str(sdf))) is None:
                print('Can not load ', sdf)
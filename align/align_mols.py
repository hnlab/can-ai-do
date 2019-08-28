"""Align mols with Open3DAlign.
"""
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('mol', nargs='+')
parser.add_argument('-ref')
parser.add_argument('-o', default='align.pdb')
args = parser.parse_args()

mols = []
for mol_file in args.mol:
    mol = Chem.MolFromPDBFile(mol_file)
    if mol is None:
        continue
    else:
        mol.file_name = mol_file
        mols.append(mol)

if args.ref:
    ref = Chem.MolFromPDBFile(args.ref)
else:
    # NumAtoms = [m.GetNumAtoms() for m in mols]
    # ref = mols[NumAtoms.index(max(NumAtoms))]
    N = len(mols)
    score = np.zeros((N, N))
    for i in range(N):
        # copy probe mol
        probe = Chem.Mol(mols[i])
        for j in range(N):
            if i == j:
                score[i,j] = 0
            try:
                O3A = rdMolAlign.GetO3A(probe, mols[j])
                score[i,j] = O3A.Align()
            except ValueError as e:
                score[i,j] = 100
                print(mol.file_name)
    sum_score = score.sum(axis=0)
    ref = mols[sum_score.argmin()]
    
pdb = Chem.PDBWriter(args.o)
for mol in mols:
    try:
        O3A = rdMolAlign.GetO3A(mol, ref)
        score = O3A.Align()
        print(score)
        pdb.write(mol)
    except ValueError as e:
        print(mol.file_name)
        print(e)
pdb.close()

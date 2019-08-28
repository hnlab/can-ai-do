"""Align mols with Open3DAlign.
"""
import argparse
import numpy as np
from pathlib import Path

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, rdMolAlign

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('mol', nargs='+')
parser.add_argument('-tc', default=0.8)
parser.add_argument('-o', default='aligned')
args = parser.parse_args()
output = Path(args.o)
output.mkdir(exist_ok=True, parents=True)

mols = []
for mol_file in args.mol:
    mol = Chem.MolFromPDBFile(mol_file)
    if mol is None:
        continue
    else:
        mol.id = Path(mol_file).stem
        mols.append(mol)

fps = [AllChem.GetMorganFingerprint(m, 2) for m in mols]
simi = [DataStructs.BulkTanimotoSimilarity(i, fps) for i in fps]

pdbs = {}
# NumAtoms = [m.GetNumAtoms() for m in mols]
# ref = mols[NumAtoms.index(max(NumAtoms))]
N = len(mols)
score = np.zeros((N, N))
for i in range(N):
    # copy probe mol
    ref = mols[i]
    _break = False
    for j in range(N):
        if i == j:
            continue
        if simi[i][j] < args.tc:
            continue
        if ref.id not in pdbs:
            pdb_name = output / ref.id
            pdb_name = pdb_name.with_suffix('.pdb')
            pdbs[ref.id] = Chem.PDBWriter(str(pdb_name))
            pdbs[ref.id].write(ref)
        try:
            # probe = Chem.Mol(mols[j])
            probe = Chem.MolFromSmiles(Chem.MolToSmiles(mols[j]))
            probe = Chem.AddHs(probe)
            AllChem.EmbedMolecule(probe)
            probe = Chem.RemoveHs(probe)
            pdbs[ref.id].write(probe)
            O3A = rdMolAlign.GetO3A(probe, ref, maxIters=100)
            score = O3A.Align()
            pdbs[ref.id].write(probe)
            _break = True
            break
        except ValueError as e:
            print(e)
            print(mol[i].id, mol[j].id)
    if _break:
        break

for pdb in pdbs.values():
    pdb.close()

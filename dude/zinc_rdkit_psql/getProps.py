"""Get Props and save in smiles.
"""
import gzip
import argparse
from pathlib import Path
from datetime import datetime

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.DataStructs import BulkTanimotoSimilarity

import itertools
from multiprocessing import Pool


def getProp(mol_line):
    smiles, mol_id = mol_line.split()[:2]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    mol_id = mol_id.replace('ZINC', '')
    mol_id = int(mol_id)
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    q = Chem.GetFormalCharge(mol)
    return tuple([smiles, mol_id, mw, logp, rotb, hbd, hba, q])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--smiles", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    start = datetime.now()
    
    smiles = Path(args.smiles)
    if smiles.suffix == '.gz':
        f = gzip.open(smiles, 'rt')
    else:
        f = open(smiles, 'r')
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix == '.gz':
        smi = gzip.open(args.output, 'wt')
    else:
        smi = open(args.output, 'w')
    smi.write("smiles zinc_id mw logp rotb hbd hba q\n")
    fmt = '{} {} {:.2f} {:.2f} {} {} {} {}\n'
    pool = Pool(processes=4)
    N = 10000
    counter_smi = 0
    counter_mol = 0
    while True:
        print("{}: Saved {:8d}/{:8d} in {}".format(
            datetime.now() - start, counter_mol, counter_smi, args.output))
        results = pool.map(getProp, itertools.islice(f, N))
        counter_smi += N
        if results:
            for r in results:
                if r is None: continue
                smi.write(fmt.format(*r))
                counter_mol += 1
        else:
            break

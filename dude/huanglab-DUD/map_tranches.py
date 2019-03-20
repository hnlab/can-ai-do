"""Map ZINC12 files into tranches.
"""
__version__ = "0.1.0"
__author__ = "Jincai Yang, jincai.yang42@gmail.com"

import yaml
import argparse
import numpy as np
from time import time
from pathlib import Path

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors as D
from rdkit.Chem import rdMolDescriptors as CD

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-i",
    "--input",
    nargs='+',
    required=True,
    help="input smiles files"
)
parser.add_argument(
    "-o", "--output", required=True, help="output dir")
args = parser.parse_args()


def timer(start,end):
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds)

def get_prop_array(mol):
    mw = CD.CalcExactMolWt(mol)
    logp = Chem.Crippen.MolLogP(mol)
    rotb = D.NumRotatableBonds(mol)
    hbd = CD.CalcNumHBD(mol)
    hba = CD.CalcNumHBA(mol)
    q = Chem.GetFormalCharge(mol)
    return np.array([mw, logp, rotb, hbd, hba, q])

def map_tranche(mw, logp):
    name = 'ABCDEFGHIJK'
    mw_slice = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500]
    logp_slice = [-1, 0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    tranche = ''
    for i, mwi in enumerate(mw_slice):
        if mw <= mwi:
            tranche += name[i]
            break
    else:
        tranche += name[i+1]
    for i, logpi in enumerate(logp_slice):
        if logp <= logpi:
            tranche += name[i]
            break
    else:
        tranche += name[i]
    return tranche

start = time()
output = Path(args.output)
output.mkdir(parents=True, exist_ok=True)
name = 'ABCDEFGHIJK'
tranches_path = {}
for i in name:
    for j in name:
        t_path = output/(i+j)
        t_path.mkdir(exist_ok=True)
        tranches_path[i+j] = t_path
tranches_mols = {}
tranches_count = {}
write_count = 0
for smiles_file in args.input:
    for m in Chem.SmilesMolSupplier(str(smiles_file), titleLine=False):
        if m is None: continue
        prop = get_prop_array(m)
        t = map_tranche(*prop[:2])
        if t not in tranches_mols:
            tranches_mols[t] = []
            tranches_count[t] = 0
        tranches_mols[t].append(m)
        tranches_count[t] += 1
        if tranches_count[t] % 1000 == 0:
            smi_name = "{:04d}.smi".format(tranches_count[t]//1000)
            smi_path = tranches_path[t]/smi_name
            f = Chem.SmilesWriter(str(smi_path))
            for m in tranches_mols[t]:
                f.write(m)
                write_count += 1
                if write_count % 1000 == 0:
                    print("{} write {:10d} mols".format(timer(start,time()), write_count))
            f.close()
            tranches_mols[t] = []

for t, mols in tranches_mols.items():
    if len(mols) == 0: continue
    smi_name = "{:04d}.smi".format(tranches_count[t]//1000 + 1)
    smi_path = tranches_path[t]/smi_name
    f = Chem.SmilesWriter(str(smi_path))
    for m in mols:
        f.write(m)
        write_count += 1
        if write_count % 1000 == 0:
            print("{} write {:10d} mols".format(timer(start,time()), write_count))
    f.close()

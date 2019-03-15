"""Generating smiles decoys from ZINC.
"""
__version__ = "0.1.0"
__author__ = "Jincai Yang, jincai.yang42@gmail.com"

import yaml
import argparse
import numpy as np
from pathlib import Path

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors as D
from rdkit.Chem import rdMolDescriptors as CD
from collections import namedtuple

example_text = """Example:
    genDecoys.py -a target/actives.smi
    
"""
parser = argparse.ArgumentParser(
    description=__doc__,
    epilog=example_text,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-a",
    "--actives",
    nargs='+',
    required=True,
    help=
    "actives in smiles format, one target one file, will use Path().parts[-2] as target name"
)
parser.add_argument(
    "-n", default=50, type=int, help="number of decoys per active")
parser.add_argument(
    "-z", default=15, type=int, help="zinc version: 12 or 15, default is 15")
parser.add_argument(
    "-mw",
    default=125,
    type=int,
    help="molecular weight range, default: 125, meaning +/- 125")
parser.add_argument(
    "-logp",
    default=3.6,
    type=float,
    help="logP, default: 3.6, meaning +/- 3.6")
parser.add_argument(
    "-rotb",
    default=5,
    type=int,
    help="number of rotation bonds, default: 5, meaning +/- 5")
parser.add_argument(
    "-hbd",
    default=4,
    type=int,
    help="number of hydrogen bond donor, default: 4, meaning +/- 4")
parser.add_argument(
    "-hba",
    default=3,
    type=int,
    help="number of hydrogen bond acceptor, default: 3, meaning +/- 3")
parser.add_argument(
    "-q",
    default=2,
    type=int,
    help="net charge, default: 2, meaning +/- 2")
parser.add_argument(
    "-tc_same",
    default=0.6,
    type=float,
    help=
    "filter out similar decoys against SAME target base on Tanimoto Coefficient, default: 0.6"
)
parser.add_argument(
    "-tc_diff",
    default=0.6,
    type=float,
    help=
    "filter out similar decoys against DIFFERENT targets base on Tanimoto Coefficient, default: 1, meaning no filter"
)
parser.add_argument(
    "-o", "--output_dir", default="output", help="output dir, default: output")
args = parser.parse_args()

MAX_PROP_DIFF = np.array([
    args.mw, args.logp, args.rotb, args.hbd, args.hba, args.q])

def get_props(mol):
    mw = CD.CalcExactMolWt(mol)
    logp = Chem.Crippen.MolLogP(mol)
    rotb = D.NumRotatableBonds(mol)
    hbd = CD.CalcNumHBD(mol)
    hba = CD.CalcNumHBA(mol)
    q = Chem.GetFormalCharge(mol)
    return np.array([mw, logp, rotb, hbd, hba, q])


def next_smiles_file(mw, logp):
    # map mw and logp to tranche name
    name = 'ABCDEFGHIJK':
    mw_slice = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500]
    logp_slice = [-1, 0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    tranche_init = ''
    for i, mwi in enumerate(mw_slice):
        if mw <= mwi:
            name_mw = name[max(0,i-1):i+2]
            break
    else:
        tranche_mw = name[max(0,i-1):i+2]
    tranche_init += name[i]
    for i, logpi in enumerate(logp_slice):
        if logp <= logpi:
            name_p += name[max(0,i-1):i+2]
            break
    else:
        name_p += name[max(0,i-1):i+2]
    tranche_init += name[i]
    
    # return smiles files in match tranche first.
    tranche = ZINC_PATH/tranche_init
    for smi in tranche.glob("*/*.smi"):
        yield smi
    
    # return smiles files in neighbor tranches randomly.
    smi_files = []
    for n_mw in name_mw:
        for n_p in name_p:
            if n_mw + n_p == tranche_init:
                continue
            tranche = ZINC_PATH/(n_mw + n_p)
            for smi in tranche.glob("*/*.smi"):
                smi_files.append(smi)
    for i in np.random.permutation(len(smi_files)):
        yield smi_files[i]

# step 1: select decoys in range


# step 2.1: remove decoys with SAME ID against all ligands in same target
targets = []
actives = []
actives_fps = []
actives_props = []
for a_file in args.actives:
    a_file = Path(a_file)
    if len(a_file.parts) > 1:
        target = a_file.parts[-2]
    else:
        target = a_file.stem
    mols = [m for m in Chem.SmilesMolSupplier(str(a_file), titleLine=False) if m is not None]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024) for m in mols]
    props = [get_prop(m) for m in mols]
    targets.append(target)
    actives.append(mols)
    actives_fps.append(fps)
    actives_props.append(props)

# write actives


DONE = [False for i in targets]
decoys = [[] for i in targets]
decoys_fps = [[] for i in targets]
decoys_props = [[] for i in targets]
decoys_count = [[0 for a in t_as] for t_as in actives]
decoys_done = [[False for a in t_as] for t_as in actives]
while not all(DONE):
    for ti, target in enumerate(targets):
        if DONE[ti]: continue
        if all(decoys_done[ti]):
            DONE[ti] = True
            continue
        a_fps = actives_fps[ti]
        a_props = actives_props[ti]
        for ai in range(len(a_fps)):
            if decoys_done[ti][ai] == True:
                continue
            a_fp = a_fps[ai]
            a_prop = a_props[ai]
            for m in next_smiles_file(*prop[:2]):
                fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
                if DataStructs.TanimotoSimilarity(a_fp, fp) > args.tc:
                    continue
                prop = get_prop(m)
                diff = np.abs(prop-a_prop)
                if diff > MAX_PROP_DIFF:
                    continue
                a_simi = DataStructs.BulkTanimotoSimilarity(fp, a_fps)
                if max(a_simi) > args.tc:
                    continue
                _continue = False
                for tj, d_fps in decoys_fps:
                    if ti == tj:
                        max_tc = args.tc_same
                    else:
                        max_tc = args.tc_diff
                    if max_tc == 1: continue
                    for start in range(0, len(d_fps), 100):
                        fps = d_fps[start:start + 100]
                        d_simi = DataStructs.BulkTanimotoSimilarity(fp, fps)
                        if max(d_simi) > args.tc_same:
                            _continue = True
                            break
                    if _continue: break
                if _continue: continue
                decoys[ti].appned(m)
                decoys_fps[ti].append(fp)
                decoys_count[ti][ai] += 1
                if decoys_count[ti][ai] == args.num_decoys:
                    decoys_done[ti][ai] = True
                break
            else:
                # run out of zinc mols
                decoys_done[ti][ai] = True


# step 2.2: remove decoys with SAME ID in all ligands (Tc > 0.25) in same target
# step 2.3: remove decoys similar to decoys (Tc > 0.6) in same target

output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)
for i, target in enumerate(targets):
    tdir = output_dir / target
    tdir.mkdir(exist_ok=True)
    a_file = tdir / "actives_final.smi"
    d_file = tdir / "decoys_final.smi"
    a_smi = Chem.SmilesWriter(str(a_file), includeHeader=False)
    d_smi = Chem.SmilesWriter(str(d_file), includeHeader=False)
    for a in actives[i]:
        a_smi.write(a)
    for d in decoys[i]:
        d_smi.write(d)

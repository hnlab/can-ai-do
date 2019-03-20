"""Generating smiles decoys from ZINC.
"""
__version__ = "0.1.1"
__author__ = "Jincai Yang, jincai.yang42@gmail.com"

import yaml
import argparse
import numpy as np
from pathlib import Path
from datetime import datetime
start = datetime.now()

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors as D
from rdkit.Chem import rdMolDescriptors as CD

example_text = """Example:
    genDecoys.py -a target/actives.smi -z zinc_path -o output
    
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
    "-z", "--zinc_path", required=True, help="ZINC path, a dir for ZINC15, or file for ZINC12.")
parser.add_argument(
    "-n", "--num_decoys", default=50, type=int, help="number of decoys per active")
parser.add_argument(
    "-p", "--probability", default=1, type=float, help="Probability of accepting decoys, DUDE use 0.067")
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
    "-tc",
    default=0.35,
    type=float,
    help=
    "find dissimilar decoys base on Tanimoto Coefficient, default: 0.35"
)
parser.add_argument(
    "-tc_same",
    default=0.6,
    type=float,
    help=
    "filter out similar decoys against SAME target base on Tanimoto Coefficient, default: 0.6"
)
parser.add_argument(
    "-tc_diff",
    default=1,
    type=float,
    help=
    "filter out similar decoys against DIFFERENT targets base on Tanimoto Coefficient, default: 1, meaning no filter"
)
parser.add_argument(
    "-o", "--output", required=True, help="output dir")
args = parser.parse_args()


MAX_PROP_DIFF = np.array([
    args.mw, args.logp, args.rotb, args.hbd, args.hba, args.q])

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

def extend_tranches(tranche):
    "yield neighbor tranches randomly"
    name = 'ABCDEFGHIJK'
    name_mw, name_logp = tranche
    i = name.find(name_mw)
    j = name.find(name_logp)
    neighbor_tranches = []
    for n_mw in name[max(0, i-1):i+2]:
        for n_logp in name[max(0, j-1):j+2]:
            if n_mw + n_logp == tranche: continue
            neighbor_tranches.append(n_mw + n_logp)
    for t in np.random.permutation(neighbor_tranches):
        yield t

def tranche_supplier(zinc_path, tranche):
    zinc_path = Path(zinc_path)
    assert zinc_path.exists()
    if zinc_path.is_file():
        for m in Chem.SmilesMolSupplier(str(zinc_path), titleLine=False):
            if m is not None:
                yield m
    else:
        if tranche.upper() == "ALL":
            # yield mol from all smiles files
            for smi in zinc_path.rglob("*.smi"):
                for m in Chem.SmilesMolSupplier(str(smi)):
                    if m is not None:
                        yield m
        else:
            # yield mol from smiles files in match tranche first.
            tranche_path = zinc_path/tranche
            smi_files = list(tranche_path.rglob("*.smi"))
            for smi in np.random.permutation(smi_files):
                for m in Chem.SmilesMolSupplier(str(smi)):
                    if m is not None:
                        yield m
            # yield mol from smiles files in neighbor tranches randomly.
            for t in extend_tranches(tranche):
                t_path = zinc_path/t
                smi_files = list(t_path.rglob("*.smi"))
                for smi in np.random.permutation(smi_files):
                    for m in Chem.SmilesMolSupplier(str(smi)):
                        if m is not None:
                            yield m

def similar(fp, fps, max_tc, step=128):
    if max_tc >= 1:
        return False
    for idx in range(0, len(fps), step):
        simi = DataStructs.BulkTanimotoSimilarity(fp, fps[idx:idx + step])
        if max(simi) > max_tc:
            return True
    return False
# step 1: select decoys in range

output = Path(args.output)
output.mkdir(parents=True, exist_ok=True)
# output.mkdir(parents=True)

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
    props = []
    props = [get_prop_array(m) for m in mols]
    targets.append(target)
    actives.append(mols)
    actives_fps.append(fps)
    actives_props.append(props)
    print("{} loading {:7s} {:4d} actives from {}".format(datetime.now(), target, len(mols), a_file))

# write actives and group actives into tranches
decoys_smi = []
actives_tranches = {}
zinc_path = Path(args.zinc_path)
for i, target in enumerate(targets):
    tdir = output / target
    tdir.mkdir(exist_ok=True)
    a_file = tdir / "actives_final.smi"
    d_file = tdir / "decoys_final.smi"
    a_smi = Chem.SmilesWriter(str(a_file), includeHeader=False)
    for j, a in enumerate(actives[i]):
        a_smi.write(a)
        if zinc_path.is_file():
            tranche = 'ALL'
        else:
            prop = actives_props[i][j]
            tranche = map_tranche(*prop[:2])
        if tranche not in actives_tranches:
            actives_tranches[tranche] = []
        actives_tranches[tranche].append((i,j))
    a_smi.close()
    d_smi = Chem.SmilesWriter(str(d_file), includeHeader=False)
    decoys_smi.append(d_smi)
with open(output/"config.yaml", "w") as f:
    yaml.dump(args, f)

mol_count = 0
decoys_fps = [[] for i in targets]
decoys_props = [[] for i in targets]
discard_ids = [set() for i in targets]
global_discard_ids = set()
decoys_count = [[0 for a in t_a] for t_a in actives]
total_decoys = sum([args.num_decoys for t in actives for i in t])
actives_faild = [[] for i in targets]
for tranche, group_idx in actives_tranches.items():
    group_idx = np.array(group_idx)
    no_progress = 0
    for m in tranche_supplier(zinc_path, tranche):
        no_progress += 1
        if no_progress >= 10000:
            for t, i in group_idx:
                if decoys_count[t][i] < args.num_decoys:
                    actives_faild[t].append(actives[t][i])
            break
        if all([decoys_count[t][i] >= args.num_decoys for t, i in group_idx]):
            break
        if mol_count % 1000 == 0:
            sum_decoys = sum([i for t in decoys_count for i in t])
            print("{} generate {:8d}/{} from {:10d} zinc mols".format(
                datetime.now(), sum_decoys, total_decoys, mol_count))
            # print([decoys_count[t][i] for t, i in group_idx])
        mol_count += 1
        _id = m.GetProp("_Name")
        if args.tc_diff < 1 and _id in global_discard_ids:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
        prop = get_prop_array(m)
        for ti, ai in np.random.permutation(group_idx):
            if decoys_count[ti][ai] >= args.num_decoys:
                continue
            if _id in discard_ids[ti]:
                continue
            if np.random.rand() > args.probability:
                continue
            a_fp = actives_fps[ti][ai]
            if DataStructs.TanimotoSimilarity(fp, a_fp) > args.tc:
                continue
            a_prop = actives_props[ti][ai]
            diff = np.abs(prop-a_prop)
            if np.any(diff > MAX_PROP_DIFF):
                continue
            if similar(fp, actives_fps[ti], args.tc):
                discard_ids[ti].add(_id)
                continue
            if args.tc_same < 1 and similar(fp, decoys_fps[ti], args.tc_same):
                discard_ids[ti].add(_id)
                continue
            if args.tc_diff < 1:
                _continue = False
                for tj, d_fps in enumerate(decoys_fps):
                    if ti == tj:
                        continue
                    if similar(fp, d_fps, args.tc_diff):
                        global_discard_ids.add(_id)
                        _continue = True
                        break
                if _continue: continue
            decoys_smi[ti].write(m)
            decoys_fps[ti].append(fp)
            decoys_props[ti].append(prop)
            decoys_count[ti][ai] += 1
            discard_ids[ti].add(_id)
            if args.tc_diff < 1:
                global_discard_ids.add(_id)
            no_progress = 0
sum_decoys = sum([i for t in decoys_count for i in t])
print("{} generate {:8d}/{} from {:10d} zinc mols".format(
    datetime.now(), sum_decoys, total_decoys, mol_count))

for smi in decoys_smi:
    smi.close()

print("\nSummary:")
for t, target in enumerate(targets):
    a_fail_mols = actives_faild[t]
    a_fail = len(a_fail_mols)
    a_total = len(actives[t])
    a_done = a_total - a_fail
    d_count = sum(decoys_count[t])
    d_total = a_total * args.num_decoys
    print("{:7} generate {:6d}/{:6d} decoys for {:4}/{:4} actives, {} failures.".format(
        target, d_count,d_total, a_done, a_total, a_fail))
    if a_fail > 0:
        smi_name = output/target/"failed_actives.smi"
        smi = Chem.SmilesWriter(str(smi_name))
        for m in a_fail_mols:
            smi.write(m)
timedelta = datetime.now() - start
print("\nResult saved in {}\nTime elapsed {}".format(output, timedelta))

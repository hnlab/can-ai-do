"""Deduping decoys with same zinc id in same target.
"""
import yaml
import argparse
import numpy as np
from pathlib import Path

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from collections import namedtuple

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "-t", "--target_dir", nargs='+', required=True, help="target dir names")
parser.add_argument(
    "-c", "--config", required=True, help="config in yaml format")
parser.add_argument(
    "-d",
    "--decoys_glob",
    default="decoys/ligand*/*property_matched_decoys.txt",
    help=
    "decoys glob pattern, default: decoys/ligand*/*property_matched_decoys.txt"
)
parser.add_argument(
    "--output_name",
    default="decoys_final.ism",
    help="output base name, default: decoys_final.ism")
parser.add_argument(
    "-o", "--output_dir", default="output", help="output dir, default: output")
args = parser.parse_args()

# config.yaml
# DIFF_CUTOFF:
#     MW: 125
#     LogP: 3.6
#     RotB: 5
#     HBD: 4
#     HBA: 3
#     Q: 2
# DECOYS_PER_ACTIVE: 50
# ACTIVE_DECOY_TC: 0.25
# DD_TC_SAME_TARGET: 0.6
# DD_TC_CROSS_TARGET: ~
with open(args.config) as f:
    config = yaml.load(f)
diff = config["DIFF_CUTOFF"]
PROPERTY_DIFF_CUTOFF = np.array([
    diff['MW'], diff['LogP'], diff['RotB'], diff['HBD'], diff['HBA'], diff['Q']
])
DECOYS_PER_ACTIVE = config["DECOYS_PER_ACTIVE"]
ACTIVE_DECOY_TC = config["ACTIVE_DECOY_TC"]
DD_TC_SAME_TARGET = config["DD_TC_SAME_TARGET"]
DD_TC_CROSS_TARGET = config["DD_TC_CROSS_TARGET"]
# max Tc between actives decoys in same target


class Active(
        namedtuple('Active',
                   ["SMILES", "ID", "MW", "LogP", "RotB", "HBD", "HBA", "Q"])):
    __slots__ = ()

    @classmethod
    def from_str(self, mol_str):
        assert mol_str[:6] == 'LIGAND'
        data_types = (None, str, str, float, float, int, int, int, int)
        data = []
        for dtype, value in zip(data_types, mol_str.split()):
            if dtype is None:
                continue
            else:
                data.append(dtype(value))
        return self(*data)

    def to_array(self):
        return np.array(
            [self.MW, self.LogP, self.RotB, self.HBD, self.HBA, self.Q])


class Decoy(
        namedtuple('Decoy', [
            "SMILES", "ID", "MW", "LogP", "RotB", "HBD", "HBA", "Q",
            "Tc_to_Lig"
        ])):
    __slots__ = ()

    @classmethod
    def from_str(self, mol_str):
        assert mol_str[:5] == 'DECOY'
        data_types = (None, None, str, str, float, float, int, int, int, int,
                      None, float)
        data = []
        for dtype, value in zip(data_types, mol_str.split()):
            if dtype is None:
                continue
            else:
                data.append(dtype(value))
        return self(*data)

    def to_array(self):
        return np.array(
            [self.MW, self.LogP, self.RotB, self.HBD, self.HBA, self.Q])


# step 1: select decoys in range
# step 2.1: remove decoys with SAME ID against all ligands in same target
actives = {}
candidate_decoys = {}
for t in args.target_dir:
    print(t)
    t = Path(t)
    t_actives = []
    t_decoys = []
    actives[t.name] = t_actives
    candidate_decoys[t.name] = t_decoys
    uniq_ids = set()
    for decoys_file in t.glob(args.decoys_glob):
        print(decoys_file)
        f = open(decoys_file)
        for line in f:
            if line[:5] == "DECOY":
                decoy = Decoy.from_str(line)
                if decoy.ID in uniq_ids:
                    continue
                if decoy.Tc_to_Lig > ACTIVE_DECOY_TC:
                    continue
                property_diff = np.abs(decoy.to_array() - active.to_array())
                if np.any(property_diff > PROPERTY_DIFF_CUTOFF):
                    continue
                t_decoys.append(decoy)
            elif line[:6] == "LIGAND":
                active = Active.from_str(line)
                t_actives.append(active)
            else:
                continue
        f.close()

# step 2.2: remove decoys with SAME ID in all ligands (Tc > 0.25) in same target
# step 2.3: remove decoys similar to decoys (Tc > 0.6) in same target

Target = namedtuple('Target', [
    'name', 'actives', 'actives_fps', 'candidate_decoys', 'random_idx',
    'decoys', 'decoys_fps', 'decoys_count'
])
targets = []
for name, t_decoys in candidate_decoys.items():
    t_actives = actives[name]
    t_actives_fps = []
    for a in t_actives:
        m = Chem.MolFromSmiles(a.SMILES)
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
        t_actives_fps.append(fp)
    decoys_count = [0 for i in t_actives]
    random_idx = np.random.permutation(len(t_decoys))
    target = Target(name, t_actives, t_actives_fps, t_decoys, random_idx, [],
                    [], decoys_count)
    targets.append(target)

DECOYS_PER_ACTIVE = 50
DONE = {t.name: False for t in targets}
idx = {t.name: 0 for t in targets}
count = 0
while any([DONE[t.name] == False for t in targets]):
    count += 1
    if count % 1000 == 0:
        print("filter & assign decoys cycle {:7d}:".format(count))
        for t in targets:
            num_decoys = DECOYS_PER_ACTIVE*len(t.decoys_count)
            active_done = sum([i==DECOYS_PER_ACTIVE for i in t.decoys_count])
            print("target: {:>5s}, {:5d}/{} decoys assigned, {:3d}/{} actives done.".format(
                t.name, sum(t.decoys_count), num_decoys, active_done, len(t.decoys_count)))
        print()
    for ti in targets:
        name = ti.name
        if DONE[name]: continue
        if idx[name] >= len(ti.candidate_decoys):
            DONE[name] = True
            continue
        if min(ti.decoys_count) >= DECOYS_PER_ACTIVE:
            DONE[name] = True
            continue
        # use random idx in candidate decoys
        # to achieve random pick in matched decoys for each active
        # thus final decoys need to compare fps can be minimum,
        # because we can filter then with property first.
        decoy_idx = ti.random_idx[idx[name]]
        idx[name] += 1
        decoy = ti.candidate_decoys[decoy_idx]

        match_active_idx = None
        for ai, active in enumerate(ti.actives):
            property_diff = np.abs(decoy.to_array() - active.to_array())
            if np.any(property_diff > PROPERTY_DIFF_CUTOFF):
                continue
            if match_active_idx is None:
                match_active_idx = ai
                continue
            if ti.decoys_count[ai] < ti.decoys_count[match_active_idx]:
                match_active_idx = ai
        if ti.decoys_count[match_active_idx] >= DECOYS_PER_ACTIVE:
            continue

        # filter decoy with fp similarity
        m = Chem.MolFromSmiles(decoy.SMILES)
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
        a_simi = DataStructs.BulkTanimotoSimilarity(fp, ti.actives_fps)
        if max(a_simi) > ACTIVE_DECOY_TC:
            continue
        _continue = False
        for tj in targets:
            if ti.name == tj.name:
                max_cutoff = DD_TC_SAME_TARGET
            else:
                if DD_TC_CROSS_TARGET is None or DD_TC_CROSS_TARGET >= 1:
                    # no need to check
                    continue
                max_cutoff = DD_TC_CROSS_TARGET
            for start in range(0, len(tj.decoys_fps), 100):
                fps = tj.decoys_fps[start:start + 100]
                d_simi = DataStructs.BulkTanimotoSimilarity(fp, fps)
                if max(d_simi) > max_cutoff:
                    _continue = True
                    break
            if _continue: break
        if _continue: continue
        ti.decoys.append(decoy)
        ti.decoys_count[match_active_idx] += 1
        ti.decoys_fps.append(fp)

targets_dict = {t.name: t for t in targets}
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)
for tdir in args.target_dir:
    tdir = Path(tdir)
    target = targets_dict[tdir.name]
    out_tdir = output_dir / tdir.name
    out_tdir.mkdir(exist_ok=True)
    out_file = out_tdir / args.output_name
    with open(out_file, 'w') as f:
        for decoy in target.decoys:
            f.write("{} {}\n".format(decoy.SMILES, decoy.ID))

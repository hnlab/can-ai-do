import os
import sys
import json
import argparse
import numpy as np
import scipy.sparse as sp
import itertools
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem

zinc_file_names = []
zinc_drug_like_nums = []
print("count ZINC mols ...")
with open("wc_line.txt") as f:
    for line in f:
        line = line.strip()
        line_num, file_name = line.split()
        zinc_file_names.append(file_name)
        zinc_drug_like_nums.append(int(line_num)-1) # smiles with title
zinc_drug_like_nums = np.array(zinc_drug_like_nums)
zinc_file_names = np.array(zinc_file_names)
valid_smi = zinc_drug_like_nums > 0
zinc_drug_like_nums = zinc_drug_like_nums[valid_smi]
zinc_file_names = zinc_file_names[valid_smi]

print("{} drug like mols.".format(sum(zinc_drug_like_nums)))

datadir="../all"
fold_file="../2split/crossFamilySplit/family5fold.json"
with open(fold_file) as f:
    folds = json.load(f)
targets = np.array([name for k, f in folds.items() for name in f])
decoy_nums = []
for name in targets:
    file_name = os.path.join(datadir, name)
    decoy_dude = os.path.join(file_name, 'decoys_final.ism')
    with open(decoy_dude) as f:
        decoy_nums.append(sum(1 for line in f)) # smiles without title

np.random.seed(123)
indices = np.random.choice(sum(zinc_drug_like_nums), sum(decoy_nums))
indices = np.sort(indices)
all_decoy_temp = tempfile.TemporaryFile('w+')
count = 0
write_count = 0
for num_zinc, zinc_file in zip(zinc_drug_like_nums, zinc_file_names):
    tmp_idx = indices[np.logical_and(indices>=count, indices<count+num_zinc)]
    tmp_idx -= count - 1
    count += num_zinc
    if len(tmp_idx) == 0: continue
    f = open(zinc_file)
    idx_i = 0
    for i, line in enumerate(f):
        if idx_i == len(tmp_idx): continue
        if tmp_idx[idx_i] <= i:
            smiles, name = line.rstrip().split()
            m = Chem.MolFromSmiles(smiles)
            if m is None:
                print(smiles)
                continue
            all_decoy_temp.write(line)
            write_count += 1
            idx_i += 1
    f.close()
    print('picked {:8d}/{}'.format(write_count, sum(decoy_nums)))

start = 0
indices = np.random.permutation(sum(decoy_nums))
for name, decoy_num in zip(targets, decoy_nums):
    all_decoy_temp.seek(0)
    file_name = os.path.join(datadir, name)
    decoy_random = os.path.join(file_name, 'decoys_random.smi')
    print("{} needs {} decoys.".format(name, decoy_num))
    decoyf = open(decoy_random, 'w')
    tmp_idx = indices[start:start+decoy_num]
    start += decoy_num
    tmp_idx = np.sort(tmp_idx)
    count = 0
    write_count = 0
    idx_i = 0
    for i, line in enumerate(all_decoy_temp):
        if idx_i == len(tmp_idx): continue
        if tmp_idx[idx_i] <= i:
            decoyf.write(line) # smiles without title
            write_count += 1
            idx_i += 1
    print("picked {} decoys in {}.\n".format(write_count, decoy_random))
    decoyf.close()
all_decoy_temp.close()

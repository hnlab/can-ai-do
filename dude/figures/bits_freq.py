"""convert smi into fingerprints and then count bits freq. About 1 hours for ZINC12.
"""
import json
import argparse
import numpy as np
from tqdm import tqdm
from pathlib import Path
import multiprocessing as mp

from rdkit import Chem
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-s", "--smi", required=True, help="smiles_files")
parser.add_argument("-n", type=int, default=2048, help="default is 2048.")
args = parser.parse_args()
args.smi = Path(args.smi)


def mfp2(m):
    # radius 2 MorganFingerprint equal ECFP4
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=args.n)
    return fp


def smiles2fp(line):
    smiles, name, *_ = line.split()
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    fp = mfp2(m)
    return fp


# http://rdkit.blogspot.com/2016/02/morgan-fingerprint-bit-statistics.html
freq = [0 for i in range(args.n)]
vaild_count = 0
with open(args.smi) as f:
    total = sum(1 for line in f)
pbar = tqdm(desc='count bits freq', total=total, unit=' mol')
with open(args.smi) as f:
    # 5000000 chars about 100000 lines 
    nchar = 5000000
    lines = f.readlines(nchar)
    while lines:
        with mp.Pool() as p:
            fps = p.map(smiles2fp, lines)
        for fp in fps:
            if fp is None:
                continue
            vaild_count += 1
            for i in fp.GetOnBits():
                freq[i] += 1
        pbar.update(len(lines))
        lines = f.readlines(nchar)
pbar.close()

fps_freq = np.array(freq) / vaild_count
freq_file = args.smi.with_suffix(args.smi.suffix + '.freq.json')
with open(freq_file, 'w') as f:
    json.dump(fps_freq.tolist(), f)
print(f'fps freq saved to {freq_file}')
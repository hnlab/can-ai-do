"""plot line chart of bits freq.
"""
import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--dude_freq', required=True)
parser.add_argument('--zinc_freq', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

args.output = Path(args.output)
df = pd.read_csv(args.dude_freq)
with open(args.zinc_freq) as f:
    zinc_freq = json.load(f)
df['zinc_freq'] = zinc_freq

df = df.sort_values(by='zinc_freq', ascending=False)
df = df.reset_index(drop=True)
print(df.tail())
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(df.index, df.active_freq, color='red', label='Actives', marker='+')
ax.scatter(df.index, df.decoy_freq, color='blue', label='Decoys', marker='+')
ax.scatter(df.index, df.zinc_freq, color='black', label='ZINC', marker='+')
# ax.legend()
ax.set_ylim([-0.05, 1.05])
ax.set_xlabel('2048 bits of Morgan Fingerprint sorted by frequency in ZINC')
ax.set_ylabel('frequency')
ax2 = fig.add_axes([0.29,0.265,0.6,0.6])
ax2.scatter(df.index, df.active_freq, color='red', label='Actives', marker='+')
ax2.scatter(df.index, df.decoy_freq, color='blue', label='Decoys', marker='+')
ax2.scatter(df.index, df.zinc_freq, color='black', label='ZINC', marker='+')
ax2.legend()
ax2.set_ylim([-0.005, 0.1])
jpg = args.output.with_suffix('.all_bits' + args.output.suffix)
fig.savefig(jpg, dpi=300)
print(f'all bits saved to {jpg}')

df = df.loc[(df.mean_freq >= 0.03) & (np.abs(df.log2_FC) >= 1)]
df = df.sort_values(by='zinc_freq', ascending=False)
df = df.reset_index(drop=True)
fig, ax = plt.subplots(figsize=(8, 6))
for i, (lower, upper) in enumerate(zip(df.active_freq, df.decoy_freq)):
    if lower > upper:
        lower, upper = upper, lower
    ax.vlines(i, lower, upper, colors='silver', zorder=0)
ax.scatter(df.index, df.active_freq, color='red', label='Actives', marker='+')
ax.scatter(df.index, df.decoy_freq, color='blue', label='Decoys', marker='+')
ax.scatter(df.index, df.zinc_freq, color='black', label='ZINC', marker='+')
ax.set_xlabel(f'{len(df)} high frequency and change significant bits in DUD-E')
ax.set_ylabel('frequency')
ax.legend()
fig.savefig(args.output, dpi=300)
print(f'significant bits plot saved to {args.output}')

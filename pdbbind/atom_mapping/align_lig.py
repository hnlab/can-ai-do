#!/usr/bin/python3
"""align ligand base on protein align.
"""
import os
import random
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-l',
                    '--ligand',
                    required=True,
                    help='ligand need to be aligned')
parser.add_argument('-p',
                    '--protein',
                    required=True,
                    help='the protein which ligand binding')
parser.add_argument('-r',
                    '--reference',
                    required=True,
                    help='reference protein')
args = parser.parse_args()

residues = []
with open(args.ligand) as f:
    for line in f:
        if len(line) > 6 and line[:6] in ('ATOM  ', 'HETATM'):
            residue = line[17:20]
            residues.append(residue)
residues = set(residues)
print(residues)

args.protein = Path(args.protein).resolve()
args.ligand = Path(args.ligand).resolve()
args.reference = Path(args.reference).resolve()
lig_output = args.ligand.with_suffix('.align.pdb')
pro_output = args.protein.with_suffix('.align.pdb')

radnom_str = ''.join(random.choices("abcdefghijklmnopqrstuvwxyz", k=10))
command_file = f'/tmp/.chimera.{radnom_str}.com'

with open(command_file, 'w') as f:
    f.write(f'open {args.reference} {args.protein} {args.ligand}\n'
            f'combine #1,2\n'  #1protein + #2ligand -> #3 combination
            f'mm #0 #3\n'
            f'select #3:{":".join(residues)}\n'
            f'write relative #0 selected #3 {lig_output}\n'
            f'select invert sel\n'
            f'write relative #0 selected #3 {pro_output}\n'
            f'stop\n')

os.system(f"chimera --nogui {command_file}")
# os.system(f"chimera {command_file}")
os.remove(command_file)

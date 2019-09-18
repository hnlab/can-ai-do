#!/usr/bin/python3
"""adjust atom radius based on occupancy.
"""
import os
import random
import argparse
import numpy as np
from pathlib import Path
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "-i",
    "--input",
    required=True,
    nargs="+",
    help="input pdb file, output chimera session as input.view.py")
parser.add_argument('--ribbon',
                    action="store_true",
                    help="open file one more time and show as ribbon")
args = parser.parse_args()

radnom_str = ''.join(random.choices("abcdefghijklmnopqrstuvwxyz", k=10))
command_file = f'/tmp/.chimera.{radnom_str}.com'
with open(command_file, 'w') as f:
    for i, pdb_file in enumerate(args.input):
        pdb_file = Path(pdb_file).resolve()
        f.write(f"open {pdb_file}\n")
        if args.ribbon:
            model_id = i * 2
        else:
            model_id = i
        f.write(f"~rib #{model_id}\n")
        f.write(f"disp #{model_id}\n")
        f.write("repr bs\n")
        f.write('labelopt info %(name)s %(occupancy).2f\n')
        indices = []
        occupancies = []
        for line in open(pdb_file):
            if len(line) >= 6 and line[:6] in ("ATOM  ", "HETATM"):
                index = int(line[6:11])
                occ = float(line[54:60])
                indices.append(index)
                occupancies.append(occ)
                f.write(f'select #{model_id}@/serialNumber={index}\n')
                if occ >= 0.3:
                    f.write(f"label sel\n")
                    f.write(f"color blue ,la sel\n")
                if occ <= -0.1:
                    f.write(f"label sel\n")
                    f.write(f"color red ,la sel\n")
                if occ >= -0.01:
                    r = occ * 6 + 0.8
                else:
                    r = -occ * 6 + 0.8
                    f.write(f"setattr a color black sel\n")
                f.write(f"setattr a radius {r:.3} sel\n")
        # label highest 1 and lowest 1
        rank = np.argsort(occupancies)
        for rank_idx in (rank[0], rank[-1]):
            occ = occupancies[rank_idx]
            index = indices[rank_idx]
            f.write(f'select #{model_id}@/serialNumber={index}\n')
            if occ >= 0:
                f.write(f"label sel\n")
                f.write(f"color blue ,la sel\n")
            else:
                f.write(f"label sel\n")
                f.write(f"color red ,la sel\n")
        if args.ribbon:
            # f.write("color byhet\n")
            f.write(f"open {pdb_file}\n")
            f.write(f"~disp #{model_id+1}\n")
    f.write("setattr g display false\n")
    f.write(f"save {pdb_file}.view.py\n")
    f.write("stop\n")

    # without hetatom color with --nogui
    # os.system("chimera --nogui /tmp/.chimera.com")
os.system(f"chimera {command_file}")
os.remove(command_file)
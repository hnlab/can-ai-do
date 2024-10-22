#!/usr/bin/python3
"""adjust atom radius based on occupancy.
"""
import os
import random
import argparse
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
        f.write(f"surface #{model_id}\n")
        f.write('labelopt info %(name)s %(occupancy).2f\n')
        for line in open(pdb_file):
            if len(line) >= 6 and line[:6] in ("ATOM  ", "HETATM"):
                index = int(line[6:11])
                f.write(f'select #{model_id}@/serialNumber={index}\n')
                occ = float(line[54:60])
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
                f.write(
                    f"setattr a charge {occ:.3} sel\n")
        f.write(f'coulombic atoms #{model_id} -10 red 0 white 10 blue\n')
        if args.ribbon:
            # f.write("color byhet\n")
            f.write(f"open {pdb_file}\n")
            f.write(f"~disp #{model_id+1}\n")
    f.write("setattr g display false\n")
    f.write(f"save {pdb_file}.charge.py\n")
    f.write("stop\n")

    # without hetatom color with --nogui
    # os.system("chimera --nogui /tmp/.chimera.com")
os.system(f"chimera {command_file}")
os.remove(command_file)
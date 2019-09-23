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
        f.write(f"~rib #{model_id}\n"
                f"disp #{model_id}\n"
                f"repr bs #{model_id}\n"
                # f'labelopt info %(occupancy).2f\n'
                f'labelopt info %(name)s %(occupancy).2f\n')
        # set red sphere background to atoms with occupancy < -0.01
        bg_id = model_id + 100
        f.write(f"open #{bg_id} {pdb_file}\n"
                f"~rib #{bg_id}\n"
                f"~disp #{bg_id}\n"
                f"~bond #{bg_id}\n"
                f"select #{bg_id}@/occupancy<-0.01\n"
                f"disp sel\n"
                # f"repr sphere sel\n"
                f"repr bs sel\n"  #  ball-and-stick
                # f"setattr a radius 0.7 sel\n"
                f"color black ,a sel\n"
                f"transparency 80 sel\n")
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
                    f.write(f"color black ,la sel\n")
                if occ <= -0.1:
                    f.write(f"label sel\n")
                    f.write(f"color red ,la sel\n")
                if occ >= 0:
                    r = occ * 6 + 0.8
                else:
                    r = -occ * 6 + 0.8
                f.write(f"setattr a radius {r:.3} sel\n")
                if occ < -0.01:
                    f.write(
                        f"setattr a radius {r*1.5} #{bg_id}@/serialNumber={index}\n"
                    )
        # label highest 1 and lowest 1
        rank = np.argsort(occupancies)
        for rank_idx in (rank[0], rank[-1]):
            occ = occupancies[rank_idx]
            index = indices[rank_idx]
            f.write(f'select #{model_id}@/serialNumber={index}\n')
            if occ >= 0:
                f.write(f"label sel\n")
                f.write(f"color black ,la sel\n")
            else:
                f.write(f"label sel\n")
                f.write(f"color red ,la sel\n")
        if args.ribbon:
            # f.write("color byhet\n")
            f.write(f"open {pdb_file}\n")
            f.write(f"~disp #{model_id+1}\n")

    f.write("setattr g display false\n"
            "~select\n"
            f"save {pdb_file}.view.py\n"
            "stop\n")

    # without hetatom color with --nogui
    # os.system("chimera --nogui /tmp/.chimera.com")
os.system(f"chimera {command_file}")
os.remove(command_file)
print(f"output chimera seesion file on {pdb_file}.view.py")
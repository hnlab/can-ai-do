#!/usr/bin/python3
"""adjust atom radius based on occupancy.
"""
import os
import argparse
from pathlib import Path
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "-i",
    "--input",
    required=True,
    nargs="+",
    help="input pdb file, output chimera session as input.view.py")
args = parser.parse_args()

for pdb_file in args.input:
    pdb_file = Path(pdb_file).resolve()
    # print(pdb_file.name)
    command_file = f"/tmp/.chimera.{pdb_file.name}.com"
    with open(command_file, 'w') as f:
        f.write(f"open {pdb_file}\n")
        f.write("~rib\n")
        f.write("disp\n")
        f.write("repr bs\n")
        for line in open(pdb_file):
            if len(line) >= 6 and line[:6] in ("ATOM  ", "HETATM"):
                index = int(line[6:11])
                occ = float(line[54:60])
                if occ >= -0.01:
                    r = occ * 8 + 0.8
                else:
                    r = -occ * 0.8 + 0.8
                    f.write(f"setattr a color black @/serialNumber={index}\n")
                f.write(f"setattr a radius {r:.3} @/serialNumber={index}\n")

        # f.write("color byhet\n")
        f.write(f"open {pdb_file}\n")
        f.write(f"~disp #1\n")
        f.write("setattr g display false\n")
        f.write(f"save {pdb_file}.view.py\n")
        f.write("stop\n")

    # without hetatom color with --nogui
    # os.system("chimera --nogui /tmp/.chimera.com")
    os.system(f"chimera {command_file}")
    os.remove(command_file)
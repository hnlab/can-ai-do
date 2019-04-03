"""Split all files into small parts.
"""
import gzip
import shutil
import tempfile
import argparse
import subprocess
from pathlib import Path
from datetime import datetime as dt
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-s", "--smiles", nargs='+', required=True)
parser.add_argument("-l", "--length", default=1000000, type=int, help="number of lines in each part, default 1,000,000")
parser.add_argument("-o", "--output_dir", required=True)
args = parser.parse_args()


def page_generator_from_files(smiles_files, length):
    page = []
    count = 0
    for smi in smiles_files:
        print("{}: loading {}".format(dt.now(), smi))
        smi = Path(smi)
        if smi.suffix == '.gz':
            f = gzip.open(smi, 'rt')
        else:
            f = open(smi, 'r')
        for i in f:
            if count < length:
                page.append(i)
                count += 1
            else:
                yield page
                page = []
                count = 0
    if page:
        yield page

start = dt.now()

output = Path(args.output_dir)
file_count = 0
for page in page_generator_from_files(args.smiles, args.length):
    file_count += 1
    # save under dirs separately for job submitting.
    file_dir = output / 'd.{:04d}'.format(file_count)
    file_dir.mkdir(parents=True, exist_ok=True)
    file_path = file_dir / "{}.{:04d}.smi".format(output.stem, file_count)
    with open(file_path, 'w') as f:
        f.writelines(page)
    print("{}: saved {}".format(dt.now(), file_path))
print("Total elapsed time: {}".format(dt.now()-start))

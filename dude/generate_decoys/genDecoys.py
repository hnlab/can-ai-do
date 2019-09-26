"""Generating decoys from database based on active smiles.
"""
import json
import gzip
import random
import argparse
import numpy as np
from pathlib import Path
from datetime import datetime as dt

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

import psycopg2
import psycopg2.extras
psycopg2.extensions.set_wait_callback(psycopg2.extras.wait_select)

from tqdm import tqdm
from multiprocessing.dummy import Pool

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-a", "--actives_file", required=True)
parser.add_argument("-n",
                    default=50,
                    type=int,
                    help="n decoys per acive. default 50")
parser.add_argument("-d", "--dbname", default="zinc")
parser.add_argument("-s", "--schema", default="zinc")
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
parser.add_argument("-P", "--processes", type=int)
parser.add_argument("-m", "--match_level", type=int, default=3)
parser.add_argument("-o", "--output", required=True, help="output dir")
parser.add_argument("--drug_like", action="store_true")
args = parser.parse_args()

start = dt.now()

args.subset = 'full'
if args.drug_like:
    args.subset = 'drug_like'


def getProp(mol_line):
    smiles, mol_id = mol_line.split()[:2]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    q = Chem.GetFormalCharge(mol)
    return tuple([mol_id, smiles, mw, logp, rotb, hbd, hba, q])


def random_name(length=10):
    import random
    import string
    s = string.ascii_letters + string.digits
    return ''.join(random.choices(s, k=length))


# It is long but simple, select num decoys based range from narrow to wide
# DEFAULTS FROM MYSINGER - these are used if property not specified
PROP_DIFF = np.array([
    [-1, 20, 35, 50, 65, 80, 100],  #"MWT_RANGES_LOWER"
    [20, 35, 50, 65, 80, 100, 125],  #"MWT_RANGES_UPPER"
    [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6],  #"LOGP_RANGES"
    [1, 2, 2, 3, 3, 4, 5],  #"RB_RANGES"
    [0, 0, 1, 1, 2, 2, 3],  #"HBD_RANGES"
    [0, 1, 2, 2, 3, 3, 4],  #"HBA_RANGES"
    [0, 0, 0, 0, 0, 1, 2],  #"CHG_RANGES"
])


def generate_decoys(active):
    _connect = psycopg2.connect(host=args.host,
                                dbname=args.dbname,
                                port=args.port)
    _cursor = _connect.cursor()
    _cursor.execute(
        f"CREATE TEMP TABLE decoys (like {args.schema}.{args.subset});"
        f"SET enable_seqscan TO OFF;")
    mol_id, smiles, mw, logp, rotb, hbd, hba, q = active
    for i, level in enumerate(PROP_DIFF.T):
        if i >= args.match_level:
            continue
        mw_lower, mw_upper, logp_diff, rotb_diff, hbd_diff, hba_diff, q_diff = level
        _cursor.execute(
            f"INSERT INTO decoys\n"
            f"SELECT * \n"
            f"  FROM {args.schema}.{args.subset}\n"
            f" WHERE ABS (mw - {mw}) > {mw_lower}\n"
            f"   AND ABS (mw - {mw}) <= {mw_upper}\n"
            f"   AND ABS (logp - {logp}) <= {logp_diff}\n"
            f"   AND ABS (rotb - {rotb}) <= {rotb_diff}\n"
            f"   AND ABS (hbd - {hbd}) <= {hbd_diff}\n"
            f"   AND ABS (hba - {hba}) <= {hba_diff}\n"
            f"   AND ABS (q - {q}) <= {q_diff}\n"
            # f" ORDER BY db.mfp2 <%> morganbv_fp(mol_from_smiles('{smiles}'::cstring)) DESC\n"
            f" LIMIT 3000 - (SELECT COUNT(*) FROM decoys);")
    _cursor.execute(
        f"INSERT INTO {args.schema}.{decoys_table}\n"
        f"SELECT zinc_id, smiles, mw, logp, rotb, hbd, hba, q, max(tanimoto_sml(D.mfp2,A.mfp2))\n"
        f"  FROM decoys AS D\n"
        f" CROSS JOIN {args.schema}.{actives_table}_fps AS A\n"
        f" GROUP BY zinc_id, smiles, mw, logp, rotb, hbd, hba, q\n"
        f"    ON CONFLICT (zinc_id) DO NOTHING\n")
    _connect.commit()
    _connect.close()


connect = psycopg2.connect(host=args.host, dbname=args.dbname, port=args.port)
cursor = connect.cursor()

# loading actives
actives_file = Path(args.actives_file)
actives_props = []
if len(actives_file.parts) > 2:
    # path/to/target/actives_final.smi
    target = actives_file.parts[-2]
else:
    # target.smi
    target = actives_file.stem
if actives_file.suffix == '.gz':
    f = gzip.open(actives_file, 'rt')
else:
    f = open(actives_file, 'r')
for p in map(getProp, f):
    # props: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    if p is not None:
        actives_props.append(p)
print(
    f"{dt.now()}: loading {target} {len(actives_props)} actives from {actives_file}"
)

# save actives into database
tmp_name = f"tmp_{random_name()}"
actives_table = f"{tmp_name}_{target}_actives"
decoys_table = f"{tmp_name}_{target}_decoys"
cursor.execute(f"CREATE TABLE IF NOT EXISTS {args.schema}.{actives_table} ("
               f" mol_id text PRIMARY KEY,"
               f" smiles text,"
               f"     mw real,"
               f"   logp real,"
               f"   rotb smallint,"
               f"    hbd smallint,"
               f"    hba smallint,"
               f"      q smallint)")

psycopg2.extras.execute_values(
    cursor,
    f'INSERT INTO {args.schema}.{actives_table} VALUES %s ON CONFLICT (mol_id) DO NOTHING',
    actives_props,
    template=None,
    page_size=128)
cursor.execute(f"SELECT morganbv_fp(mol_from_smiles(smiles::cstring)) as mfp2"
               f"  INTO {args.schema}.{actives_table}_fps"
               f"  FROM {args.schema}.{actives_table};")

cursor.execute(f"CREATE TABLE IF NOT EXISTS {args.schema}.{decoys_table} ("
               f"zinc_id integer PRIMARY KEY,"
               f" smiles text,"
               f"     mw real,"
               f"   logp real,"
               f"   rotb smallint,"
               f"    hbd smallint,"
               f"    hba smallint,"
               f"      q smallint,"
               f" max_tc real);")
connect.commit()

# get decoys for each active
output = Path(args.output)
output.mkdir(parents=True, exist_ok=True)
with open(output / 'args.json', 'w') as f:
    json.dump(vars(args), f, sort_keys=True, indent=4)

tdir = output / target
tdir.mkdir(exist_ok=True)
actives_file = tdir / "actives_final.smi"
decoys_file = tdir / "decoys_final.smi"
# prop: mol_id, smiles, mw, logp, rotb, hbd, hba, q
decoys_file = tdir / "decoys_final.smi"
with open(actives_file, 'w') as f:
    for p in actives_props:
        # p: mol_id, smiles, mw, logp, rotb, hbd, hba, q
        f.write(f'{p[1]} {p[0]}\n')

pool = Pool(processes=args.processes)
print(f"generating decoys against target {target}:")
for decoys in tqdm(pool.imap_unordered(generate_decoys, actives_props),
                   desc=f"Decoys for {target}",
                   total=len(actives_props),
                   smoothing=0):
    pass
pool.close()

cursor.execute(
    f"SELECT * \n"
    f"  FROM {args.schema}.{decoys_table}\n"
    f"ORDER BY max_tc\n"
    f"LIMIT (SELECT (COUNT(*) / 4) FROM {args.schema}.{decoys_table})")
decoys = cursor.fetchall()

cursor.execute(f"DROP TABLE IF EXISTS {args.schema}.{actives_table};\n"
               f"DROP TABLE IF EXISTS {args.schema}.{actives_table}_fps;\n"
               f"DROP TABLE IF EXISTS {args.schema}.{decoys_table};\n")
connect.close()

print(f"max_tc between actives and decoys {max([i[-1] for i in decoys])}")
decoys_list = [[] for i in actives_props]
decoys_num = np.zeros(len(actives_props))
decoys_unfilled = np.ones(len(actives_props), dtype=bool)
# active: mol_id, smiles, mw, logp, rotb, hbd, hba, q
actives_props_np = np.array([i[2:8] for i in actives_props])
for decoy in decoys:
    for i, level in enumerate(PROP_DIFF.T):
        if i >= args.match_level:
            continue
        # mw_lower, mw_upper, logp_diff, rotb_diff, hbd_diff, hba_diff, q_diff = level
        diff_upper = level[1:]
        # decoy: mol_id, smiles, mw, logp, rotb, hbd, hba, q, max_tc
        diff = np.abs(actives_props_np - np.abs(decoy[2:8]))
        match_mask = np.all(diff <= diff_upper, axis=1)
        match_mask = match_mask & decoys_unfilled
        if sum(match_mask) == 0:
            continue
        match_ind = np.flatnonzero(match_mask)
        match = match_ind[np.argmin(decoys_num[match_ind])]
        decoys_list[match].append((decoy[1], decoy[0]))
        decoys_num[match] += 1
        if decoys_num[match] >= 750:
            decoys_unfilled[match] = False
        break

with open(decoys_file, 'w') as f:
    for ds in decoys_list:
        if len(ds) > 50:
            ds = random.choices(ds, k=50)
        for decoy in ds:
            smiles, mol_id = decoy
            f.write(f"{smiles} ZINC{mol_id}\n")
print(f"Total elapsed time: {dt.now()-start}")

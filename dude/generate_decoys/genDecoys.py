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
parser.add_argument("-u", "--user")
parser.add_argument("-s", "--schema", default="zinc")
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
parser.add_argument("-P", "--processes", type=int)
parser.add_argument("-m", "--match_level", type=int, default=5)
parser.add_argument("-o", "--output", required=True, help="output dir")
parser.add_argument("--drug_like", action="store_true")
parser.add_argument("--heavyMW500", action="store_true")
args = parser.parse_args()

start = dt.now()

args.subset = 'full'
if args.drug_like:
    args.subset = 'drug_like'


def getProp(mol):
    mol_id = mol.GetProp("_Name")
    if mol is None: return None
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    mw = Descriptors.MolWt(mol)
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
PROP_DIFF = PROP_DIFF.T


def generate_decoys(active):
    _connect = psycopg2.connect(host=args.host,
                                dbname=args.dbname,
                                user=args.user,
                                port=args.port)
    _cursor = _connect.cursor()
    _cursor.execute(f"SET rdkit.tanimoto_threshold = 0.35;"
                    f"SET enable_seqscan TO OFF;")
    mol_id, smiles, mw, logp, rotb, hbd, hba, q = active
    for i, level in enumerate(PROP_DIFF):
        if i > args.match_level:
            continue
        mw_lower, mw_upper, logp_diff, rotb_diff, hbd_diff, hba_diff, q_diff = level
        _cursor.execute(
            f"DROP TABLE IF EXISTS decoys;"
            f"CREATE TEMP TABLE decoys (like {args.schema}.{args.subset});"
            f"INSERT INTO decoys\n"
            f"SELECT * \n"
            f"  FROM {args.schema}.{args.subset} D\n"
            f" WHERE ABS (mw - {mw}) > {mw_lower}\n"
            f"   AND ABS (mw - {mw}) <= {mw_upper}\n"
            f"   AND ABS (logp - {logp}) <= {logp_diff}\n"
            f"   AND ABS (rotb - {rotb}) <= {rotb_diff}\n"
            f"   AND ABS (hbd - {hbd}) <= {hbd_diff}\n"
            f"   AND ABS (hba - {hba}) <= {hba_diff}\n"
            f"   AND ABS (q - {q}) <= {q_diff}\n"
            f"   AND NOT EXISTS ("
            f"       SELECT FROM {args.schema}.{actives_table}_fps AS A\n"
            f"        WHERE A.mfp2 % D.mfp2)\n"
            # f"ORDER BY RANDOM()\n"
            f" LIMIT GREATEST(0, 100 - (SELECT COUNT(*) FROM {args.schema}.{decoys_table} WHERE mol_id = '{mol_id}'));\n"
            f"\n"
            f"INSERT INTO {args.schema}.{decoys_table}\n"
            f"SELECT zinc_id, smiles, mw, logp, rotb, hbd, hba, q, '{mol_id}'\n"
            f"  FROM decoys AS D\n"
            f"    ON CONFLICT (zinc_id) DO NOTHING\n")
        _connect.commit()
    _connect.close()


connect = psycopg2.connect(host=args.host,
                           dbname=args.dbname,
                           user=args.user,
                           port=args.port)
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
    f = gzip.open(actives_file, 'r')
else:
    f = open(actives_file, 'r')

if any([
        actives_file.match(i)
        for i in ('*.smi', '*.smi.gz', '*.ism', '*.ism.gz')
]):
    supplier = Chem.SmilesMolSupplierFromText(f.read())
elif any([actives_file.match(i) for i in ('*sdf', '*.sdf.gz')]):
    supplier = Chem.ForwardSDMolSupplier(f)
else:
    raise Exception(f"Can't read molecule from {actives_file}")

from collections import Counter
mol_count = Counter()
for p in map(getProp, supplier):
    # props: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    if p is not None:
        p = list(p)
        mol_id = p[0]
        mol_count[mol_id] += 1
        p[0] = f"{mol_id}_{mol_count[mol_id]}"
        actives_props.append(p)
print(f"{dt.now()}: loading {len(actives_props)} actives from {actives_file}")

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
               f" mol_id text);")
connect.commit()

# get decoys for each active
pool = Pool(processes=args.processes)
for decoys in tqdm(pool.imap_unordered(generate_decoys, actives_props),
                   desc=f"Decoys for {target}",
                   total=len(actives_props),
                   smoothing=0):
    pass
pool.close()

cursor.execute(
    f"SELECT zinc_id, smiles, mw, logp, rotb, hbd, hba, q, mol_id\n"
    f"  FROM (SELECT d.*,\n"
    f"               ROW_NUMBER() OVER (PARTITION BY mol_id) AS r\n"
    f"          FROM {args.schema}.{decoys_table} d) D\n"
    f"WHERE D.r <= 50;")
decoys = cursor.fetchall()

cursor.execute(f"DROP TABLE IF EXISTS {args.schema}.{actives_table};\n"
               f"DROP TABLE IF EXISTS {args.schema}.{actives_table}_fps;\n"
               f"DROP TABLE IF EXISTS {args.schema}.{decoys_table};\n")
connect.close()

print(f"Total number of decoys: {len(decoys)}")
output = Path(args.output)
tdir = output / target
tdir.mkdir(parents=True, exist_ok=True)
actives_file = tdir / "actives_final.smi"
decoys_file = tdir / "decoys_final.smi"
with open(actives_file, 'w') as f:
    for p in actives_props:
        # p: mol_id, smiles, mw, logp, rotb, hbd, hba, q
        f.write(f'{p[1]} {p[0]}\n')
with open(decoys_file, 'w') as f:
    for decoy in decoys:
        smiles = decoy[1]
        zinc_id = decoy[0]
        f.write(f"{smiles} ZINC{zinc_id}\n")

print(f"Total elapsed time: {dt.now()-start}\n")
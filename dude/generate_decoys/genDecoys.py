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
    _cursor.execute(
        f"CREATE TEMP TABLE decoys (like {args.schema}.{args.subset});"
        f"SET enable_seqscan TO OFF;")
    mol_id, smiles, mw, logp, rotb, hbd, hba, q = active
    for i, level in enumerate(PROP_DIFF):
        if i > args.match_level:
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

for p in map(getProp, supplier):
    # props: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    if p is not None:
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
               f" max_tc real);")
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


def assign_decoys(actives_props, decoys, target_dir, heavyMW500=False):
    actives_file = tdir / "actives_final.smi"
    decoys_file = tdir / "decoys_final.smi"
    # prop: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    with open(actives_file, 'w') as f:
        for p in actives_props:
            # p: mol_id, smiles, mw, logp, rotb, hbd, hba, q
            f.write(f'{p[1]} {p[0]}\n')
    decoys_list = [[] for i in actives_props]
    decoys_num = np.zeros(len(actives_props), dtype=int)
    decoys_unfilled = np.ones(len(actives_props), dtype=bool)
    # active: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    actives_props_np = np.array([i[2:8] for i in actives_props])
    for decoy in decoys:
        if heavyMW500:
            m = Chem.MolFromSmiles(decoy[1])
            if Descriptors.HeavyAtomMolWt(m) > 500:
                continue
        for i, level in enumerate(PROP_DIFF):
            if i > args.match_level:
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
            if decoys_num[match] >= 500:
                decoys_unfilled[match] = False
            break

    print(f"decoys info:\n"
          f"#actoves assigned < {args.n} "
          f"decoys: {sum(decoys_num < args.n)}/{len(decoys_num)}")
    for i in np.flatnonzero(decoys_num < args.n):
        print(f"{decoys_num[i]:3d} decoys for active: {actives_props[i]}")

    with open(decoys_file, 'w') as f:
        for ds in decoys_list:
            if len(ds) > args.n:
                ds = random.choices(ds, k=args.n)
            for decoy in ds:
                smiles, mol_id = decoy
                f.write(f"{smiles} ZINC{mol_id}\n")


output = Path(args.output)
tdir = output / target
tdir.mkdir(parents=True, exist_ok=True)
assign_decoys(actives_props, decoys, tdir, heavyMW500=False)

if args.heavyMW500:
    tdir = output.with_suffix('.heavyMW500') / target
    tdir.mkdir(parents=True, exist_ok=True)
    print("Limiting the molecular weight of heavy atoms <= 500")
    assign_decoys(actives_props, decoys, tdir, heavyMW500=True)

print(f"Total elapsed time: {dt.now()-start}\n")
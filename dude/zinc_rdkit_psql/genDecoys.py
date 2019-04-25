"""Generating decoys from database based on active smiles.
"""
import random
import argparse
import subprocess
from pathlib import Path
from datetime import datetime as dt

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.DataStructs import BulkTanimotoSimilarity

import psycopg2
import psycopg2.extras
psycopg2.extensions.set_wait_callback(psycopg2.extras.wait_select)

from tqdm import tqdm
from multiprocessing import Pool

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-a", "--actives", nargs='+', required=True)
parser.add_argument(
    "-tc",
    default=0.35,
    type=float,
    help=
    "max similar (tc, 0-1) between actives and decoys agaist same target, default 0.35"
)
parser.add_argument(
    "-x",
    default=0,
    type=float,
    help=
    "min similar (tc, 0-1) between actives and decoys agaist different targets, default 0.0"
)
parser.add_argument(
    "-X",
    default=1,
    type=float,
    help=
    "max similar (tc, 0-1) between decoys and decoys agaist different targets, default 1"
)
parser.add_argument(
    "-n",
    default=50,
    type=int,
    help=
    "n decoys per acive. default 50"
)
parser.add_argument(
    "-N",
    default=50,
    type=int,
    help=
    "N candidate decoys per acive, random select n decoys from N candidate decoys, N >= n. default: N = n"
)
parser.add_argument("-d", "--dbname", required=True)
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
parser.add_argument("-o", "--output", required=True, help="output dir")
args = parser.parse_args()

if args.N is None:
    args.N = args.n

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


def props_generator_from_files(smiles_files):
    for smi in smiles_files:
        print("{}: loading {}".format(dt.now(), smi))
        smi = Path(smi)
        if smi.suffix == '.gz':
            f = gzip.open(smi, 'rt')
        else:
            f = open(smi, 'r')
        for props in map(getProp, f):
            if props is not None:
                yield props


# It is long but simple, select num decoys based range from narrow to wide

### DEFAULTS FROM MYSINGER - these are used if property not specified
# PROP_DIFF = np.array([
#     [ 20,  35,  50,  65,  80, 100, 125],#"MWT_RANGES"
#     [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6],#"LOGP_RANGES"
#     [  1,   2,   2,   3,   3,   4,   5],#"RB_RANGES"
#     [  0,   0,   1,   1,   2,   2,   3],#"HBD_RANGES"
#     [  0,   1,   2,   2,   3,   3,   4],#"HBA_RANGES"
#     [  0,   0,   0,   0,   0,   1,   2],#"CHG_RANGES"
#     ])

# SET TRANSACTION ISOLATION LEVEL SERIALIZABLE;
decoys_query = """
DROP TABLE IF EXISTS decoys;
CREATE TEMP TABLE decoys (
    zinc_id integer PRIMARY KEY,
    smiles text,
    fp bfp,
    target text);

set rdkit.tanimoto_threshold={tc};

INSERT INTO decoys
SELECT zinc_id, smiles, mfp2, '{target}'
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE ABS (mw - {mw}) <= 20.0
   AND ABS (logp - {logp}) <= 0.4
   AND ABS (rotb - {rotb}) in (0, 1)
   AND hbd = {hbd}
   AND hba = {hba}
   AND q = {q}
   AND NOT mfp2 % morganbv_fp(mol_from_smiles('{smiles}'::cstring))
   AND NOT EXISTS (SELECT 1 FROM {job}_d_fps AS D
                     WHERE D.zid = zinc_id
                        OR tanimoto_sml(D.fp,mfp2) > {X})
   AND NOT EXISTS (SELECT 1 FROM {job}_a_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND (NOT EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}')
         OR EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}'
                       AND tanimoto_sml(fp,mfp2) > {x}))
    --   OR EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}' AND fp%mfp2))
 LIMIT {num} - (SELECT COUNT(*) FROM decoys)
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles, mfp2, '{target}'
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE ABS (mw - {mw}) <= 35.0
   AND ABS (mw - {mw}) > 20.0
   AND ABS (logp - {logp}) <= 0.8
   AND ABS (rotb - {rotb}) in (0, 1, 2)
   AND hbd = {hbd}
   AND ABS (hba - {hba}) in  (0, 1)
   AND q = {q}
   AND NOT mfp2 % morganbv_fp(mol_from_smiles('{smiles}'::cstring))
   AND NOT EXISTS (SELECT 1 FROM {job}_d_fps AS D
                     WHERE D.zid = zinc_id
                        OR tanimoto_sml(D.fp,mfp2) > {X})
   AND NOT EXISTS (SELECT 1 FROM {job}_a_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND (NOT EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}')
         OR EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}'
                       AND tanimoto_sml(fp,mfp2) > {x}))
 LIMIT {num} - (SELECT COUNT(*) FROM decoys)
ON CONFLICT (zinc_id) DO NOTHING;
    
INSERT INTO decoys
SELECT zinc_id, smiles, mfp2, '{target}'
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE ABS (mw - {mw}) <= 50.0
   AND ABS (mw - {mw}) > 35.0
   AND ABS (logp - {logp}) <= 1.2
   AND ABS (rotb - {rotb}) in (0, 1, 2)
   AND ABS (hbd - {hbd}) in (0, 1)
   AND ABS (hba - {hba}) in (0, 1, 2)
   AND q = {q}
   AND NOT mfp2 % morganbv_fp(mol_from_smiles('{smiles}'::cstring))
   AND NOT EXISTS (SELECT 1 FROM {job}_d_fps AS D
                     WHERE D.zid = zinc_id
                        OR tanimoto_sml(D.fp,mfp2) > {X})
   AND NOT EXISTS (SELECT 1 FROM {job}_a_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND (NOT EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}')
         OR EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}'
                       AND tanimoto_sml(fp,mfp2) > {x}))
 LIMIT {num} - (SELECT COUNT(*) FROM decoys)
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles, mfp2, '{target}'
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE ABS (mw - {mw}) <= 65.0
   AND ABS (mw - {mw}) > 50.0
   AND ABS (logp - {logp}) <= 1.8
   AND ABS (rotb - {rotb}) in (0, 1, 2, 3)
   AND ABS (hbd - {hbd}) in (0, 1)
   AND ABS (hba - {hba}) in (0, 1, 2)
   AND q = {q}
   AND NOT mfp2 % morganbv_fp(mol_from_smiles('{smiles}'::cstring))
   AND NOT EXISTS (SELECT 1 FROM {job}_d_fps AS D
                     WHERE D.zid = zinc_id
                        OR tanimoto_sml(D.fp,mfp2) > {X})
   AND NOT EXISTS (SELECT 1 FROM {job}_a_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND (NOT EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}')
         OR EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}'
                       AND tanimoto_sml(fp,mfp2) > {x}))
 LIMIT {num} - (SELECT COUNT(*) FROM decoys)
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles, mfp2, '{target}'
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE ABS (mw - {mw}) <= 80.0
   AND ABS (mw - {mw}) > 65.0
   AND ABS (logp - {logp}) <= 2.4
   AND ABS (rotb - {rotb}) in (0, 1, 2, 3)
   AND ABS (hbd - {hbd}) in (0, 1, 2)
   AND ABS (hba - {hba}) in (0, 1, 2, 3)
   AND q = {q}
   AND NOT mfp2 % morganbv_fp(mol_from_smiles('{smiles}'::cstring))
   AND NOT EXISTS (SELECT 1 FROM {job}_d_fps AS D
                     WHERE D.zid = zinc_id
                        OR tanimoto_sml(D.fp,mfp2) > {X})
   AND NOT EXISTS (SELECT 1 FROM {job}_a_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND (NOT EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}')
         OR EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}'
                       AND tanimoto_sml(fp,mfp2) > {x}))
 LIMIT {num} - (SELECT COUNT(*) FROM decoys)
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles, mfp2, '{target}'
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE ABS (mw - {mw}) <= 100.0
   AND ABS (mw - {mw}) > 80.0
   AND ABS (logp - {logp}) <= 3.0
   AND ABS (rotb - {rotb}) in (0, 1, 2, 3, 4)
   AND ABS (hbd - {hbd}) in (0, 1, 2)
   AND ABS (hba - {hba}) in (0, 1, 2, 3)
   AND ABS (q - {q}) in (0, 1)
   AND NOT mfp2 % morganbv_fp(mol_from_smiles('{smiles}'::cstring))
   AND NOT EXISTS (SELECT 1 FROM {job}_d_fps AS D
                     WHERE D.zid = zinc_id
                        OR tanimoto_sml(D.fp,mfp2) > {X})
   AND NOT EXISTS (SELECT 1 FROM {job}_a_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND (NOT EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}')
         OR EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}'
                       AND tanimoto_sml(fp,mfp2) > {x}))
 LIMIT {num} - (SELECT COUNT(*) FROM decoys)
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles, mfp2, '{target}'
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE ABS (mw - {mw}) <= 125.0
   AND ABS (mw - {mw}) > 100.0
   AND ABS (logp - {logp}) <= 3.6
   AND ABS (rotb - {rotb}) in (0, 1, 2, 3, 4, 5)
   AND ABS (hbd - {hbd}) in (0, 1, 2, 3)
   AND ABS (hba - {hba}) in (0, 1, 2, 3, 4)
   AND ABS (q - {q}) in (0, 1, 2)
   AND NOT mfp2 % morganbv_fp(mol_from_smiles('{smiles}'::cstring))
   AND NOT EXISTS (SELECT 1 FROM {job}_d_fps AS D
                     WHERE D.zid = zinc_id
                        OR tanimoto_sml(D.fp,mfp2) > {X})
   AND NOT EXISTS (SELECT 1 FROM {job}_a_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND (NOT EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}')
         OR EXISTS (SELECT 1 FROM {job}_a_fps WHERE target <> '{target}'
                       AND tanimoto_sml(fp,mfp2) > {x}))
 LIMIT {num} - (SELECT COUNT(*) FROM decoys)
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO {job}_d_fps
SELECT zinc_id, fp, target
  FROM decoys
  LIMIT {num}
ON CONFLICT (zid) DO NOTHING;

-- SELECT smiles, zinc_id, mw, logp, rotb, hbd, hba, q
--   FROM decoys JOIN dud.props using (zinc_id)
--  LIMIT {num};

SELECT smiles, zinc_id FROM decoys LIMIT {num};

"""
# AND zinc_id NOT IN (SELECT zinc_id FORM decoys)
start = dt.now()

# loading actives
targets = []
actives_props = []
for a_file in args.actives:
    a_file = Path(a_file)
    if len(a_file.parts) > 1:
        target = a_file.parts[-2]
    else:
        target = a_file.stem
    targets.append(target)
    props = []
    total_mol = 0
    if a_file.suffix == '.gz':
        f = gzip.open(a_file, 'rt')
    else:
        f = open(a_file, 'r')
    for p in map(getProp, f):
        # props: mol_id, smiles, mw, logp, rotb, hbd, hba, q
        if p is not None:
            props.append(p)
    actives_props.append(props)
    print("{} loading {:7s} {:4d} actives from {}".format(
        dt.now() - start, target, len(props), a_file))

connect = psycopg2.connect(host=args.host, dbname=args.dbname, port=args.port)
cursor = connect.cursor()

job_name = 'job' + str(
    abs(hash(''.join(args.actives) + str(dt.now()))))
# init tables for saving actives
init_actives = """
    CREATE INDEX IF NOT EXISTS prop_idx ON dud.props (mw, logp, rotb, hba, hbd, q);

    DROP TABLE IF EXISTS {job}_a_fps;
    DROP TABLE IF EXISTS {job}_a;
    DROP TABLE IF EXISTS {job}_d_fps;

    CREATE TABLE {job}_a (
        name text,
        smiles text,
        target text);
    CREATE TABLE {job}_a_fps (
        name text,
        fp bfp,
        target text);
    CREATE TABLE {job}_d_fps (
        zid int PRIMARY KEY,
        fp bfp,
        target text);

    CREATE INDEX IF NOT EXISTS dfps_idx ON {job}_d_fps USING gist(fp);
    CREATE INDEX IF NOT EXISTS afps_idx ON {job}_a_fps USING gist(fp);
    CREATE INDEX IF NOT EXISTS d_target_idx ON {job}_d_fps (target);
    CREATE INDEX IF NOT EXISTS a_target_idx ON {job}_a_fps (target);
    """.format(job=job_name)
cursor.execute(init_actives)
connect.commit()

# generate fps into table actives_fps
for i, target in enumerate(targets):
    # props: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    insert_query = 'INSERT INTO {job}_a VALUES %s'.format(
        job=job_name)
    actives_values = [(p[0], p[1], target) for p in actives_props[i]]
    psycopg2.extras.execute_values(cursor,
                                   insert_query,
                                   actives_values,
                                   template=None,
                                   page_size=100)
    cursor.execute("""
        INSERT INTO {job}_a_fps
        SELECT name, morganbv_fp(mol_from_smiles(smiles::cstring)) as fp, target
        FROM {job}_a;
        """.format(job=job_name))
    connect.commit()

def generate_decoys(kwargs):
    _t = dt.now()
    _connect = psycopg2.connect(host=args.host,
                                dbname=args.dbname,
                                port=args.port)
    _cursor = _connect.cursor()
    _cursor.execute(decoys_query.format(**kwargs))
    _decoys = _cursor.fetchall()
    _connect.commit()
    _connect.close()
    print("Time for query {} decoys for 1 active: {}".format(len(_decoys), dt.now() - _t))
    return _decoys


# get decoys for each actives
output = Path(args.output)
output.mkdir(parents=True, exist_ok=True)
for i, target in enumerate(targets):
    tdir = output / target
    tdir.mkdir(exist_ok=True)
    a_file = tdir / "actives_final.smi"
    d_file = tdir / "decoys_final.smi"
    # prop: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    d_file = tdir / "decoys_final.smi"
    with open(a_file, 'w') as f:
        for p in actives_props[i]:
            # p: mol_id, smiles, mw, logp, rotb, hbd, hba, q
            f.write('{} {}\n'.format(p[1], p[0]))
    f = open(d_file, 'w')
    job_kwargs = []
    for p in actives_props[i]:
        mol_id, smiles, mw, logp, rotb, hbd, hba, q = p
        job_kwargs.append({
            'target': target,
            'job': job_name,
            'num': args.N,
            'mw': mw,
            'logp': logp,
            'rotb': rotb,
            'hbd': hbd,
            'hba': hba,
            'q': q,
            'tc': args.tc,
            'x': args.x,
            'X': args.X,
            'smiles':smiles,
        })
    pool = Pool()
    N = len(job_kwargs)
    print("generating decoys for {} actives against target {}:".format(
        N, target))
    for decoys in tqdm(pool.imap_unordered(generate_decoys, job_kwargs),
                       total=N):
        if len(decoys) > args.n:
            decoys = random.sample(decoys, args.n)
        for smiles, name in decoys:
            f.write('{} {}\n'.format(smiles, name))
    pool.close()
    f.close()

rm_actives = """
    DROP TABLE IF EXISTS {job}_a_fps;
    DROP TABLE IF EXISTS {job}_d_fps;
    DROP TABLE IF EXISTS {job}_a;
    """.format(job=job_name)
cursor.execute(rm_actives)
connect.commit()
connect.close()
# subprocess.call(['pg_ctl', '-o', port_option, '-D', args.dbpath, '-l', 'logfile', 'stop'])
# print("you can connect to db using:\n" +
# "pg_ctl -o '-F -p {}' -D {} -l logfile start\n".format(args.port, args.dbpath) +
# "psql -p {} {}\n".format(args.port, args.dbname))
print("Total elapsed time: {}".format(dt.now() - start))

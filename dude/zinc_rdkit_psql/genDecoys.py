"""Generating decoys from database based on active smiles.
"""
import json
import random
import argparse
import numpy as np
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
    default=35,
    type=int,
    help="max similar between actives and decoys agaist same target, default 35"
)
parser.add_argument(
    "-x",
    default=35,
    type=int,
    help=
    "min similar between actives and decoys agaist different targets, default 35"
)
parser.add_argument("-n",
                    default=50,
                    type=int,
                    help="n decoys per acive. default 50")
parser.add_argument("-N",
                    type=int,
                    help="N candidate decoys per acive, default N = n")
parser.add_argument("-d", "--dbname", required=True)
parser.add_argument("-s", "--schema", required=True)
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
parser.add_argument("-P", "--processes", type=int)
parser.add_argument("--remove_old_simi", action="store_true")
parser.add_argument("-m", "--match_level", type=int, default=3)
parser.add_argument("-o", "--output", required=True, help="output dir")
args = parser.parse_args()

start = dt.now()

if args.N is None:
    args.N = args.n

connect = psycopg2.connect(host=args.host, dbname=args.dbname, port=args.port)
cursor = connect.cursor()


# loading actives
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

job_name = 'job' + str(abs(hash(''.join(args.actives) + str(dt.now()))))
# init tables for saving actives
init_actives = """
CREATE TABLE {job}_a (
    name text,
    smiles text,
    target text);
""".format(job=job_name, schema=args.schema)
cursor.execute(init_actives)
connect.commit()

# generate fps into table actives_fps
for i, target in enumerate(targets):
    # props: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    insert_query = 'INSERT INTO {job}_a VALUES %s'.format(job=job_name)
    actives_values = [(p[0], p[1], target) for p in actives_props[i]]
    psycopg2.extras.execute_values(cursor,
                                   insert_query,
                                   actives_values,
                                   template=None,
                                   page_size=100)
    cursor.execute("""
        SELECT morganbv_fp(mol_from_smiles(smiles::cstring)) as fp
        INTO {job}_a_fps_{target}
        FROM {job}_a
        WHERE target = '{target}';

        CREATE INDEX IF NOT EXISTS afps_idx_{target}
        ON {job}_a_fps_{target} USING gist(fp);
        """.format(job=job_name, target=target))
    connect.commit()

### get a subset mols which are similar to actives with tanimoto threshold.
simi_kwargs = []
for i, target in enumerate(targets):
    kwargs = {'schema': args.schema, 'simi': args.tc, 'target': target}
    simi_kwargs.append(kwargs)
    if args.x > 0 and args.x != args.tc:
        kwargs = kwargs.copy()
        kwargs['simi'] = args.x
        simi_kwargs.append(kwargs)

if args.remove_old_simi:
    OLD_SIMI_TABLES = []
else:
    cursor.execute("""
    SELECT table_name
      FROM information_schema.tables
     WHERE table_name ~ 'simi'
       AND table_schema = '{schema}'
    """.format(schema=args.schema))
    OLD_SIMI_TABLES = set([i[0] for i in cursor.fetchall()])

cursor.execute("""
SELECT table_name
  FROM information_schema.tables
 WHERE table_name ~ 'part'
   AND table_schema = '{schema}'
""".format(schema=args.schema))
part_tables = [i[0] for i in cursor.fetchall()]

create_table_simi = """
DROP TABLE IF EXISTS {schema}.simi{simi}_{target};
CREATE TABLE IF NOT EXISTS {schema}.simi{simi}_{target} (
            zinc_id integer,
            smiles text,
            mw real,
            logp real,
            rotb smallint,
            hbd smallint,
            hba smallint,
            q smallint)
"""
simi_part_kwargs = []
for kwargs in simi_kwargs:
    table = 'simi{simi}_{target}'.format(**kwargs)
    if table in OLD_SIMI_TABLES:
        continue
    cursor.execute(create_table_simi.format(**kwargs))
    for part in part_tables:
        d = kwargs.copy()
        d['part'] = part
        d['simi_float'] = kwargs['simi'] / 100.
        d['job'] = job_name
        simi_part_kwargs.append(d)
connect.commit()

get_simi_decoys_query = """
CREATE INDEX IF NOT EXISTS fps_idx_{part} ON {schema}.{part} USING gist(mfp2);
SET rdkit.tanimoto_threshold = {simi_float};

INSERT INTO {schema}.simi{simi}_{target}
SELECT DISTINCT ON (zinc_id) zinc_id, smiles, mw, logp, rotb, hbd, hba, q
  FROM {schema}.{part} db
 CROSS JOIN {job}_a_fps_{target} a
 WHERE a.fp % db.mfp2;
"""


def get_simi(kwargs):
    _connect = psycopg2.connect(host=args.host,
                                dbname=args.dbname,
                                port=args.port)
    _cursor = _connect.cursor()
    _t = dt.now()
    # print(get_simi_decoys_query.format(**kwargs))
    _cursor.execute(get_simi_decoys_query.format(**kwargs))
    print("Time for get simi {} for {:>10s} from {:>8s}: {}".format(
        kwargs['simi'], kwargs['target'], kwargs['part'],
        dt.now() - _t))
    _connect.commit()
    _connect.close()


if OLD_SIMI_TABLES:
    print('\nUsing old simi tables\n')
N = len(simi_part_kwargs)
# mix different targets for better tqdm timing.
simi_part_kwargs = sorted(simi_part_kwargs, key=lambda i: i['part'])
if N > 0:
    print('get {} more simi tables'.format(N))
    pool = Pool(processes=args.processes)
    for _ in tqdm(pool.imap_unordered(get_simi, simi_part_kwargs),
                  total=N,
                  smoothing=0):
        pass
    pool.close()

create_index_on_zinc_id = """
CREATE INDEX IF NOT EXISTS simi{simi}_{target}_pkey
ON {schema}.simi{simi}_{target} (zinc_id);
"""
for kwargs in simi_kwargs:
    cursor.execute(create_index_on_zinc_id.format(**kwargs))
connect.commit()

rm_actives = """
    DROP TABLE IF EXISTS {job}_a;
    DROP TABLE IF EXISTS {job}_a_fps;
    """.format(job=job_name)
cursor.execute(rm_actives)
connect.commit()

# It is long but simple, select num decoys based range from narrow to wide

### DEFAULTS FROM MYSINGER - these are used if property not specified
PROP_DIFF = np.array([
    [-1, 20, 35, 50, 65, 80, 100],  #"MWT_RANGES_LOWER"
    [20, 35, 50, 65, 80, 100, 125],  #"MWT_RANGES_UPPER"
    [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6],  #"LOGP_RANGES"
    [1, 2, 2, 3, 3, 4, 5],  #"RB_RANGES"
    [0, 0, 1, 1, 2, 2, 3],  #"HBD_RANGES"
    [0, 1, 2, 2, 3, 3, 4],  #"HBA_RANGES"
    [0, 0, 0, 0, 0, 1, 2],  #"CHG_RANGES"
])

props_match_block = """ABS (mw - {{mw}}) > {d_mw_l}
AND ABS (mw - {{mw}}) <= {d_mw_u}
AND ABS (logp - {{logp}}) <= {d_logp}
AND ABS (rotb - {{rotb}}) <= {d_rotb}
AND ABS (hbd - {{hbd}}) <= {d_hbd}
AND ABS (hba - {{hba}}) <= {d_hba}
AND ABS (q - {{q}}) <= {d_q}
"""
diff_keys = ('d_mw_l', 'd_mw_u', 'd_logp', 'd_rotb', 'd_hbd', 'd_hba', 'd_q')
MATCH_LEVEL_BLOCK = []
for i, level in enumerate(PROP_DIFF.T):
    # '{{}}'.format() ==> '{}'
    block = props_match_block.format(**dict(zip(diff_keys, level)))
    MATCH_LEVEL_BLOCK.append(block)


def generate_decoys(kwargs):
    _t = dt.now()
    _connect = psycopg2.connect(host=args.host,
                                dbname=args.dbname,
                                port=args.port)
    _cursor = _connect.cursor()
    _cursor.execute("""CREATE TEMP TABLE decoys (
    zinc_id integer PRIMARY KEY,
    smiles text);
    SET enable_seqscan TO OFF;
    """)
    for level, match in enumerate(MATCH_LEVEL_BLOCK):
        query = 'INSERT INTO decoys\n'
        select_subquerys = []
        for ti in kwargs['targets']:
            if ti == kwargs['target']: continue
            q = (
                '(SELECT zinc_id, smiles\n' +
                'FROM {{schema}}.simi{{x}}_{ti} other\n'.format(ti=ti) +
                'WHERE NOT EXISTS (\n' +
                '    SELECT FROM {schema}.simi{tc}_{target} self WHERE self.zinc_id = other.zinc_id)\n'
                'AND ' + match + 'ORDER BY RANDOM() LIMIT {num})\n')
            select_subquerys.append(q)
        query += 'UNION ALL\n'.join(select_subquerys) + '\n'
        query += 'ON CONFLICT (zinc_id) DO NOTHING;\n'
        query += 'SELECT smiles, zinc_id FROM decoys ORDER BY RANDOM() LIMIT {num};'
        # print(query)
        # print(query.format(**kwargs))
        _cursor.execute(query.format(**kwargs))
        _decoys = _cursor.fetchall()
        if len(_decoys) >= kwargs['num'] or level >= args.match_level:
            break
    _connect.commit()
    _connect.close()
    # print("Time for query {:3d} decoys for 1 active in {:>10s}: {}".format(
    #     len(_decoys), kwargs['target'],
    #     dt.now() - _t))
    return _decoys


# get decoys for each actives
output = Path(args.output)
output.mkdir(parents=True, exist_ok=True)
with open(output / 'args.json', 'w') as f:
    json.dump(vars(args), f, sort_keys=True, indent=4)
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
            'targets': targets,
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
            'smiles': smiles,
            'schema': args.schema,
        })
    pool = Pool(processes=args.processes)
    N = len(job_kwargs)
    print("generating decoys for {} actives against target {}:".format(
        N, target))
    for decoys in tqdm(pool.imap_unordered(generate_decoys, job_kwargs),
                       total=N,
                       smoothing=0):
        if len(decoys) > args.n:
            decoys = random.sample(decoys, args.n)
        for smiles, name in decoys:
            f.write('{} {}\n'.format(smiles, name))
    pool.close()
    f.close()

connect.close()

print("Total elapsed time: {}".format(dt.now() - start))

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
parser.add_argument("-n",
                    default=50,
                    type=int,
                    help="n decoys per acive. default 50")
parser.add_argument("-d", "--dbname", required=True)
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
parser.add_argument("-P", "--processes", type=int)
parser.add_argument("-m", "--match_level", type=int, default=2)
parser.add_argument("-o", "--output", required=True, help="output dir")
args = parser.parse_args()

start = dt.now()

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

# init tables for saving actives
init_actives = """
DROP TABLE IF EXISTS actives;
CREATE TABLE actives (
    name text,
    smiles text,
    target text);
"""
cursor.execute(init_actives)
connect.commit()

# generate fps into table actives_fps
for i, target in enumerate(targets):
    # props: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    insert_query = 'INSERT INTO actives VALUES %s'
    actives_values = [(p[0], p[1], target) for p in actives_props[i]]
    psycopg2.extras.execute_values(cursor,
                                   insert_query,
                                   actives_values,
                                   template=None,
                                   page_size=100)
    cursor.execute("""
        DROP TABLE IF EXISTS afps_{target};
        SELECT morganbv_fp(mol_from_smiles(smiles::cstring)) as fp
        INTO afps_{target}
        FROM actives
        WHERE target = '{target}';
        CREATE INDEX IF NOT EXISTS afps_idx_{target}
        ON afps_{target} USING gist(fp);
        """.format(target=target))
    connect.commit()

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
PROP_DIFF = PROP_DIFF.T

props_match_block = """(
    ({{mw}} > mw + {d_mw_l} AND {{mw}} <= mw + {d_mw_u})
    OR
    ({{mw}} < mw - {d_mw_l} AND {{mw}} >= mw - {d_mw_u})
)
AND ABS (logp - {{logp}}) <= {d_logp}
AND ABS (rotb - {{rotb}}) <= {d_rotb}
AND ABS (hbd - {{hbd}}) <= {d_hbd}
AND ABS (hba - {{hba}}) <= {d_hba}
AND ABS (q - {{q}}) <= {d_q}
"""
diff_keys = ('d_mw_l', 'd_mw_u', 'd_logp', 'd_rotb', 'd_hbd', 'd_hba', 'd_q')
MATCH_LEVEL_BLOCK = []
for i, diff_range in enumerate(PROP_DIFF):
    # '{{}}'.format() ==> '{}'
    block = props_match_block.format(**dict(zip(diff_keys, diff_range)))
    MATCH_LEVEL_BLOCK.append(block)

cursor.execute("""
SELECT table_name
  FROM information_schema.tables
 WHERE table_name ~ 'zinc_split1m'
   AND table_schema = 'public';
""")
part_names = [i[0] for i in cursor.fetchall()]

PART_MW_RANGES = {}
for part in part_names:
    cursor.execute("SELECT MIN(mw), MAX(mw) FROM {}".format(part))
    mw_min, mw_max = cursor.fetchone()
    PART_MW_RANGES[part] = (mw_min, mw_max)


def range_overlop(range1, range2):
    range1 = tuple(range1)
    range2 = tuple(range2)
    sum_range = (range1[1] - range1[0]) + (range2[1] - range2[0])
    union_range = max(range1 + range2) - min(range1 + range2)
    overlap = sum_range - union_range
    if overlap < 0:
        overlap = 0
    return overlap


def sample_parts(mw, level):
    mw_diff_low, mw_diff_up = PROP_DIFF[level][:2]
    if mw_diff_low < 0:
        mw_diff_low = 0
    candidate_ranges = [(mw - mw_diff_up, mw - mw_diff_low),
                        (mw + mw_diff_low, mw + mw_diff_up)]
    candidate_parts = []
    weights = []
    for name, part_mw_range in PART_MW_RANGES.items():
        overlap = 0
        for _range in candidate_ranges:
            overlap += range_overlop(_range, part_mw_range)
        if overlap > 0:
            candidate_parts.append(name)
            # weights.append(overlap)
            part_mid = sum(part_mw_range) / 2
            mid_dist2 = min(1, abs(mw - part_mid))**2
            weights.append(1 / mid_dist2)
    weights = np.array(weights) / sum(weights)
    N = len(candidate_parts)
    return np.random.choice(candidate_parts, size=N, replace=False, p=weights)


def generate_decoys(kwargs):
    _t = dt.now()
    _connect = psycopg2.connect(host=args.host,
                                dbname=args.dbname,
                                port=args.port)
    _cursor = _connect.cursor()
    _decoys = []
    num = args.n * 15
    for level, match in enumerate(MATCH_LEVEL_BLOCK):
        if level > args.match_level:
            continue
        selected_parts = sample_parts(kwargs['mw'], level)
        for part in selected_parts:
            query = '''
                SET enable_seqscan TO OFF;
                SET rdkit.tanimoto_threshold = {tc};
                SELECT smiles, zinc_id
                FROM {part}
                WHERE NOT EXISTS (
                    SELECT FROM afps_{target} a
                    WHERE a.fp % mfp2)
                AND {match}
                -- ORDER BY db.mfp2 <%> morganbv_fp(mol_from_smiles('{smiles}'::cstring)) DESC
                LIMIT {num};
                '''.format(match=match,
                           tc=args.tc / 100.,
                           part=part,
                           num=15 * args.n - len(_decoys),
                           **kwargs)
            # print(query.format(**kwargs))
            _cursor.execute(query.format(**kwargs))
            _decoys.extend(_cursor.fetchall())
            _connect.commit()
            # print(part, len(_decoys), PART_MW_RANGES[part], kwargs['mw'])
            if len(_decoys) >= args.n * 15:
                _connect.close()
                return _decoys
    _connect.close()
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
    active_kwargs = []
    for active in actives_props[i]:
        mol_id, smiles, mw, logp, rotb, hbd, hba, q = active
        active_kwargs.append({
            'target': target,
            'mw': mw,
            'logp': logp,
            'rotb': rotb,
            'hbd': hbd,
            'hba': hba,
            'q': q,
            'smiles': smiles,
        })
    pool = Pool(processes=args.processes)
    N = len(active_kwargs)
    print("generating decoys for {} actives against target {}:".format(
        N, target))
    # for decoys in tqdm(map(generate_decoys, active_kwargs)): # for debug
    for decoys in tqdm(pool.imap_unordered(generate_decoys, active_kwargs),
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

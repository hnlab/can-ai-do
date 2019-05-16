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
parser.add_argument("-s", "--schema", required=True)
parser.add_argument("-S", "--schema_part")
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
parser.add_argument("-P", "--processes", type=int)
parser.add_argument("--remove_old_simi", action="store_true")
parser.add_argument("-m", "--match_level", type=int, default=3)
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

cursor.execute("""
SELECT table_name
  FROM information_schema.tables
 WHERE table_name ~ 'part'
   AND table_schema = '{schema}'
""".format(schema=args.schema))
part_tables = [i[0] for i in cursor.fetchall()]

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
    zinc_id integer,
    smiles text);
    SET enable_seqscan TO OFF;
    """)
    np.random.seed()
    for level, match in enumerate(MATCH_LEVEL_BLOCK):
        if level >= args.match_level:
            continue
        for part in np.random.permutation(part_tables):
            query = '''
                SET rdkit.tanimoto_threshold = {tc};
                INSERT INTO decoys
                SELECT zinc_id, smiles
                FROM {schema}.{part} db
                WHERE NOT EXISTS (
                    SELECT FROM {schema}.{part} ex
                    WHERE ex.mfp2 % morganbv_fp(mol_from_smiles('{smiles}'::cstring))
                    AND db.zinc_id = ex.zinc_id
                )
                AND {match}
                ORDER BY db.mfp2 <%> morganbv_fp(mol_from_smiles('{smiles}'::cstring)) DESC
                LIMIT 15 * {num} - (SELECT COUNT(*) FROM decoys);
                
                SELECT smiles, zinc_id FROM decoys ORDER BY RANDOM() LIMIT {num};
                '''.format(match=match,
                           tc=args.tc / 100.,
                           part=part,
                           schema=args.schema,
                           num=args.n,
                           **kwargs)
            # print(query.format(**kwargs))
            _cursor.execute(query.format(**kwargs))
            _decoys = _cursor.fetchall()
            if len(_decoys) >= args.n:
                _connect.commit()
                _connect.close()
                return _decoys
    _connect.commit()
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

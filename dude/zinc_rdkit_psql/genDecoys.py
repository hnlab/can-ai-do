"""Load mols from smiles into database.
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

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-a", "--actives", nargs='+', required=True)
parser.add_argument("-D", "--dbpath", required=True)
parser.add_argument("-d", "--dbname", required=True)
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
parser.add_argument(
    "-o", "--output", required=True, help="output dir")
args = parser.parse_args()

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
    return tuple([mol_id,smiles, mw, logp, rotb, hbd, hba, q])

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
    smiles text);

set rdkit.tanimoto_threshold=0.35;

INSERT INTO decoys
SELECT zinc_id, smiles
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE (SELECT COUNT(*) FROM decoys) <= {num}
   AND ABS (mw - {mw}) <= 20.0 
   AND ABS (logp - {logp}) <= 0.4
   AND ABS (rotb - {rotb}) in (0, 1)
   AND hbd = {hbd}
   AND hba = {hba}
   AND q = {q}
   AND NOT EXISTS (SELECT 1 FROM actives_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND EXISTS (SELECT 1 FROM actives_fps
                WHERE fp%mfp2 AND target <> '{target}')
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE (SELECT COUNT(*) FROM decoys) <= {num}
   AND ABS (mw - {mw}) <= 35.0
   AND ABS (mw - {mw}) > 20.0
   AND ABS (logp - {logp}) <= 0.8
   AND ABS (rotb - {rotb}) in (0, 1, 2)
   AND hbd = {hbd}
   AND ABS (hba - {hba}) in  (0, 1)
   AND q = {q}
   AND NOT EXISTS (SELECT 1 FROM actives_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND EXISTS (SELECT 1 FROM actives_fps
                WHERE fp%mfp2 AND target <> '{target}')
ON CONFLICT (zinc_id) DO NOTHING;
    
INSERT INTO decoys
SELECT zinc_id, smiles
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE (SELECT COUNT(*) FROM decoys) <= {num}
   AND ABS (mw - {mw}) <= 50.0
   AND ABS (mw - {mw}) > 35.0
   AND ABS (logp - {logp}) <= 1.2
   AND ABS (rotb - {rotb}) in (0, 1, 2)
   AND ABS (hbd - {hbd}) in (0, 1)
   AND ABS (hba - {hba}) in (0, 1, 2)
   AND q = {q}
   AND NOT EXISTS (SELECT 1 FROM actives_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND EXISTS (SELECT 1 FROM actives_fps
                WHERE fp%mfp2 AND target <> '{target}')
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE (SELECT COUNT(*) FROM decoys) <= {num}
   AND ABS (mw - {mw}) <= 65.0
   AND ABS (mw - {mw}) > 50.0
   AND ABS (logp - {logp}) <= 1.8
   AND ABS (rotb - {rotb}) in (0, 1, 2, 3)
   AND ABS (hbd - {hbd}) in (0, 1)
   AND ABS (hba - {hba}) in (0, 1, 2)
   AND q = {q}
   AND NOT EXISTS (SELECT 1 FROM actives_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND EXISTS (SELECT 1 FROM actives_fps
                WHERE fp%mfp2 AND target <> '{target}')
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE (SELECT COUNT(*) FROM decoys) <= {num}
   AND ABS (mw - {mw}) <= 80.0
   AND ABS (mw - {mw}) > 65.0
   AND ABS (logp - {logp}) <= 2.4
   AND ABS (rotb - {rotb}) in (0, 1, 2, 3)
   AND ABS (hbd - {hbd}) in (0, 1, 2)
   AND ABS (hba - {hba}) in (0, 1, 2, 3)
   AND q = {q}
   AND NOT EXISTS (SELECT 1 FROM actives_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND EXISTS (SELECT 1 FROM actives_fps
                WHERE fp%mfp2 AND target <> '{target}')
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE (SELECT COUNT(*) FROM decoys) <= {num}
   AND ABS (mw - {mw}) <= 100.0
   AND ABS (mw - {mw}) > 80.0
   AND ABS (logp - {logp}) <= 3.0
   AND ABS (rotb - {rotb}) in (0, 1, 2, 3, 4)
   AND ABS (hbd - {hbd}) in (0, 1, 2)
   AND ABS (hba - {hba}) in (0, 1, 2, 3)
   AND ABS (q - {q}) in (0, 1)
   AND NOT EXISTS (SELECT 1 FROM actives_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND EXISTS (SELECT 1 FROM actives_fps
                WHERE fp%mfp2 AND target <> '{target}')
ON CONFLICT (zinc_id) DO NOTHING;

INSERT INTO decoys
SELECT zinc_id, smiles
  FROM dud.props JOIN dud.fps using (zinc_id)
 WHERE (SELECT COUNT(*) FROM decoys) <= {num}
   AND ABS (mw - {mw}) <= 125.0
   AND ABS (mw - {mw}) > 100.0
   AND ABS (logp - {logp}) <= 3.6
   AND ABS (rotb - {rotb}) in (0, 1, 2, 3, 4, 5)
   AND ABS (hbd - {hbd}) in (0, 1, 2, 3)
   AND ABS (hba - {hba}) in (0, 1, 2, 3, 4)
   AND ABS (q - {q}) in (0, 1, 2)
   AND NOT EXISTS (SELECT 1 FROM actives_fps
                    WHERE fp%mfp2 AND target = '{target}')
   AND EXISTS (SELECT 1 FROM actives_fps
                WHERE fp%mfp2 AND target <> '{target}')
ON CONFLICT (zinc_id) DO NOTHING;

SELECT smiles, zinc_id FROM decoys;

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

# init a temp postdatabase
subprocess.call(['initdb', '-D', args.dbpath])
port_option = "-F -p {}".format(args.port)
code = subprocess.call(['pg_ctl', '-o', port_option, '-D', args.dbpath,'-l','logfile', 'start'])
# sleep 3 seconds to wait pg_ctl start succeed.
import time
time.sleep(5)
connect = psycopg2.connect(host=args.host, dbname=args.dbname, port=args.port)
cursor = connect.cursor()


# init tables for saving actives
init_actives = """
    DROP TABLE IF EXISTS actives_fps;
    DROP TABLE IF EXISTS actives;
    CREATE TABLE actives (
        name text PRIMARY KEY,
        smiles text,
        target text);
    CREATE TABLE actives_fps (
        name text PRIMARY KEY,
        fp bfp,
        target text);
    """
cursor.execute(init_actives)
# connect.commit()

# generate fps into table actives_fps
for i, target in enumerate(targets):
    # props: mol_id, smiles, mw, logp, rotb, hbd, hba, q
    insert_query = 'INSERT INTO actives VALUES %s ON CONFLICT (name) DO NOTHING'
    actives_values = [(p[0], p[1], target) for p in actives_props[i]]
    psycopg2.extras.execute_values(
        cursor, insert_query, actives_values, template=None, page_size=100)
    cursor.execute("""
        INSERT INTO actives_fps
        SELECT name, morganbv_fp(mol_from_smiles(smiles::cstring)) as fp, target
        FROM actives;
        CREATE INDEX IF NOT EXISTS afps_idx ON actives_fps USING gist(fp);
        """)
    # connect.commit()

t_idx = dt.now()
create_index = """
DROP INDEX IF EXISTS mw_logp_idx;
CREATE INDEX IF NOT EXISTS prop_idx ON dud.props (mw, logp, rotb, hba, hbd, q);
"""
# CREATE INDEX IF NOT EXISTS actives_fps_idx ON actives_fps USING gist(fp);
cursor.execute(create_index)
connect.commit()
print("Time for create index mw_logp_idx: {}".format(dt.now() - t_idx))
# get decoys for each actives
output = Path(args.output)
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
    for p in actives_props[i]:
        mol_id, smiles, mw, logp, rotb, hbd, hba, q = p
        t_q = dt.now()
        cursor.execute(decoys_query.format(
            target=target, num=750, mw=mw, logp=logp, rotb=rotb, hbd=hbd, hba=hba, q=q))
        print("Time for query decoys for {}: {}".format(mol_id, dt.now() - t_q))
        decoys = cursor.fetchall()
        print(decoys)
        for smiles, name in random.sample(decoys, min(50,len(decoys))):
            f.write('{} {}\n'.format(smiles, name))
        connect.commit()
    f.close()

rm_actives = """
    DROP TABLE IF EXISTS actives_fps;
    DROP TABLE IF EXISTS actives;
    """
cursor.execute(rm_actives)
connect.commit()
connect.close()
subprocess.call(['pg_ctl', '-o', port_option, '-D', args.dbpath, '-l', 'logfile', 'stop'])
print("you can connect to db using:\n" +
"pg_ctl -o '-F -p {}' -D {} -l logfile start\n".format(args.port, args.dbpath) +
"psql -p {} {}\n".format(args.port, args.dbname))
print("Total elapsed time: {}".format(dt.now()-start))

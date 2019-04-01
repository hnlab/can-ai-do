"""Load mols from smiles into database.
"""
import gzip
import argparse
from pathlib import Path
from datetime import datetime as dt

import psycopg2

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-s", "--smiles", nargs='+', required=True)
parser.add_argument("-d", "--dbname", required=True)
parser.add_argument("-u", "--users", required=True)
args = parser.parse_args()

start = dt.now()
print("Job start at {}\n".format(start))
conn = psycopg2.connect("host=localhost dbname={} user={}".format(
    args.dbname, args.users))
cur = conn.cursor()
cur.execute("""
    CREATE TABLE IF NOT EXISTS props (
        id serial,
        smiles text NOT NULL,
        zinc_id integer,
        mw real,
        logp real,
        rotb smallint,
        hbd smallint,
        hba smallint,
        q smallint)
    """)
for smi in args.smiles:
    smi = Path(smi)
    if smi.suffix == '.gz':
        f = gzip.open(smi, 'rt')
    else:
        f = open(smi, 'r')
    header = next(f)
    cur.copy_from(
        f,
        'props',
        columns=('smiles', 'zinc_id', 'mw', 'logp', 'rotb', 'hbd', 'hba', 'q'),
        sep=' ')
    f.close()
    conn.commit()
    print("{}: loaded {}".format(dt.now()-start, smi))

# generate DUD schema from https://www.rdkit.org/docs/Cartridge.html#loading-chembl
sql = """
CREATE EXTENSION IF NOT EXISTS rdkit;
CREATE SCHEMA dud;

SELECT * INTO dud.mols
FROM (SELECT id, mol_from_smiles(smiles::cstring) m  FROM props) tmp 
WHERE m IS NOT NULL;

CREATE INDEX molidx ON dud.mols USING gist(m);
ALTER TABLE dud.mols ADD PRIMARY KEY (id);

SELECT id, morganbv_fp(m) as mfp2 INTO dud.fps FROM dud.mols;
CREATE INDEX fps_mfp2_idx ON dud.fps USING gist(mfp2);
ALTER TABLE dud.fps ADD PRIMARY KEY (id);

CREATE OR REPLACE FUNCTION get_mfp2_neighbors(smiles text)
RETURNS table (id integer, m mol, similarity double precision) AS $$
    SELECT id, m, tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)), mfp2) AS similarity
    FROM dud.fps JOIN dud.mols USING (id)
    WHERE morganbv_fp(mol_from_smiles($1::cstring))%mfp2
    ORDER BY morganbv_fp(mol_from_smiles($1::cstring))<%>mfp2;
$$ language sql stable;
"""
start_dud = dt.now()
print("\nStart create schema dud at {}".format(start_dud))
cur.execute(sql)
conn.commit()
conn.close()
print("schema dud create in {}\n".format(dt.now()-start_dud))
print("Job end at {}".format(dt.now()))
print("Total elapsed time: {}".format(dt.now()-start))

"""Combine sql.gz into database.
"""
import time
import shutil
import argparse
import subprocess
from pathlib import Path
from datetime import datetime as dt

import psycopg2

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-s", "--sql", nargs='+', required=True)
parser.add_argument("-D", "--dbpath")
parser.add_argument("-d", "--dbname", required=True, help="zinc12 or zinc15")
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
args = parser.parse_args()

start = dt.now()

if args.dbpath is not None:
    # init a temp postdatabase
    subprocess.call(['initdb', '-D', args.dbpath])
    port_option = "-F -p {}".format(args.port)
    code = subprocess.call([
        'pg_ctl', '-o', port_option, '-D', args.dbpath, '-l', 'logfile',
        'start'
    ])
    # sleep 3 seconds to wait pg_ctl start succeed.
    time.sleep(5)
    subprocess.call(['createdb', '-p', args.port, args.dbname])
    time.sleep(3)

connect = psycopg2.connect(host=args.host, dbname=args.dbname, port=args.port)
cursor = connect.cursor()

# init a schema dud
init_dud = """
    CREATE EXTENSION IF NOT EXISTS rdkit;
    CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA public;
    CREATE SCHEMA IF NOT EXISTS dud;
    CREATE TABLE IF NOT EXISTS dud.props (
        zinc_id integer PRIMARY KEY,
        smiles text,
        mw real,
        logp real,
        rotb smallint,
        hbd smallint,
        hba smallint,
        q smallint);
    CREATE TABLE IF NOT EXISTS dud.fps (
        zinc_id integer PRIMARY KEY,
        mfp2 bfp);
    CREATE TABLE IF NOT EXISTS dud.mols (
        zinc_id integer PRIMARY KEY,
        m mol);
    DROP TABLE IF EXISTS raw.fps, raw.mols, raw.props;
    DROP SCHEMA IF EXISTS raw;
    """
cursor.execute(init_dud)
connect.commit()
# load smiles and props into temp database
# using zinc_id as primary key, do nothing when same zinc_id.
for i, sql_file in enumerate(args.sql):
    sql_file = Path(sql_file).absolute()
    restore_sql = 'gunzip < {} | psql -p {} {}'.format(sql_file, args.port,
                                                       args.dbname)
    subprocess.check_call(restore_sql, shell=True)
    insert_query = """
    INSERT INTO dud.props (SELECT * FROM raw.props)
        ON CONFLICT (zinc_id) DO NOTHING;
    INSERT INTO dud.fps (SELECT * FROM raw.fps)
        ON CONFLICT (zinc_id) DO NOTHING;
    INSERT INTO dud.mols (SELECT * FROM raw.mols)
        ON CONFLICT (zinc_id) DO NOTHING;
    DROP TABLE raw.fps, raw.mols, raw.props;
    DROP SCHEMA raw;
    """
    cursor.execute(insert_query)
    connect.commit()
    print("{}: Loading {:4d}/{} {}".format(dt.now() - start, i, len(args.sql),
                                           sql_file))

create_indexes = """
CREATE INDEX IF NOT EXISTS fps_mfp2_idx ON dud.fps USING gist(mfp2);
CREATE INDEX IF NOT EXISTS molidx ON dud.mols USING gist(m);
CREATE INDEX IF NOT EXISTS prop_idx ON dud.props (mw, logp, rotb, hba, hbd, q);
"""
cursor.execute(create_indexes)
connect.commit()
connect.close()

print("Total elapsed time: {}".format(dt.now() - start))

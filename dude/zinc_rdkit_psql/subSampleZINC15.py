"""Generating decoys from database based on active smiles.
"""
import argparse
from pathlib import Path
from datetime import datetime as dt

import psycopg2
import psycopg2.extras
psycopg2.extensions.set_wait_callback(psycopg2.extras.wait_select)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-s', '--schema', help="new schema name")
parser.add_argument('-f',
                    '--fraction',
                    default=16646056. / 1361729919,
                    type=float,
                    help="fraction of zinc15")
parser.add_argument("-d", "--dbname", required=True)
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
args = parser.parse_args()

start = dt.now()

if args.schema is None:
    args.schema = 'sub{:03d}'.format(int(args.fraction * 1000))

sub_sample = """
CREATE EXTENSION IF NOT EXISTS rdkit;
CREATE SCHEMA {schema};

CREATE TABLE IF NOT EXISTS {schema}.props (
    zinc_id integer PRIMARY KEY,
    smiles text,
    mw real,
    logp real,
    rotb smallint,
    hbd smallint,
    hba smallint,
    q smallint);

CREATE TABLE IF NOT EXISTS {schema}.fps (
    zinc_id integer PRIMARY KEY,
    mfp2 bfp);

CREATE TEMP TABLE sub_ids (
    zinc_id integer PRIMARY KEY);

INSERT INTO sub_ids
SELECT zinc_id FROM dud.props WHERE RANDOM() <= {fraction};

INSERT INTO {schema}.props
SELECT t.* FROM dud.props AS t
 WHERE t.zinc_id IN (SELECT zinc_id FROM sub_ids);

INSERT INTO {schema}.fps
SELECT t.* FROM dud.fps AS t
 WHERE t.zinc_id IN (SELECT zinc_id FROM sub_ids);

CREATE INDEX fps_mfp2_idx ON {schema}.fps USING gist(mfp2);
CREATE INDEX prop_idx ON {schema}.props (mw, logp, rotb, hba, hbd, q);
""".format(schema=args.schema, fraction=args.fraction)

connect = psycopg2.connect(host=args.host, dbname=args.dbname, port=args.port)
cursor = connect.cursor()
print('executing sql:\n' + sub_sample + '\nplease wait ...')
cursor.execute(sub_sample)
connect.commit()
connect.close()
print("Total elapsed time: {}".format(dt.now() - start))

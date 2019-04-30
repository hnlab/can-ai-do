"""Generating decoys from database based on active smiles.
"""
import argparse
from pathlib import Path
from datetime import datetime as dt

import psycopg2
import psycopg2.extras
psycopg2.extensions.set_wait_callback(psycopg2.extras.wait_select)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-s', '--schema', help="new schema name", required=True)
parser.add_argument('-S', '--size', default=1000000, type=int)
parser.add_argument("-d", "--dbname", required=True)
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
args = parser.parse_args()

start = dt.now()

split_query = """
CREATE EXTENSION IF NOT EXISTS rdkit;
CREATE SCHEMA IF NOT EXISTS {schema};

SELECT props.*, fps.mfp2 INTO {schema}.part{part:03d}
FROM dud.props JOIN dud.fps USING (zinc_id)
WHERE zinc_id > {max_id}
ORDER BY zinc_id
LIMIT {size};

CREATE OR REPLACE FUNCTION check_max_id() RETURNS int AS $$
BEGIN
    IF EXISTS (SELECT 1 FROM {schema}.part{part:03d}) THEN
        RETURN (SELECT MAX(zinc_id) FROM {schema}.part{part:03d});
    ELSE
        DROP TABLE {schema}.part{part:03d};
        RETURN NULL;
    END IF;
END $$ LANGUAGE plpgsql;

SELECT check_max_id();
"""

connect = psycopg2.connect(host=args.host, dbname=args.dbname, port=args.port)
cursor = connect.cursor()
print('Splitting into schema {} with size {}:'.format(args.schema, args.size))

part_count = 0
max_id = -1
while max_id is not None:
    print("{} working on part: {}.part{:03d}".format(dt.now() - start,
                                                     args.schema, part_count))
    cursor.execute(
        split_query.format(schema=args.schema,
                           part=part_count,
                           size=args.size,
                           max_id=max_id))
    max_id = cursor.fetchone()[0]
    part_count += 1
connect.commit()
connect.close()
print("Total elapsed time: {}".format(dt.now() - start))

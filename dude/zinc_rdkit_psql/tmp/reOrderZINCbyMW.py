"""UNION ALL zinc mols ORDER BY mw and split into parts again."""
import json
import gzip
import random
import argparse
import numpy as np
from tqdm import tqdm
from pathlib import Path
from datetime import datetime as dt

from multiprocessing import Pool

import psycopg2
import psycopg2.extras
psycopg2.extensions.set_wait_callback(psycopg2.extras.wait_select)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-d", "--dbname", required=True)
parser.add_argument("-s", "--schema", required=True)
parser.add_argument("--host", default='localhost')
parser.add_argument("-p", "--port", default='5432')
parser.add_argument("-n", "--name", default='zinc')
args = parser.parse_args()

start = dt.now()

connect = psycopg2.connect(host=args.host, dbname=args.dbname, port=args.port)
cursor = connect.cursor()

# cursor.execute("""
# SELECT table_name
#   FROM information_schema.tables
#  WHERE table_name ~ 'part'
#    AND table_schema = '{schema}'
# """.format(schema=args.schema))
# part_tables = [i[0] for i in cursor.fetchall()]

# part_tables = part_tables
# select_from = []
# for part in part_tables:
#     select_from.append('SELECT * FROM {schema}.{part}\n'.format(
#         schema=args.schema, part=part))
# union_all = 'UNION ALL\n'.join(select_from)

# union_by_mw = """
# -- DROP TABLE IF EXISTS zinc_by_mw;
# SELECT row_number() OVER (ORDER BY mw) AS mw_row_num, * INTO zinc_by_mw FROM (
# {union_all}
# ) AS tmp
# ORDER BY mw;
# CREATE INDEX zinc_by_mw_idx ON zinc_by_mw (mw_row_num);
# """
# cursor.execute(union_by_mw.format(union_all=union_all))
# connect.commit()

# for part in part_tables:
#     cursor.execute('DROP TABLE {schema}.{part}'.format(schema=args.schema,
#                                                        part=part))
#     connect.commit()

select_by_mw = """
SET enable_seqscan TO OFF;
DROP TABLE IF EXISTS {part};
SELECT zinc_id, smiles, mw, logp, rotb, hbd, hba, q, mfp2 INTO {part} FROM zinc_by_mw
WHERE mw_row_num > {offset} AND mw_row_num <= {offset} + {limit}
ORDER BY mw_row_num;

CREATE INDEX {part}_mw_idx ON {part} (mw);
CREATE INDEX {part}_id_idx ON {part} (zinc_id);
CREATE INDEX {part}_fp_idx ON {part} USING gist(mfp2);
"""
limit = 1000000
cursor.execute('SELECT max(mw_row_num) FROM zinc_by_mw')
N = cursor.fetchone()[0]
if N % limit == 0:
    npart = N // limit
else:
    npart = N // limit + 1

select_query = []
for i in range(npart):
    part = '{name}_split1m_{i:05d}'.format(name=args.name, i=i)
    offset = limit * i
    query = select_by_mw.format(part=part, limit=limit, offset=offset)
    select_query.append(query)


def run_query_help(query):
    _connect = psycopg2.connect(host=args.host,
                                dbname=args.dbname,
                                port=args.port)
    _cursor = _connect.cursor()
    _cursor.execute(query)
    _connect.commit()
    _connect.close()
    return None


p = Pool()
for _ in tqdm(p.imap(run_query_help, select_query), total=npart):
    pass

connect.commit()
connect.close()
print("Total elapsed time: {}".format(dt.now() - start))

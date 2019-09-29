"""Load mols from smiles into database. 50000/min, about 5.5 hours for zinc12.
"""
import gzip
import shutil
import tempfile
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


def getProp(mol_line):
    smiles, mol_id = mol_line.split()[:2]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    mol_id = mol_id.replace('ZINC', '')
    mol_id = int(mol_id)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    q = Chem.GetFormalCharge(mol)  # can not get charge in rdkit Cartridge
    return tuple([mol_id, smiles, mw, logp, rotb, hbd, hba, q])


def props_generator_from_files(smiles_files):
    for smi in smiles_files:
        print(f"{dt.now()}: loading {smi}")
        smi = Path(smi)
        if smi.suffix == '.gz':
            f = gzip.open(smi, 'rt')
        else:
            f = open(smi, 'r')
        for props in map(getProp, f):
            if props is not None:
                yield props


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--smiles", nargs='+', required=True)
    parser.add_argument("-D",
                        "--dbpath",
                        required=True,
                        help="path of database")
    parser.add_argument("-d",
                        "--dbname",
                        default="zinc",
                        help="name of database, default=zinc")
    parser.add_argument("-n",
                        "--schema",
                        default="zinc",
                        help="name of schema, default=zinc")
    parser.add_argument("--host", default='localhost')
    parser.add_argument("-p", "--port", default='5432')
    args = parser.parse_args()

    start = dt.now()

    # init database
    dbpath = Path(args.dbpath)
    if dbpath.exists():
        raise Exception("Database already exists, may try:\n\n"
                        f"pg_ctl -o '-F -p {args.port}' -D {dbpath} stop\n"
                        f"rm -r {args.dbpath}\n")
    dbpath.mkdir(parents=True)
    subprocess.check_call(['initdb', '-D', dbpath])
    # -F : Disables fsync calls for improved performance, at the risk of data corruption in the event of a system crash.
    # -o : pass options to postgres, https://www.postgresql.org/docs/8.4/reference-server.html
    subprocess.check_call([
        'pg_ctl', '-o', f'-F -p {args.port}', '-D',
        str(dbpath), '-l', 'zinc.log', 'start'
    ])
    subprocess.check_call(['createdb', '-p', args.port, args.dbname])
    connect = psycopg2.connect(host='localhost',
                               dbname=args.dbname,
                               port=args.port)
    cursor = connect.cursor()

    # init schema and table props
    cursor.execute(f"CREATE EXTENSION IF NOT EXISTS rdkit;"
                   f"CREATE SCHEMA IF NOT EXISTS {args.schema};"
                   f"CREATE TABLE IF NOT EXISTS {args.schema}.props ("
                   f"       zinc_id integer PRIMARY KEY,"
                   f"        smiles text,"
                   f"            mw real,"
                   f"          logp real,"
                   f"          rotb smallint,"
                   f"           hbd smallint,"
                   f"           hba smallint,"
                   f"             q smallint)")

    print(f"{dt.now()}: "
          f"load smiles and props into {args.dbname}.{args.schema}.props")
    # using zinc_id as primary key, do nothing when same zinc_id.
    props_generator = props_generator_from_files(args.smiles)
    insert_query = f'INSERT INTO {args.schema}.props VALUES %s ON CONFLICT (zinc_id) DO NOTHING'
    psycopg2.extras.execute_values(cursor,
                                   insert_query,
                                   props_generator,
                                   template=None,
                                   page_size=512)

    print(f"{dt.now()}: generate fps for zinc full set")
    # follow https://www.rdkit.org/docs/Cartridge.html#loading-chembl
    cursor.execute(
        f"SELECT *, morganbv_fp(mol_from_smiles(smiles::cstring)) as mfp2"
        f"  INTO {args.schema}.full"
        f"  FROM {args.schema}.props;"
        f"CREATE INDEX fp_idx_full on {args.schema}.full USING gist(mfp2);"
        f"CREATE INDEX prop_idx_full on {args.schema}.full (mw, logp, rotb, hba, hbd, q);"
    )
    connect.commit()
    print(f"{dt.now()}: zinc full set at {args.dbname}.{args.schema}.full")
    # create zinc drug like subset
    cursor.execute(
        f"SELECT *"
        f"  INTO {args.schema}.drug_like"
        f"  FROM {args.schema}.full"
        f" WHERE mw <=500 AND mw >= 150"
        f"   AND logp <= 5 AND rotb <= 7"
        f"   AND hbd <= 5 AND hba <= 10;"
        f"CREATE INDEX fp_idx_drug on {args.schema}.drug_like USING gist(mfp2);"
        f"CREATE INDEX prop_idx_drug on {args.schema}.drug_like (mw, logp, rotb, hba, hbd, q);"
    )
    connect.commit()
    print(f"{dt.now()}: zinc drug_like set at "
          f"{args.dbname}.{args.schema}.drug_like")
    connect.close()
    print("\npostgres database had started by:\n"
          f"pg_ctl -o '-F -p {args.port}' -D {dbpath} -l zinc.log start\n")
    print(f"Total elapsed time: {dt.now() - start}")
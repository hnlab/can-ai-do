"""Load mols from smiles into database.
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

import socket
def get_free_tcp_port():
    tcp = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    tcp.bind(('', 0))
    addr, port = tcp.getsockname()
    tcp.close()
    return port

def getProp(mol_line):
    smiles, mol_id = mol_line.split()[:2]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    mol_id = mol_id.replace('ZINC', '')
    mol_id = int(mol_id)
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--smiles", nargs='+', required=True)
    parser.add_argument("-o", "--output", default='raw')
    args = parser.parse_args()
    
    start = dt.now()

    # init a temp postdatabase
    import time
    dbpath = tempfile.mkdtemp()
    dbname = 'temp'
    subprocess.call(['initdb', '-D', dbpath])

    # try 3 time until pg_ctl start at a free port.
    code = 404
    try_count = 3
    while code != 0 and try_count > 0:
        try_count -= 1
        port = str(get_free_tcp_port())
        port_option = "-F -p {}".format(port)
        try:
            code = subprocess.call(['pg_ctl', '-o', port_option, '-D', dbpath,'-l','logfile', 'start'])
        except:
            pass
    # sleep 3 seconds to wait pg_ctl start succeed.
    time.sleep(3)
    subprocess.call(['createdb', '-p', port, dbname])
    time.sleep(3)
    connect = psycopg2.connect(host='localhost', dbname=dbname, port=port)
    cursor = connect.cursor()

    # init a temp schema raw
    cursor.execute("""
        CREATE EXTENSION IF NOT EXISTS rdkit;
        CREATE SCHEMA raw;
        CREATE TABLE IF NOT EXISTS raw.props (
            zinc_id integer PRIMARY KEY,
            smiles text,
            mw real,
            logp real,
            rotb smallint,
            hbd smallint,
            hba smallint,
            q smallint)
        """)

    # load smiles and props into temp database
    # using zinc_id as primary key, do nothing when same zinc_id.
    props_generator = props_generator_from_files(args.smiles)
    insert_query = 'INSERT INTO raw.props values %s ON CONFLICT (zinc_id) DO NOTHING'
    psycopg2.extras.execute_values(
        cursor, insert_query, props_generator, template=None, page_size=100)

    # generate fps and mols
    # follow https://www.rdkit.org/docs/Cartridge.html#loading-chembl
    sql = """
    SELECT * INTO raw.mols
    FROM (SELECT zinc_id, mol_from_smiles(smiles::cstring) m FROM raw.props) tmp 
    WHERE m IS NOT NULL;
    
    SELECT zinc_id, morganbv_fp(m) as mfp2 INTO raw.fps FROM raw.mols;
    """
    cursor.execute(sql)
    connect.commit()
    connect.close()
    
    # dump data and delete database
    output = Path(args.output).with_suffix('.sql.gz')
    subprocess.call(['pg_dump', '-p', port, '--no-owner', '-Z', '9', dbname, '-f', output])
    subprocess.call(['pg_ctl', '-o', port_option, '-D', dbpath, '-l', 'logfile', 'stop'])
    shutil.rmtree(dbpath)
    print("Total elapsed time: {}".format(dt.now()-start))

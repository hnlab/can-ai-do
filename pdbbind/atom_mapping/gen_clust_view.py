import subprocess
from pathlib import Path
from multiprocessing.dummy import Pool
index = './index/INDEX_core_cluster.2013.clust'
root = './pdb'


def clust2view(clust):
    ref = clust[0]
    for id_ in clust:
        subprocess.check_call([
            'python', 'align_lig.py', '-l', f'{root}/atom.{id_}_ligand.pdb',
            '-p', f'{root}/atom.{id_}_pocket.pdb', '-r',
            f'{root}/atom.{ref}_pocket.pdb'
        ])
    subprocess.check_call(['python', 'occ_view.py', '-i'] +
                          [f'{root}/atom.{i}_ligand.align.pdb' for i in clust])

    # subprocess.check_call(['charge_view.py', '-i'] +
    #                       [f'{root}/atom.{i}_ligand.align.pdb' for i in clust])

clusters = [i.split() for i in open(index)]
p = Pool()
p.map(clust2view, clusters)
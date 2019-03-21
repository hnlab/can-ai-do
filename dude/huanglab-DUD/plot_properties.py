"""Deduping decoys with same zinc id in same target.
"""
import yaml
import argparse
import numpy as np
from pathlib import Path

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors as D
from rdkit.Chem import rdMolDescriptors as CD

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "-t", "--target_dir", nargs='+', required=True, help="target dir names")
parser.add_argument(
    "-o", "--output", help="output dir, default will save figure in each target dir.")
parser.add_argument(
    "-d",
    "--decoys_file",
    default="decoys_final.smi",
    help="decoys file name, default: decoys_final.smi" 
)
parser.add_argument(
    "-a",
    "--actives_file",
    default="actives_final.smi",
    help="actives file name, default: actives_final.smi" 
)
args = parser.parse_args()

def get_prop_array(mol):
    mw = CD.CalcExactMolWt(mol)
    logp = Chem.Crippen.MolLogP(mol)
    rotb = D.NumRotatableBonds(mol)
    hbd = CD.CalcNumHBD(mol)
    hba = CD.CalcNumHBA(mol)
    q = Chem.GetFormalCharge(mol)
    return np.array([mw, logp, rotb, hbd, hba, q])

if args.output:
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
for tdir in args.target_dir:
    tdir = Path(tdir)
    print(tdir)
    actives_file = list(tdir.glob("actives_final.*"))[0]
    decoys_file = list(tdir.glob("decoys_final.*"))[0]
    actives_props = []
    decoys_props = []
    for m in Chem.SmilesMolSupplier(str(actives_file), titleLine=False):
        if m is not None:
            actives_props.append(get_prop_array(m))
    for m in Chem.SmilesMolSupplier(str(decoys_file), titleLine=False):
        if m is not None:
            decoys_props.append(get_prop_array(m))
    actives_props = np.array(actives_props)
    decoys_props = np.array(decoys_props)
    props_name = ["mw", "logp", "rotb", "hbd", "hba", "q"]
    fig, axes = plt.subplots(
            nrows=2, ncols=3, figsize=(24, 12))
    axes = axes.flatten()
    for i, p in enumerate(props_name):
        ax = axes[i]
        if p in ["mw", "logp"]:
            sns.kdeplot(actives_props[:,i], label=p+'_active', color="blue", ax=ax)
            sns.kdeplot(decoys_props[:,i], label=p+'_decoy', color="red",  ax=ax)
        else:
            a_props = list(map(int, actives_props[:,i]))
            d_props = list(map(int, decoys_props[:,i]))
            prop_max = max(max(a_props),max(d_props))
            prop_min = min(min(a_props),min(d_props))
            prop_xticks = np.arange(prop_min, prop_max+1)
            prop_bins = prop_xticks - 0.5
            sns.distplot(a_props, bins=prop_bins, label=p+'_active', color="blue", kde=False, norm_hist=True, ax=ax)
            sns.distplot(d_props, bins=prop_bins, label=p+'_decoy', color="red", kde=False, norm_hist=True, ax=ax)
            ax.legend()
            ax.set_xticks(prop_xticks)
            # a_count = list(map(int, actives_props[:,i]))
            # d_count = list(map(int, decoys_props[:,i]))
            # ax = axes[i]
            # ax.hist(a_count)
            # ax.hist(d_count)
        ax.autoscale(enable=True, axis='y')
    fig_path = tdir/"props.png"
    fig_path = fig_path.resolve()
    fig.savefig(fig_path)
    if args.output:
        fig_path = output/ "_".join(fig_path.parts[-3:])
        fig.savefig(fig_path)
    print("figure saved at {}".format(fig_path))

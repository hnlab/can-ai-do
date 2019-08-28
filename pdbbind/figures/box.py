#%%
import json
import numpy as np
import pandas as pd
from pathlib import Path

import seaborn as sns
from matplotlib import pyplot as plt

#%%
print(Path.cwd())
root = Path('pdbbind/figures')
result = root / 'result'
from1 = result / 'result.runEarlyStopFrom1stEpoch.csv'

# use data from same_ligand replace 'protein' only in from1.
# `binding` in same_ligand in better.
# more details in https://github.com/hnlab/deepchem/pull/1#issuecomment-510721589
same_lig = result / 'result.same_ligand.csv'

df_from1 = pd.read_csv(from1)
df_same_lig = pd.read_csv(same_lig)
print(df_same_lig[:3])
#%%

protein_alone = df_same_lig.loc[df_same_lig['component'] == 'binding']
protein_alone.loc[protein_alone['component'] ==
                  'binding', 'component'] = 'protein'
print(len(df_same_lig))
print(len(protein_alone))

others = df_from1.loc[df_from1['component'] != 'protein']
print(len(df_from1))
print(len(others))

df = pd.concat([others, protein_alone])
print(len(df))
print(df[:3])

#%%
components = ['binding', 'ligand', 'protein']
cpn_names = ['binding complexes', 'ligands alone', 'proteins alone']
subsets = ['core', 'refined', 'general_PL']
splits = ['random', 'scaffold', 'seq']
split_names = ['random splitting', 'scaffold splitting', 'sequence splitting']
metrics = ['pearson_r2_score', 'mean_absolute_error']
metric_names = ['pearson $R^2$', 'MAE']

#%%
for metric, metric_name in zip(metrics, metric_names):
    for split, split_name in zip(splits, split_names):
        d = df.loc[df['set'] == 'test']
        d = d.loc[d['split'] == split]
        d = d.loc[d['metric'] == metric]
        # print(d)
        # showfliers=False hide outliers
        fig, ax = plt.subplots()
        sns.boxplot(
            x='subset',
            y='value',
            hue='component',
            hue_order=components,
            order=subsets,
            data=d,
            showfliers=False,
            palette="Set2",
            ax=ax,
        )
        # dodge split points
        ax = sns.swarmplot(
            x='subset',
            y='value',
            hue='component',
            hue_order=components,
            order=subsets,
            data=d,
            dodge=True,
            linewidth=1.5,
            edgecolor='gray',
            palette="Set2",
            # palette=['black', 'black'],
            ax=ax,
        )
        # Get the handles and labels. For this example it'll be 2 tuples
        # of length 4 each.
        handles, labels = ax.get_legend_handles_labels()
        # When creating the legend, only use the first two elements
        # to effectively remove the last two.
        # l = ax.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        loc = 'best'
        if metric == 'pearson_r2_score':
            loc = 'lower center'
        l = ax.legend(handles[0:3], split_names, frameon=False, loc=loc)
        ax.set_xlabel('')
        ax.set_ylabel(metric_name)
        ax.set_title(
            f'ACNN performance on test sets of PDBbind subsets\nwith {split_name}'
        )
        # ax.set_ylim([0.4,1])
        fig.savefig(root / '.'.join((split, metric, 'png')), dpi=300)

#%%
for metric, metric_name in zip(metrics, metric_names):
    for component, component_name in zip(components, cpn_names):
        d = df.loc[df['set'] == 'test']
        d = d.loc[d['component'] == component]
        d = d.loc[d['metric'] == metric]
        # print(d)
        # showfliers=False hide outliers
        fig, ax = plt.subplots()
        sns.boxplot(
            x='subset',
            y='value',
            hue='split',
            hue_order=splits,
            order=subsets,
            data=d,
            showfliers=False,
            palette="Set2",
            ax=ax,
        )
        # dodge split points
        ax = sns.swarmplot(
            x='subset',
            y='value',
            hue='split',
            hue_order=splits,
            order=['core', 'refined', 'general_PL'],
            data=d,
            dodge=True,
            linewidth=1.5,
            edgecolor='gray',
            palette="Set2",
            # palette=['black', 'black'],
            ax=ax,
        )
        # Get the handles and labels. For this example it'll be 2 tuples
        # of length 4 each.
        handles, labels = ax.get_legend_handles_labels()
        # When creating the legend, only use the first two elements
        # to effectively remove the last two.
        # l = ax.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        loc = 'best'
        if metric == 'pearson_r2_score':
            loc = 'lower center'
        l = ax.legend(handles[0:3], split_names, frameon=False, loc=loc)
        ax.set_xlabel('')
        ax.set_ylabel(metric_name)
        ax.set_title(
            f'ACNN performance on test sets of PDBbind subsets\ncontaining {component_name} '
        )
        # ax.set_ylim([0.4,1])
        fig.savefig(root / '.'.join((component, metric, 'png')), dpi=300)

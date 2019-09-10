#%%
import json
import numpy as np
import pandas as pd
from pathlib import Path

import seaborn as sns
from matplotlib import pyplot as plt

#%%
print(Path.cwd())
root = Path('pdbbind/pred_core.figures')
result = root / 'result'
csv = result / 'pred_core2015.csv'
df = pd.read_csv(csv)
df.tail()

#%%
components = ['binding', 'ligand', 'protein']
cpn_names = ['binding complexes', 'ligands alone', 'proteins alone']
subsets = ['refined', 'general_PL']
splits = ['random', 'fingerprint', 'sequence']
split_names = ['random splitting', 'fingerprint splitting', 'sequence splitting']
# metrics = ['pearson_r2_score', 'mean_absolute_error']
# metric_names = ['pearson $R^2$', 'MAE']

# MAE no meaning in here.
metrics = ['pearson_r2_score']
metric_names = ['pearson $R^2$']
version = 2015
# version = 2018
#%%
for metric, metric_name in zip(metrics, metric_names):
    for split, split_name in zip(splits, split_names):
        d = df.loc[df['set'] == 'test']
        d = d.loc[d['split'] == split]
        d = d.loc[d['metric'] == metric]
        d = d.loc[d['version'] == version]
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
            loc = 'upper left'
        l = ax.legend(handles[0:3], cpn_names, frameon=False, loc=loc)
        ax.set_xlabel('')
        ax.set_ylabel(metric_name)
        ax.set_title(
            f'ACNN performance on test sets of PDBbind {version} subsets\nwith {split_name}'
        )
        ax.set_ylim([-0.05,1.05])
        fig.savefig(root / f"{version}.{split}.{metric}.png", dpi=300)

#%%
for metric, metric_name in zip(metrics, metric_names):
    for component, component_name in zip(components, cpn_names):
        d = df.loc[df['set'] == 'test']
        d = d.loc[d['component'] == component]
        d = d.loc[d['metric'] == metric]
        d = d.loc[d['version'] == version]
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
            loc = 'upper left'
        l = ax.legend(handles[0:3], split_names, frameon=False, loc=loc)
        ax.set_xlabel('')
        ax.set_ylabel(metric_name)
        ax.set_title(
            f'ACNN performance on test sets of PDBbind {version} subsets\ncontaining {component_name} '
        )
        ax.set_ylim([-0.05,1.05])
        fig.savefig(root / f"{version}.{component}.{metric}.png", dpi=300)

#%%
d = df.loc[df['set'] == 'test']
d = d.loc[d['metric'] == 'pearson_r2_score']
pt = pd.pivot_table(d, index=['version', 'split', 'component'], columns=['subset'], values=['value'])
pt.to_csv(root/'pivot_table.csv')
pt
#%%
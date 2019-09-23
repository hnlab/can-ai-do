#%%
import json
import numpy as np
import pandas as pd
from pathlib import Path

import seaborn as sns
from matplotlib import pyplot as plt

#%%
print(Path.cwd())
root = Path('pdbbind/self.figures')
result = root / 'result'
# csv = result / 'pred_self.csv' # seq simi 0.5, fp simi 0.5

csv = result / 'self.scaffold.csv'

df = pd.read_csv(csv)

df.tail()

#%%
components = ['binding', 'ligand', 'protein']
cpn_names = ['binding complexes', 'ligands alone', 'proteins alone']
subsets = ['core', 'refined', 'general_PL']
# splits = ['random', 'fingerprint', 'sequence']
# split_names = ['random splitting', 'fingerprint splitting', 'sequence splitting']

splits = ['random', 'scaffold', 'sequence']
split_names = ['random splitting', 'scaffold splitting', 'sequence splitting']
# metrics = ['pearson_r2_score', 'mean_absolute_error']
# metric_names = ['pearson $R^2$', 'MAE']

# MAE no meaning in here.
metrics = ['pearson_r2_score']
metric_names = ['pearson $R^2$']
# version = 2015
version = 2018
sns.palplot(sns.color_palette('Set2'))
sns.palplot(sns.color_palette('Paired', 9))

colors = sns.color_palette("Set2", 9)
sns.palplot(sns.light_palette(colors[0], 3))
sns.palplot(sns.dark_palette(colors[0], 3))

colors_2d = []
for c in colors[:3]:
    colors_2d.extend(sns.light_palette(c, 5)[-3:])

sns.palplot(colors_2d)
colors_2d = np.reshape(colors_2d, (3, 3, 4))

colors_2d = sns.color_palette('Paired', 9)
colors_2d = sns.color_palette()[:9]
colors_2d = np.reshape(colors_2d, (3, 3, 3))

#%%
for metric, metric_name in zip(metrics, metric_names):
    for split, split_name, cs in zip(splits, split_names, colors_2d):
        d = df.loc[df['set'] == 'test']
        d = d.loc[d['split'] == split]
        d = d.loc[d['metric'] == metric]
        d = d.loc[d['version'] == version]
        # print(d)
        # showfliers=False hide outliers
        fig, ax = plt.subplots()
        cs_dict = {k: v for k, v in zip(components, cs)}
        sns.boxplot(
            x='subset',
            y='value',
            hue='component',
            hue_order=components,
            order=subsets,
            data=d,
            showfliers=False,
            # color=colors[0],
            # palette="Set2",
            palette=cs_dict,
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
            # palette="Set2",
            palette=cs_dict,
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
            loc = 'lower right'
        l = ax.legend(handles[0:3], cpn_names, frameon=False, loc=loc)
        ax.set_xlabel('')
        ax.set_ylabel(metric_name)
        ax.set_title(
            f'ACNN performance on test sets of PDBbind {version} subsets\nwith {split_name}'
        )
        ax.set_ylim([-0.05, 1.05])
        fig.savefig(root / f"{version}.{split}.{metric}.png", dpi=300)

#%%
colors_2dT = np.transpose(colors_2d, (1, 0, 2))
for metric, metric_name in zip(metrics, metric_names):
    for component, component_name, cs in zip(components, cpn_names,
                                             colors_2dT):
        d = df.loc[df['set'] == 'test']
        d = d.loc[d['component'] == component]
        d = d.loc[d['metric'] == metric]
        d = d.loc[d['version'] == version]
        # print(d)
        cs_dict = {k: v for k, v in zip(splits, cs)}
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
            palette=cs_dict,
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
            palette=cs_dict,
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
            loc = 'lower right'
        l = ax.legend(handles[0:3], split_names, frameon=False, loc=loc)
        ax.set_xlabel('')
        ax.set_ylabel(metric_name)
        ax.set_title(
            f'ACNN performance on test sets of PDBbind {version} subsets\ncontaining {component_name} '
        )
        ax.set_ylim([-0.05, 1.05])
        fig.savefig(root / f"{version}.{component}.{metric}.png", dpi=300)

#%%
d = df.loc[df['set'] == 'test']
d = d.loc[d['metric'] == 'pearson_r2_score']
pt = pd.pivot_table(d,
                    index=['version', 'split', 'component'],
                    columns=['subset'],
                    values=['value'],
                    aggfunc=[np.mean, np.std])
pt.to_csv(root / 'pivot_table.csv')
pt
#%%
#%%
import json
import numpy as np
import pandas as pd
from pathlib import Path

import seaborn as sns
from matplotlib import pyplot as plt

#%%
print(Path.cwd())
root = Path('dude/figures')
files = {
    'DUD-E': 'result/DUDE.json',
    'drug-like': 'result/directT8Z12N62DL.json',
    'no-limit': 'result/directT8Z12N62.json',
    # 'no-limit-ZINC15': 'result/directT8Z15N62.json'
}
results = {}
for k, v in files.items():
    with open(root / v) as f:
        results[k] = json.load(f)
print(json.dumps(results['DUD-E'], indent=2))

#%%
data = []
for dataset, result in results.items():
    for fold in result:
        if 'fold' not in fold:
            continue
        for feature in ('fp', 'prop'):
            for metric in ('EF1', 'ROC'):
                value = fold[feature][metric]
                data.append((dataset, fold, feature, metric, value))
df = pd.DataFrame(
    columns=('dataset', 'fold', 'feature', 'metric', 'value'), data=data)

#%%
d = df.loc[df['metric'] == 'ROC']
# showfliers=False hide outliers
ax = sns.boxplot(
    x='dataset',
    y='value',
    hue='feature',
    hue_order=['prop','fp'],
    data=d,
    showfliers=False,
    palette="Set2",
)
# dodge split points
ax = sns.swarmplot(
    x='dataset',
    y='value',
    hue='feature',
    hue_order=['prop','fp'],
    data=d,
    dodge=True,
    linewidth=1.5,
    edgecolor='gray',
    palette="Set2",
    # palette=['black', 'black'],
)
# Get the handles and labels. For this example it'll be 2 tuples
# of length 4 each.
handles, labels = ax.get_legend_handles_labels()
# When creating the legend, only use the first two elements
# to effectively remove the last two.
# l = ax.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
l = ax.legend(handles[0:2], ['PROP+RF','FP+RF'], frameon=False)
ax.set_xlabel('Diverse subsets (8 targets)')
ax.set_ylabel('AUC')
# ax.set_ylim([0.4,1])
plt.savefig(root/'AUC.png', dpi=300)

#%%
d = df.loc[df['metric'] == 'EF1']
# showfliers=False hide outliers
ax = sns.boxplot(
    x='dataset',
    y='value',
    hue='feature',
    hue_order=['prop','fp'],
    data=d,
    showfliers=False,
    palette="Set2",
)
# dodge split points
ax = sns.swarmplot(
    x='dataset',
    y='value',
    hue='feature',
    hue_order=['prop','fp'],
    data=d,
    dodge=True,
    linewidth=1.5,
    edgecolor='gray',
    palette="Set2",
    # palette=['black', 'black'],
)
# Get the handles and labels. For this example it'll be 2 tuples
# of length 4 each.
handles, labels = ax.get_legend_handles_labels()
# When creating the legend, only use the first two elements
# to effectively remove the last two.
# l = ax.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# frameon=False, no box around legend
l = ax.legend(handles[0:2], ['PROP+RF','FP+RF'], frameon=False)
ax.set_xlabel('Diverse subsets (8 targets)')
ax.set_ylabel('EF1%')
plt.savefig(root/'EF1.png', dpi=300)
#%%

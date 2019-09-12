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
df = pd.DataFrame(columns=('dataset', 'fold', 'feature', 'metric', 'value'),
                  data=data)

#%%
d = df.loc[df['feature'] == 'prop']
roc = d.loc[df['metric'] == 'ROC']
ef1 = d.loc[df['metric'] == 'EF1']
fig, ax = plt.subplots()
ax2 = ax.twinx()
width = 0.4
colors = sns.color_palette("Set2")
sns.boxplot(
    x='dataset',
    y='value',
    data=ef1,
    ax=ax,
    showfliers=False,
    positions=[-0.2, 0.8, 1.8],
    width=0.35,
    color=colors[0],
)

sns.swarmplot(
    x='dataset',
    y='value',
    data=ef1,
    ax=ax,
    dodge=True,
    linewidth=1.5,
    edgecolor='gray',
    color=colors[0],
    positions=[-0.2, 0.8, 1.8],
)

sns.boxplot(
    x='dataset',
    y='value',
    data=roc,
    ax=ax2,
    showfliers=False,
    positions=[0.2, 1.2, 2.2],
    width=0.35,
    # palette="Set2",
    color=colors[1],
)
sns.swarmplot(
    x='dataset',
    y='value',
    data=roc,
    ax=ax2,
    dodge=True,
    linewidth=1.5,
    edgecolor='gray',
    color=colors[1],
    positions=[0.2, 1.2, 2.2],
)

ax.set_xlabel('')
ax.set_ylabel('EF1')
ax.set_ylim([-6, 62])
ax2.set_ylabel('ROC')
ax2.set_ylim([0.46, 1.02])

ax.set_title("PROP+RF performance on three Diverse sets (8 targets)")

# https://github.com/0ut0fcontrol/seaborn/blob/master/seaborn/categorical.py
def add_legend_data(ax, color, label):
    """Add a dummy patch object so we can get legend data."""
    rect = plt.Rectangle([0, 0],
                         0,
                         0,
                         edgecolor='gray',
                         facecolor=color,
                         label=label)
    ax.add_patch(rect)
add_legend_data(ax, colors[0], 'EF1')
add_legend_data(ax, colors[1], 'ROC')
ax.legend(frameon=False)
fig.savefig(root/'PROP.png', dpi=300)
#%%
#%%
d = df.loc[df['feature'] == 'fp']
roc = d.loc[df['metric'] == 'ROC']
ef1 = d.loc[df['metric'] == 'EF1']
fig, ax = plt.subplots()
ax2 = ax.twinx()
width = 0.4
colors = sns.color_palette("Set2")
sns.boxplot(
    x='dataset',
    y='value',
    data=ef1,
    ax=ax,
    showfliers=False,
    positions=[-0.2, 0.8, 1.8],
    width=0.35,
    color=colors[0],
)

sns.swarmplot(
    x='dataset',
    y='value',
    data=ef1,
    ax=ax,
    dodge=True,
    linewidth=1.5,
    edgecolor='gray',
    color=colors[0],
    positions=[-0.2, 0.8, 1.8],
)

sns.boxplot(
    x='dataset',
    y='value',
    data=roc,
    ax=ax2,
    showfliers=False,
    positions=[0.2, 1.2, 2.2],
    width=0.35,
    # palette="Set2",
    color=colors[1],
)
sns.swarmplot(
    x='dataset',
    y='value',
    data=roc,
    ax=ax2,
    dodge=True,
    linewidth=1.5,
    edgecolor='gray',
    color=colors[1],
    positions=[0.2, 1.2, 2.2],
)

ax.set_xlabel('')
ax.set_ylabel('EF1')
ax.set_ylim([-6, 62])
ax2.set_ylabel('ROC')
ax2.set_ylim([0.46, 1.02])

ax.set_title("FP+RF performance on three Diverse sets (8 targets)")

# https://github.com/0ut0fcontrol/seaborn/blob/master/seaborn/categorical.py
def add_legend_data(ax, color, label):
    """Add a dummy patch object so we can get legend data."""
    rect = plt.Rectangle([0, 0],
                         0,
                         0,
                         edgecolor='gray',
                         facecolor=color,
                         label=label)
    ax.add_patch(rect)
add_legend_data(ax, colors[0], 'EF1')
add_legend_data(ax, colors[1], 'ROC')
ax.legend(frameon=False, loc='upper left')
fig.savefig(root/'FP.png', dpi=300)

#%%

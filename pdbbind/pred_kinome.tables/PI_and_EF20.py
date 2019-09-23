#%%
import numpy as np
import pandas as pd
from pathlib import Path

#%%
root = Path(
    "/pubhome/jcyang/git/deepchem/examples/pdbbind/withHs.tmp/kinome.jobs/")


def PI(E, P):
  w_sum = 0.
  wc_sum = 0.
  for i, (ei, pi) in enumerate(zip(E, P)):
    for j, (ej, pj) in enumerate(zip(E, P)):
      if j > i:
        de = ej - ei
        dp = pj - pi
        w = abs(de)
        # print(w, dp)
        if abs(dp) < 1e-5:
          c = 0
        elif de / dp > 0:
          c = 1
        else:
          c = -1
        w_sum += w
        wc_sum += w * c
  return wc_sum / w_sum


#%%
def csv2EF20(csv):
  df = pd.read_csv(
      csv,
      sep='-|,',
      engine='python',
      names=['kinase', 'chain', 'ligand', 'pK', 'pK_pred'],
      header=0)
  ligands = df.ligand.unique()
  EF20s = []
  PIs = []
  for lig in ligands:
    df_lig = df[df.ligand == lig]
    N = len(df_lig)
    n = int(N * 0.2)
    # print(df_lig.tail())
    df_lig = df_lig.sort_values(by='pK_pred', ascending=False)
    # print(df_lig.tail())
    head20 = df_lig[:n]
    t = sum(head20.pK > 7)
    T = sum(df_lig.pK > 7)
    EF20 = (t / n) / (T / N)
    pi = PI(df_lig.pK, df_lig.pK_pred)
    # print(f"{lig} enrich {t}/{n} from {T}/{N} with EF20 {EF20}, PI {pi}")
    yield lig, EF20, pi


#%%
all_EF20s = []
mean_scores = None
count = 0
data = []
for csv in root.glob("job*15*bind*/kinome.csv"):
  print(csv)
  # job.2015.refined.binding.random.555/kinome.csv
  seed = Path(csv).parent.name.split('.')[-1]
  for ligand, EF20, pi in csv2EF20(csv):
    data.append((ligand, EF20, pi, seed))

df = pd.DataFrame(data, columns=['ligand', 'EF20', 'PI', 'seed'])

pt = pd.pivot_table(df,
                    index=['ligand'],
                    values=['EF20','PI'],
                    aggfunc=[np.mean, np.std],
                    margins=True,
                    margins_name='mean')
pt.to_csv(root / 'pivot_table.csv')
pt

#%%

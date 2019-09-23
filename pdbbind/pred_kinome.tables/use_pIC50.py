#%%
import pandas as pd
from pathlib import Path
root = Path("/pubhome/jcyang/git/deepchem/examples/pdbbind/withHs.tmp/kinome.jobs/")
result = pd.read_table(root/"kinome-wan-with_rank.tsv")

result.tail()

#%%
import re
# '-51.982(Rank-38/pIC50-5)'
p = re.compile(r'([-.+\d]+).*pIC50-([-+.\d]+)')
pIC50s = pd.DataFrame()
binding_energies = pd.DataFrame()
for column_name in result:
    if '(' not in column_name:
        pIC50s[column_name] = result[column_name]
        binding_energies[column_name] =result[column_name]
    else:
        pIC50_column = []
        binding_column = []
        for value in result[column_name]:
            try:
                matches = re.findall(p,value)
                binding, pIC50 = matches[0]
            except TypeError as e:
                binding, pIC50 = value, value
            pIC50_column.append(float(pIC50))
            binding_column.append(float(binding))
        pIC50s[column_name] = pIC50_column
        binding_energies[column_name] = binding_column

pIC50s.to_csv(root/'pIC50.tsv', sep='\t', index=False, na_rep='NA')
binding_energies.to_csv(root/'energies.tsv', sep='\t', index=False, na_rep='NA')
#%%


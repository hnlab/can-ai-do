# 1.Download ZINC12 all smiles;
- 1. goto http://zinc12.docking.org/subsets/all-purchasable
- 2. click download button
- 3. click [All](http://zinc12.docking.org/db/bysubset/6/6_p0.smi.gz) under "select Reference (pH 7)"
- 4. download file named "6_p0.smi.gz", 142 MB

# 2. drug-like in ZINC12
```
p.mwt <= 500 and p.mwt >= 150 and p.xlogp <= 5 and p.rb <=7 and p.psa < 150 and p.n_h_donors <= 5 and p.n_h_acceptors <= 10
```

# 3. Download [DUD-E dataset](http://dude.docking.org/)
```bash
  bash 1download_dude.sh
```

# 4. result
result of diverse subsets. 8-fold cross validation.
![EF1](figures/EF1.png)
![AUC](figures/AUC.png)


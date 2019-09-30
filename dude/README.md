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
## 4.1 performance on DUD-E

| subset         | PROP EF1 | PROP AUC | FP EF1 | FP AUC |
| :------------- | -------: | -------: | -----: | -----: |
| D8_DUDE(ism)   |    21.26 |     0.81 |  10.39 |   0.79 |
| D8_DUDE(sdf)   |     4.79 |     0.73 |   5.57 |   0.59 |
| D8_rebulit     |     4.48 |     0.57 |   4.51 |   0.72 |
| Full(3-fold)*  |    11.66 |     0.66 |  14.57 |   0.85 |
| Full(10-fold)^ |    20.50 |     0.78 |      - |      - |

\* 3-fold is clustering cross target cross validation.

^ 10-fold is random cross target cross validation.

## 4.2 distribution

Should use **.sdf** format in DUD-E for **.ism** lake of charge information.

Diverse use **.ism**, lake of charge information
![d8_ism](figures/dist/distribution_diverse_DUDE_ism.jpg)

Diverse use **.sdf**, match better, very less decoys with HeavyAtomMolWt > 500.
![d8_sdf](figures/dist/distribution_diverse_DUDE_sdf.jpg)

Diverse rebuilt, also less decoys with HighAtomMolWt > 500 for ZINC12 has less than 1% of molecules with HighAtomMolWt > 500.
![d8_rebuilt](figures/dist/distribution_diverse_full.jpg)

Full DUD-E use **.ism**
![dude_ism](figures/dist/distribution_DUDE_ism.jpg)

Full DUD-E use **.sdf**
![dude_sdf](figures/dist/distribution_DUDE_sdf.jpg)
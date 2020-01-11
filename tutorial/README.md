# Tutorial
This tutorial will help you to reproduce the result in the paper.
> Test env:   
> CentOS 7  
> Intel(R) Xeon(R) Gold 5117 CPU @ 2.00GHz (14 cores)   
> No GPU


## 1. Setup conda environment
Install [miniconda](https://docs.conda.io/en/latest/miniconda.html#) and create environment use command:
```bash
# cpu
conda create -n test -c deepchem -c rdkit -c conda-forge -c omnia deepchem=2.3.0 tqdm seaborn rdkit scikit-learn parallel openbabel python=3.6
# gpu
conda create -n test -c deepchem -c rdkit -c conda-forge -c omnia deepchem-gpu=2.3.0 tqdm seaborn rdkit scikit-learn parallel openbabel python=3.6
conda activate test

# conda env export > environment.yml
```

## 2. Random forest on DUD-E

### 2.1 Download DUD-E
```bash
cd RF_DUD-E
bash 1download_dude.sh
```

### 2.2 Fast test by training models using properties only.   
It will convert mol2_files into fingerprints at first time.
```bash
# DUD-E|Cross-Class CV
python cross_target_RF.py -f family3fold.json --prop_only -o result/prop_only
# ignoring rdkit warnings
# output: {'prop': {'ROC': 0.6905003110510383, 'EF1': 13.47612326096607}}
```

### 2.3 Train models on fingerprints (about 2 hours for 1 job).
```bash
# DUD-E|Cross-Class CV
python cross_target_RF.py -o result/family3fold.origin -f family3fold.json

# DUD-E|Random CV
# read target names from family3fold.json and then random split into 3 folds.
python cross_target_RF.py -o result/random3fold.origin -f family3fold.json --random_fold

# DUD-E(MW<=500)|Cross-Class CV
python cross_target_RF.py -o result/family3fold.MW500 -f family3fold.json --MW500

# DUD-E(MW<=500)|Random CV
python cross_target_RF.py -o result/random3fold.MW500 -f family3fold.json --MW500 --random_fold 
```

## 3. ACNN on PDBbind
### 3.1 Download PDBbind v2015:
```
cd ACNN_PDBbind
bash get_pdbbind.sh
```

### 3.2 Fast test by training models on core set, random split.
Train on `core training subset` (binding complex), test on `core test subset` (binding complex):
```bash
python ACNN.py -result_dir result/binding_core_core_random
# ignoring tensorflow warnings and info.

# output:
# ...
# performance of model best on validation dataset:
# {
#   "train": {
#     "pearson_r2_score": 0.41563833250684845,
#     "mean_absolute_error": 1.3624139727078952
#   },
#   "valid": {
#     "pearson_r2_score": 0.08187588216669629,
#     "mean_absolute_error": 2.250938710664448
#   },
#   "test": {
#     "pearson_r2_score": 0.02674488967403431,
#     "mean_absolute_error": 2.1726850719451907
#   }
# }
# Elapsed time 0:05:36.309268
```

### 3.3 Train on `refined set` (ligand alone, core set removed), test on `core set` (ligand alone):

Remove samples similar to core set based on protein sequence similarity (identity > 40%), 2036 samples left.
```bash
python ACNN.py -component ligand -subset refined -test_core 2015 -result_dir result/ligand_refined_core_sequence -clust_file similarity/protein_sequence/v2015.simi0.4.clust.json

# output:
# ...
# ... 18 epoch
# ...
# performance of model best on validation dataset:
# {
#   "train": {
#     "pearson_r2_score": 0.693115165450781,
#     "mean_absolute_error": 1.532309637413275
#   },
#   "valid": {
#     "pearson_r2_score": 0.5149077363761936,
#     "mean_absolute_error": 1.8515285022586003
#   },
#   "test": {
#     "pearson_r2_score": 0.7576200667942937,
#     "mean_absolute_error": 1.7107928460928112
#   }
# }
# Elapsed time 0:53:55.791942
```

Remove samples similar to core set based on ligand scaffold similarity (scaffold Tc > 0.8).   
2049 samples left, subsample to 2036 samples.
```bash
python ACNN.py -component ligand -subset refined -test_core 2015 -result_dir result/ligand_refined_core_scaffold -clust_file similarity/ligand_scaffold/v2015.simi0.8.clust.json -reload
# use `-reload` to reuse the featurized dataset   
# featurized: convert pdb files into atomic numbers, coordinates, and neighbor-list.

# output:
# ...
# ... 23 epoch
# ...
# performance of model best on validation dataset:
# {
#   "train": {
#     "pearson_r2_score": 0.8418531080338748,
#     "mean_absolute_error": 0.9304247710538223
#   },
#   "valid": {
#     "pearson_r2_score": 0.8314677703696665,
#     "mean_absolute_error": 1.3423126808802286
#   },
#   "test": {
#     "pearson_r2_score": 0.767228342273758,
#     "mean_absolute_error": 0.826871400588598
#   }
# }
# Elapsed time 1:18:57.334186
```

For comparison, we randomly subsample the refined set to same size (train + valid = 2036 samples).
```bash
python ACNN.py -component ligand -subset refined -test_core 2015 -result_dir result/ligand_refined_core_random -reload

# output:
# ...
# ... 13 epoch
# ...
# performance of model best on validation dataset:
# {
#   "train": {
#     "pearson_r2_score": 0.7764282043985961,
#     "mean_absolute_error": 1.3956259161884608
#   },
#   "valid": {
#     "pearson_r2_score": 0.6838086693503531,
#     "mean_absolute_error": 1.5956438344132664
#   },
#   "test": {
#     "pearson_r2_score": 0.7416788840309323,
#     "mean_absolute_error": 1.4271982015952085
#   }
# }
# Elapsed time 0:38:43.298831
```

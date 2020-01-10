## Setup environment
install [miniconda](https://docs.conda.io/en/latest/miniconda.html#) and create environment use command:
```bash
conda create -n test -c deepchem -c rdkit -c conda-forge -c omnia deepchem=2.3.0 tqdm seaborn rdkit scikit-learn python=3.6
conda activate test
# conda env export > environment.yml
```

## Random forest on DUD-E

Download DUD-E
```bash
cd RF_DUD-E
bash 1download_dude.sh
```

Fast test env by train models using properties only.   
It will convert mol2_files into fingerprints at first time.
```bash
python cross_target_RF.py -f family3fold.json --prop_only -o result/prop_only
# ignoring rdkit warnings
# output: {'prop': {'ROC': 0.6905003110510383, 'EF1': 13.47612326096607}}
```

Train models on fingerprints (about 2 hours for 1 job in my machine).
```bash
python cross_target_RF.py -o result/random3fold.origin -f family3fold.json --random_fold

python cross_target_RF.py -o result/family3fold.origin -f family3fold.json

python cross_target_RF.py -o result/random3fold.MW500 -f family3fold.json --MW500 --random_fold 

python cross_target_RF.py -o result/family3fold.MW500 -f family3fold.json --MW500
```

## ACNN on PDBbind

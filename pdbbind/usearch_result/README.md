# 1. convert the chain including pocket in protein from pdb to fasta
```bash
python3.6 pdbbind2fasta.py -i INDEX_core_data.2013 -d path_to_pdbbind
# output INDEX_core_data.2013.fasta 
```
# 2. cluster by [usearch](https://www.drive5.com/usearch/manual/cmd_cluster_fast.html)

```
./usearch11.0.667_i86linux32 -cluster_fast INDEX_core_data.2013.fasta -id 0.8 -uc INDEX_core_data.2013.uc
```

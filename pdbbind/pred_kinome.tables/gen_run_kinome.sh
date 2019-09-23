exe=/pubhome/jcyang/git/deepchem/examples/pdbbind/pdbbind_atomic_conv.py
data=/pubhome/cshen/docus/dataset/1-raw
save_dir=/pubhome/jcyang/tmp/withHs
ligand_split=/pubhome/jcyang/git/deepchem/examples/pdbbind/similarity/ligand_fingerprint/result
protein_split=/pubhome/jcyang/git/deepchem/examples/pdbbind/similarity/protein_sequence/result
for version in 2015 2018; do
    for subset in core refined general_PL; do
        for component in binding protein ligand; do
            for split in random sequence fingerprint; do
                if [ $split = random ]; then
                    clust=''
                elif [ $split = sequence ]; then
                    clust="-clust_file $protein_split/v${version}.simi0.5.clust.json"
                    continue
                else
                    clust="-clust_file $ligand_split/v${version}.simi0.5.clust.json"
                    continue
                fi
                for seed in 111 222 333 444 555; do
                    dir=job.$version.$subset.$component.$split.$seed
                    test -d $dir || mkdir $dir
                    echo "# export CUDA_VISIBLE_DEVICES=0" >$dir/run.sh
                    echo "source /pubhome/jcyang/anaconda3/bin/activate" >>$dir/run.sh
                    echo "conda activate deepchem" >>$dir/run.sh
                    echo 'export PYTHONPATH=/pubhome/jcyang/git/deepchem:$PYTHONPATH' >>$dir/run.sh
                    echo "python $exe -data_dir $data -save_dir $save_dir \\" >>$dir/run.sh
                    echo "-version $version -subset $subset -reload -split $split -test_kinome \\" >>$dir/run.sh
                    echo "-component $component -seed $seed $clust" >>$dir/run.sh
                done
            done
        done
    done
done
# export CUDA_VISIBLE_DEVICES=0
# source /pubhome/jcyang/anaconda3/bin/activate
# conda activate deepchem
# export PYTHONPATH=/pubhome/jcyang/git/deepchem:$PYTHONPATH
# --data_dir /pubhome/cshen/docus/dataset/1-raw

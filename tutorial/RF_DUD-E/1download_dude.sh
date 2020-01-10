#!/bin/bash
wget http://dude.docking.org/db/subsets/all/all.tar.gz
tar -xvf all.tar.gz
# missing fgfr1 decoys in all.tar.gz
cd all/fgfr1
wget wget http://dude.docking.org/targets/fgfr1/decoys_final.mol2.gz
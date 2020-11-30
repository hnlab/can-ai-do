#!/bin/bash
wget http://dude.docking.org/db/subsets/all/all.tar.gz
tar -xvf all.tar.gz

# missing fgfr1 decoys in all.tar.gz
cd all
wget http://dude.docking.org/targets/fgfr1/fgfr1.tar.gz
tar -xf fgfr1.tar.gz

cd fgfr1/
rm crystal_ligand.mol2
wget http://dude.docking.org/targets/fgfr1/docking/xtal-lig.mol2 -O crystal_ligand.mol2
rm receptor.pdb
wget http://dude.docking.org/targets/fgfr1/docking/grids/rec.crg -O receptor.pdb

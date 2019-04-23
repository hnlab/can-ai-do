"""
Two deepchem splitter for PDBbind database.
"""
from __future__ import division
from __future__ import unicode_literals

import random

__author__ = "Bharath Ramsundar, Aneesh Pappu "
__copyright__ = "Copyright 2016, Stanford University"
__license__ = "MIT"

import tempfile
import numpy as np
import pandas as pd
import itertools
import os
import deepchem as dc
from deepchem.data import DiskDataset
from deepchem.data import NumpyDataset
from deepchem.splits import Splitter


class SequenceSplitter(Splitter):
    """
    Class for doing data splits based on clustering of protein sequence.
    Need uclust file from UCLUST

    O(N**2) algorithm
  """

    def split(self,
              dataset,
              uclust=None,
              seed=None,
              frac_train=.8,
              frac_valid=.1,
              frac_test=.1,
              log_every_n=1000):
        """
        Splits proteins into train/validation/test by sequence clustering.
    """
        # load uclust file
        all_clust_nums = []
        labels = []
        with open(uclust) as f:
            for line in f:
                fields = line.split()
                all_clust_nums.append(int(fields[1]))
                labels.append(fields[8])

        # cluster index of dataset ids
        ids = dataset.ids
        inds_clust = {}
        for i, code in enumerate(ids):
            try:
                ind = labels.index(code)
                num = all_clust_nums[ind]
            except ValueError:
                print("Warning: {} not in clust file, asign it to clust -1".
                      format(code))
                num = -1
                continue
            finally:
                if num not in inds_clust:
                    inds_clust[num] = []
                inds_clust[num].append(i)

        # re numbering cluster by size
        inds = []
        for clust in sorted(inds_clust.values(), key=len, reverse=True):
            inds.extend(clust)

        np.testing.assert_almost_equal(frac_train + frac_valid + frac_test, 1.)
        data_len = len(dataset)
        train_cutoff = int(frac_train * data_len)
        valid_cutoff = int((frac_train + frac_valid) * data_len)
        train_inds = inds[:train_cutoff]
        valid_inds = inds[train_cutoff:valid_cutoff]
        test_inds = inds[valid_cutoff:]
        return train_inds, valid_inds, test_inds


def ClusterFps(fps, cutoff=0.2):
    # (ytz): this is directly copypasta'd from Greg Landrum's clustering example.
    dists = []
    nfps = len(fps)
    from rdkit import DataStructs
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])
    from rdkit.ML.Cluster import Butina
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs


class ButinaSplitter(Splitter):
    """
    Class for doing data splits based on the butina clustering of a bulk tanimoto
    fingerprint matrix.
    """

    def split(self,
              dataset,
              seed=None,
              frac_train=None,
              frac_valid=None,
              frac_test=None,
              log_every_n=1000,
              cutoff=0.2,
              pdbbind_path=None):
        """
        Splits internal compounds into train and validation based on the butina
        clustering algorithm. This splitting algorithm has an O(N^2) run time, where N
        is the number of elements in the dataset. The dataset is expected to be a classification
        dataset.

        This algorithm is designed to generate validation data that are novel chemotypes.

        Note that this function entirely disregards the ratios for frac_train, frac_valid,
        and frac_test. Furthermore, it does not generate a test set, only a train and valid set.

        Setting a small cutoff value will generate smaller, finer clusters of high similarity,
        whereas setting a large cutoff value will generate larger, coarser clusters of low similarity.
        """
        print("Performing butina clustering with cutoff of", cutoff)
        mols = []
        inds = []
        from rdkit import Chem
        from pathlib import Path
        if pdbbind_path is not None:
            pdbbind_path = Path(pdbbind_path)

        # Morgan Fingerprint need sanitize=True, so use ligand.pdb coverted by babel from ligand.mol2
        # http://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints

        # FingerprintMols.FingerprintMol is RDKFingerprint, not error when sanitize=False, but the cluster results are very different to Morgan Fingerprint.
        # http://www.rdkit.org/docs/GettingStartedInPython.html#topological-fingerprints
        for ind, _id in enumerate(dataset.ids):
            if pdbbind_path:
                pdb = pdbbind_path / _id / (_id + '_ligand.pdb')
                mol = Chem.MolFromPDBFile(str(pdb))
                # sdf = pdbbind_path / _id / (_id + '_ligand.sdf')
                # mol = next(Chem.SDMolSupplier(str(sdf), sanitize=False))
            else:
                mol = Chem.MolFromSmiles(_id)
            if mol is None:
                print(
                    "Warning: rdkit fail to load mol {}, will assign it to training set"
                    .format(_id))
                inds.append(ind)
                continue
            mols.append(mol)
        n_mols = len(mols)
        from rdkit.Chem import AllChem
        fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024) for x in mols]

        # from rdkit.Chem.Fingerprints import FingerprintMols
        # fps = [FingerprintMols.FingerprintMol(x) for x in mols]

        scaffold_sets = ClusterFps(fps, cutoff=cutoff)
        scaffold_sets = sorted(scaffold_sets, key=lambda x: -len(x))

        for c_idx, cluster in enumerate(scaffold_sets):
            inds.extend(cluster)

        np.testing.assert_almost_equal(frac_train + frac_valid + frac_test, 1.)
        data_len = len(dataset)
        train_cutoff = int(frac_train * data_len)
        valid_cutoff = int((frac_train + frac_valid) * data_len)
        train_inds = inds[:train_cutoff]
        valid_inds = inds[train_cutoff:valid_cutoff]
        test_inds = inds[valid_cutoff:]
        return train_inds, valid_inds, test_inds


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--index', required=True)
    parser.add_argument('-p', '--pdbbind_path', required=True)
    parser.add_argument(
        '-u',
        '--uclust',
        required=True,
        help=
        "uclust output, format: https://www.drive5.com/usearch/manual/opt_uc.html"
    )
    args = parser.parse_args()

    def read_index(index_file):
        codes = []
        pKs = []
        with open(index_file) as f:
            for i in f:
                if i[0] == '#': continue
                code, reso, year, pK, *others = i.split()
                codes.append(code)
                pKs.append(float(pK))
        return codes, pKs

    codes, pKs = read_index(args.index)
    dataset = NumpyDataset(pKs, ids=codes)

    print("Test SequenceSplitter")
    splitter = SequenceSplitter()
    train, valid, test = splitter.train_valid_test_split(dataset,
                                                         uclust=args.uclust)
    # print(train.ids)
    # print(valid.ids)
    # print(test.ids)
    lens = (len(train), len(valid), len(test))
    print("Total:{} train/valid/test: {} {} {}, sum:{}\n".format(
        len(dataset), *lens, sum(lens)))

    print("Test ButinaSplitter")
    splitter = ButinaSplitter()
    train, valid, test = splitter.train_valid_test_split(
        dataset, cutoff=0.2, pdbbind_path=args.pdbbind_path)
    # print(train.ids)
    # print(valid.ids)
    # print(sorted(test.ids))
    lens = (len(train), len(valid), len(test))
    print("Total:{} train/valid/test: {} {} {}, sum:{}\n".format(
        len(dataset), *lens, sum(lens)))

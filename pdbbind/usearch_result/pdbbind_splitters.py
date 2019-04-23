"""
SequenceSplitter: A deepchem splitter based on clustering of protein sequence.
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

    def __init__(self, uclust_file, *args, **kwargs):
        self.uclust_file = uclust_file
        super(SequenceSplitter, self).__init__(*args, **kwargs)

    def split(self,
              dataset,
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
        with open(self.uclust_file) as f:
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
        valid_cutoff = train_cutoff + int(frac_valid * data_len)
        train_inds = inds[:train_cutoff]
        valid_inds = inds[train_cutoff:valid_cutoff]
        test_inds = inds[valid_cutoff:]
        return train_inds, valid_inds, test_inds


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--index', required=True)
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

    splitter = SequenceSplitter(args.uclust)
    train, valid, test = splitter.train_valid_test_split(dataset)
    print(train.ids)
    print(valid.ids)
    print(test.ids)
    print("total:{} train/valid/test: {} {} {}".format(len(dataset),
                                                       len(train), len(valid),
                                                       len(test)))

"""
PDBBind dataset loader.
"""

from __future__ import division
from __future__ import unicode_literals

import logging
import multiprocessing
import os
import re
import time

import deepchem
import numpy as np
import pandas as pd
import tarfile
from deepchem.feat import rdkit_grid_featurizer as rgf
from deepchem.feat.atomic_coordinates import ComplexNeighborListFragmentAtomicCoordinates
from deepchem.feat.graph_features import AtomicConvFeaturizer
from deepchem.splits import FingerprintSplitter
from deepchem.splits import Splitter

logger = logging.getLogger(__name__)


def featurize_pdbbind(data_dir=None, feat="grid", subset="core"):
  """Featurizes pdbbind according to provided featurization"""
  tasks = ["-logKd/Ki"]
  data_dir = deepchem.utils.get_data_dir()
  pdbbind_dir = os.path.join(data_dir, "pdbbind")
  dataset_dir = os.path.join(pdbbind_dir, "%s_%s" % (subset, feat))

  if not os.path.exists(dataset_dir):
    deepchem.utils.download_url(
        'http://deepchem.io.s3-website-us-west-1.amazonaws.com/featurized_datasets/core_grid.tar.gz'
    )
    deepchem.utils.download_url(
        'http://deepchem.io.s3-website-us-west-1.amazonaws.com/featurized_datasets/full_grid.tar.gz'
    )
    deepchem.utils.download_url(
        'http://deepchem.io.s3-website-us-west-1.amazonaws.com/featurized_datasets/refined_grid.tar.gz'
    )
    if not os.path.exists(pdbbind_dir):
      os.system('mkdir ' + pdbbind_dir)
    deepchem.utils.untargz_file(
        os.path.join(data_dir, 'core_grid.tar.gz'), pdbbind_dir)
    deepchem.utils.untargz_file(
        os.path.join(data_dir, 'full_grid.tar.gz'), pdbbind_dir)
    deepchem.utils.untargz_file(
        os.path.join(data_dir, 'refined_grid.tar.gz'), pdbbind_dir)

  return deepchem.data.DiskDataset(dataset_dir), tasks


def load_pdbbind_grid(split="random",
                      featurizer="grid",
                      subset="core",
                      reload=True):
  """Load PDBBind datasets. Does not do train/test split"""
  if featurizer == 'grid':
    dataset, tasks = featurize_pdbbind(feat=featurizer, subset=subset)

    splitters = {
        'index': deepchem.splits.IndexSplitter(),
        'random': deepchem.splits.RandomSplitter(),
        'time': deepchem.splits.TimeSplitterPDBbind(dataset.ids)
    }
    splitter = splitters[split]
    train, valid, test = splitter.train_valid_test_split(dataset)

    transformers = []
    for transformer in transformers:
      train = transformer.transform(train)
    for transformer in transformers:
      valid = transformer.transform(valid)
    for transformer in transformers:
      test = transformer.transform(test)

    all_dataset = (train, valid, test)
    return tasks, all_dataset, transformers

  else:
    data_dir = deepchem.utils.get_data_dir()
    if reload:
      save_dir = os.path.join(
          data_dir, "pdbbind_" + subset + "/" + featurizer + "/" + str(split))

    dataset_file = os.path.join(data_dir, subset + "_smiles_labels.csv")

    if not os.path.exists(dataset_file):
      deepchem.utils.download_url(
          'http://deepchem.io.s3-website-us-west-1.amazonaws.com/datasets/' +
          subset + "_smiles_labels.csv")

    tasks = ["-logKd/Ki"]
    if reload:
      loaded, all_dataset, transformers = deepchem.utils.save.load_dataset_from_disk(
          save_dir)
      if loaded:
        return tasks, all_dataset, transformers

    if featurizer == 'ECFP':
      featurizer = deepchem.feat.CircularFingerprint(size=1024)
    elif featurizer == 'GraphConv':
      featurizer = deepchem.feat.ConvMolFeaturizer()
    elif featurizer == 'Weave':
      featurizer = deepchem.feat.WeaveFeaturizer()
    elif featurizer == 'Raw':
      featurizer = deepchem.feat.RawFeaturizer()

    loader = deepchem.data.CSVLoader(
        tasks=tasks, smiles_field="smiles", featurizer=featurizer)
    dataset = loader.featurize(dataset_file, shard_size=8192)
    transformers = [
        deepchem.trans.NormalizationTransformer(
            transform_y=True, dataset=dataset)
    ]

    for transformer in transformers:
      dataset = transformer.transform(dataset)
    df = pd.read_csv(dataset_file)

    if split == None:
      return tasks, (dataset, None, None), transformers

    splitters = {
        'index': deepchem.splits.IndexSplitter(),
        'random': deepchem.splits.RandomSplitter(),
        'scaffold': deepchem.splits.ScaffoldSplitter(),
        'time': deepchem.splits.TimeSplitterPDBbind(np.array(df['id']))
    }
    splitter = splitters[split]
    train, valid, test = splitter.train_valid_test_split(dataset)

    if reload:
      deepchem.utils.save.save_dataset_to_disk(save_dir, train, valid, test,
                                               transformers)

    return tasks, (train, valid, test), transformers


class FingerprintSplitter4Pdbbind(FingerprintSplitter):
  """
    Class for doing data splits based on the fingerprints of small molecules
    O(N**2) algorithm
    load small moleculer with SDF format
  """

  def __init__(self, pdbbind_data_folder, verbose=False):
    """Provide input information for splits."""
    self.pdbbind_data_folder = pdbbind_data_folder
    self.verbose = verbose
    super(FingerprintSplitter4Pdbbind, self).__init__(verbose)

  def split(self,
            dataset,
            seed=None,
            frac_train=.8,
            frac_valid=.1,
            frac_test=.1,
            log_every_n=1000):
    """
      Splits internal compounds into train/validation/test by fingerprint.
    """
    np.testing.assert_almost_equal(frac_train + frac_valid + frac_test, 1.)
    data_len = len(dataset)
    mols, fingerprints = [], []
    train_inds, valid_inds, test_inds = [], [], []
    from rdkit import Chem
    from rdkit.Chem.Fingerprints import FingerprintMols
    # for ind, smiles in enumerate(dataset.ids):
    for ind, pdb in enumerate(dataset.ids):
      ligand_file = os.path.join(self.pdbbind_data_folder, pdb,
                                 "%s_ligand.sdf" % pdb)
      suppl = Chem.SDMolSupplier(str(ligand_file), sanitize=False)
      mol = suppl[0]
      # mol = Chem.MolFromSmiles(smiles, sanitize=False)
      mols.append(mol)
      fp = FingerprintMols.FingerprintMol(mol)
      fingerprints.append(fp)

    distances = np.ones(shape=(data_len, data_len))
    from rdkit import DataStructs
    time1 = time.time()
    for i in range(data_len):
      for j in range(data_len):
        distances[i][j] = 1 - DataStructs.FingerprintSimilarity(
            fingerprints[i], fingerprints[j])
    time2 = time.time()
    print("[%s] Splitting based on Fingerprint of ligands took %0.3f s\n" %
          (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), time2 - time1))

    train_cutoff = int(frac_train * len(dataset))
    valid_cutoff = int(frac_valid * len(dataset))

    # Pick the mol closest to everything as the first element of training
    closest_ligand = np.argmin(np.sum(distances, axis=1))
    train_inds.append(closest_ligand)
    cur_distances = [float('inf')] * data_len
    self.update_distances(closest_ligand, cur_distances, distances, train_inds)
    for i in range(1, train_cutoff):
      closest_ligand = np.argmin(cur_distances)
      train_inds.append(closest_ligand)
      self.update_distances(closest_ligand, cur_distances, distances,
                            train_inds)

    # Pick the closest mol from what is left
    index, best_dist = 0, float('inf')
    for i in range(data_len):
      if i in train_inds:
        continue
      dist = np.sum(distances[i])
      if dist < best_dist:
        index, best_dist = i, dist
    valid_inds.append(index)

    leave_out_indexes = train_inds + valid_inds
    cur_distances = [float('inf')] * data_len
    self.update_distances(index, cur_distances, distances, leave_out_indexes)
    for i in range(1, valid_cutoff):
      closest_ligand = np.argmin(cur_distances)
      valid_inds.append(closest_ligand)
      leave_out_indexes.append(closest_ligand)
      self.update_distances(closest_ligand, cur_distances, distances,
                            leave_out_indexes)

    # Test is everything else
    for i in range(data_len):
      if i in leave_out_indexes:
        continue
      test_inds.append(i)
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


class ButinaSplitter4pdbbind(Splitter):
  """
    Class for doing data splits based on the butina clustering of a bulk tanimoto
    fingerprint matrix.
  """

  def __init__(self, pdbbind_path, reweight=True, *args, **kwargs):
    self.pdbbind_path = pdbbind_path
    self.reweight = reweight
    self.ids_weight = {}
    super(ButinaSplitter4pdbbind, self).__init__(*args, **kwargs)

  def train_valid_test_split(self,
                             dataset,
                             train_dir=None,
                             valid_dir=None,
                             test_dir=None,
                             frac_train=.8,
                             frac_valid=.1,
                             frac_test=.1,
                             seed=None,
                             log_every_n=1000,
                             verbose=True,
                             **kwargs):
    """
      Splits self into train/validation/test sets.

      Returns Dataset objects.
    """
    train_inds, valid_inds, test_inds = self.split(
        dataset,
        seed=seed,
        frac_train=frac_train,
        frac_test=frac_test,
        frac_valid=frac_valid,
        log_every_n=log_every_n,
        **kwargs)
    import tempfile
    if train_dir is None:
      train_dir = tempfile.mkdtemp()
    if valid_dir is None:
      valid_dir = tempfile.mkdtemp()
    if test_dir is None:
      test_dir = tempfile.mkdtemp()

    def reweight(dataset):

      def fn(x, y, w, ids):
        for i, _id in enumerate(ids):
          w[i][0] = np.float64(self.ids_weight[_id])
        return x, y, w, ids

      return dataset.transform(fn)

    if self.reweight:
      dataset = reweight(dataset)

    train_dataset = dataset.select(train_inds, train_dir)
    if frac_valid != 0:
      valid_dataset = dataset.select(valid_inds, valid_dir)
    else:
      valid_dataset = None
    test_dataset = dataset.select(test_inds, test_dir)

    return train_dataset, valid_dataset, test_dataset

  def split(self,
            dataset,
            seed=None,
            frac_train=None,
            frac_valid=None,
            frac_test=None,
            log_every_n=1000,
            cutoff=0.2):
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
    from rdkit import Chem
    from pathlib import Path
    pdbbind_path = Path(self.pdbbind_path)

    mols = []
    inds_for_split = []
    inds_ids = {}
    inds_in_mols_to_dataset = {}

    for ind, _id in enumerate(dataset.ids):
      inds_ids[ind] = _id

      # FingerprintMols.FingerprintMol is RDKFingerprint, not error when sanitize=False, but the cluster
      # results are very different to Morgan Fingerprint.
      # http://www.rdkit.org/docs/GettingStartedInPython.html#topological-fingerprints
      # sdf = pdbbind_path / _id / (_id + '_ligand.sdf')
      # mol = next(Chem.SDMolSupplier(str(sdf), sanitize=False))

      # Morgan Fingerprint need sanitize=True, so use ligand.pdb coverted by babel from ligand.mol2
      # http://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints
      pdb = pdbbind_path / _id / (_id + '_ligand.pdb')
      mol = Chem.MolFromPDBFile(str(pdb))

      if mol is None:
        print(
            "WARNING: RDKit failed to load ligand of {}, assign it to training dataset."
            .format(_id))
        inds_for_split.append(ind)
        self.ids_weight[_id] = 1.0
        continue
      mols.append(mol)
      inds_in_mols_to_dataset[len(mols) - 1] = ind

    # from rdkit.Chem.Fingerprints import FingerprintMols
    # fps = [FingerprintMols.FingerprintMol(x) for x in mols]

    from rdkit.Chem import AllChem
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024) for x in mols]

    scaffold_sets = list(ClusterFps(fps, cutoff=cutoff))
    np.random.seed(seed)
    np.random.shuffle(scaffold_sets)
    scaffold_sets = sorted(scaffold_sets, key=lambda x: -len(x))

    scaffold_sets_inds_in_dataset = [[
        inds_in_mols_to_dataset[ind_in_mols] for ind_in_mols in cluster
    ] for cluster in scaffold_sets]

    for cluster in scaffold_sets_inds_in_dataset:
      inds_for_split.extend(cluster)
      for ind in cluster:
        self.ids_weight[inds_ids[ind]] = 1.0 / len(cluster)

    inds = inds_for_split

    np.testing.assert_almost_equal(frac_train + frac_valid + frac_test, 1.)
    data_len = len(dataset)
    train_cutoff = int(frac_train * data_len)
    valid_cutoff = int((frac_train + frac_valid) * data_len)
    train_inds = inds[:train_cutoff]
    valid_inds = inds[train_cutoff:valid_cutoff]
    test_inds = inds[valid_cutoff:]
    return train_inds, valid_inds, test_inds


class ScaffoldSplitter4pdbbind(Splitter):
  """
    Class for doing data splits based on the butina clustering of a bulk tanimoto
    fingerprint matrix.
  """

  def __init__(self, pdbbind_path, reweight=True, *args, **kwargs):
    self.pdbbind_path = pdbbind_path
    self.reweight = reweight
    self.ids_weight = {}
    super(ScaffoldSplitter4pdbbind, self).__init__(*args, **kwargs)

  def train_valid_test_split(self,
                             dataset,
                             train_dir=None,
                             valid_dir=None,
                             test_dir=None,
                             frac_train=.8,
                             frac_valid=.1,
                             frac_test=.1,
                             seed=None,
                             log_every_n=1000,
                             verbose=True,
                             **kwargs):
    """
      Splits self into train/validation/test sets.

      Returns Dataset objects.
    """
    train_inds, valid_inds, test_inds = self.split(
        dataset,
        seed=seed,
        frac_train=frac_train,
        frac_test=frac_test,
        frac_valid=frac_valid,
        log_every_n=log_every_n,
        **kwargs)
    import tempfile
    if train_dir is None:
      train_dir = tempfile.mkdtemp()
    if valid_dir is None:
      valid_dir = tempfile.mkdtemp()
    if test_dir is None:
      test_dir = tempfile.mkdtemp()

    def reweight(dataset):

      def fn(x, y, w, ids):
        for i, _id in enumerate(ids):
          w[i][0] = np.float64(self.ids_weight[_id])
        return x, y, w, ids

      return dataset.transform(fn)

    if self.reweight:
      dataset = reweight(dataset)

    train_dataset = dataset.select(train_inds, train_dir)
    if frac_valid != 0:
      valid_dataset = dataset.select(valid_inds, valid_dir)
    else:
      valid_dataset = None
    test_dataset = dataset.select(test_inds, test_dir)

    return train_dataset, valid_dataset, test_dataset

  def split(self,
            dataset,
            seed=None,
            frac_train=None,
            frac_valid=None,
            frac_test=None,
            log_every_n=1000,
            cutoff=0.2):
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
    from rdkit import Chem
    from pathlib import Path
    from rdkit.Chem.Scaffolds import MurckoScaffold
    pdbbind_path = Path(self.pdbbind_path)

    mols = []
    inds_for_split = []
    inds_ids = {}
    inds_in_mols_to_dataset = {}

    for ind, _id in enumerate(dataset.ids):
      inds_ids[ind] = _id

      # FingerprintMols.FingerprintMol is RDKFingerprint, not error when sanitize=False, but the cluster
      # results are very different to Morgan Fingerprint.
      # http://www.rdkit.org/docs/GettingStartedInPython.html#topological-fingerprints
      # sdf = pdbbind_path / _id / (_id + '_ligand.sdf')
      # mol = next(Chem.SDMolSupplier(str(sdf), sanitize=False))

      # Morgan Fingerprint need sanitize=True, so use ligand.pdb coverted by babel from ligand.mol2
      # http://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints
      pdb = pdbbind_path / _id / (_id + '_ligand.pdb')
      mol = Chem.MolFromPDBFile(str(pdb))

      if mol is None:
        print(
            "WARNING: RDKit failed to load ligand of {}, assign it to training dataset."
            .format(_id))
        inds_for_split.append(ind)
        self.ids_weight[_id] = 1.0
        continue
      # use scaffold for clustering.
      # smiles = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True)
      mol = MurckoScaffold.GetScaffoldForMol(mol)
      mols.append(mol)
      inds_in_mols_to_dataset[len(mols) - 1] = ind

    # from rdkit.Chem.Fingerprints import FingerprintMols
    # fps = [FingerprintMols.FingerprintMol(x) for x in mols]

    from rdkit.Chem import AllChem
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024) for x in mols]

    scaffold_sets = list(ClusterFps(fps, cutoff=cutoff))
    np.random.seed(seed)
    np.random.shuffle(scaffold_sets)
    scaffold_sets = sorted(scaffold_sets, key=lambda x: -len(x))

    scaffold_sets_inds_in_dataset = [[
        inds_in_mols_to_dataset[ind_in_mols] for ind_in_mols in cluster
    ] for cluster in scaffold_sets]

    for cluster in scaffold_sets_inds_in_dataset:
      inds_for_split.extend(cluster)
      for ind in cluster:
        self.ids_weight[inds_ids[ind]] = 1.0 / len(cluster)

    inds = inds_for_split

    np.testing.assert_almost_equal(frac_train + frac_valid + frac_test, 1.)
    data_len = len(dataset)
    train_cutoff = int(frac_train * data_len)
    valid_cutoff = int((frac_train + frac_valid) * data_len)
    train_inds = inds[:train_cutoff]
    valid_inds = inds[train_cutoff:valid_cutoff]
    test_inds = inds[valid_cutoff:]
    return train_inds, valid_inds, test_inds


class SequenceSplitter(Splitter):
  """
    Class for doing data splits based on clustering of protein sequence.
    Need uclust file from UCLUST

    O(N**2) algorithm
  """

  def __init__(self, uclust_file, reweight=True, *args, **kwargs):
    self.uclust_file = uclust_file
    self.reweight = reweight
    self.ids_weight = {}
    super(SequenceSplitter, self).__init__(*args, **kwargs)

  def train_valid_test_split(self,
                             dataset,
                             train_dir=None,
                             valid_dir=None,
                             test_dir=None,
                             frac_train=.8,
                             frac_valid=.1,
                             frac_test=.1,
                             seed=None,
                             log_every_n=1000,
                             verbose=True,
                             **kwargs):
    """
      Splits self into train/validation/test sets.

      Returns Dataset objects.
      """
    train_inds, valid_inds, test_inds = self.split(
        dataset,
        seed=seed,
        frac_train=frac_train,
        frac_test=frac_test,
        frac_valid=frac_valid,
        log_every_n=log_every_n,
        **kwargs)
    import tempfile
    if train_dir is None:
      train_dir = tempfile.mkdtemp()
    if valid_dir is None:
      valid_dir = tempfile.mkdtemp()
    if test_dir is None:
      test_dir = tempfile.mkdtemp()

    def reweight(dataset):
      print(self.ids_weight)

      def fn(x, y, w, ids):
        for i, _id in enumerate(ids):
          w[i][0] = np.float64(self.ids_weight[_id])
        return x, y, w, ids

      return dataset.transform(fn)

    if self.reweight:
      dataset = reweight(dataset)

    train_dataset = dataset.select(train_inds, train_dir)
    if frac_valid != 0:
      valid_dataset = dataset.select(valid_inds, valid_dir)
    else:
      valid_dataset = None
    test_dataset = dataset.select(test_inds, test_dir)

    return train_dataset, valid_dataset, test_dataset

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
    uc_cluster_inds, uc_labels = [], []
    with open(self.uclust_file) as f:
      for line in f:
        # more clear even use uc_labels.index()
        if line[0] != "C":
          fields = line.split()
          uc_cluster_inds.append(int(fields[1]))
          uc_labels.append(fields[8])
        else:
          break

    # cluster index of dataset ids
    cluster_inds_dataset_inds = {}
    cluster_inds_dataset_ids = {}
    for dataset_ind, _id in enumerate(dataset.ids):
      try:
        cluster_ind = uc_cluster_inds[uc_labels.index(_id)]
      except ValueError:
        print("WARNING: {} not in cluster file, assign it to cluster -1".format(
            _id))
        cluster_ind = -1
        continue
      finally:
        if cluster_ind not in cluster_inds_dataset_inds:
          cluster_inds_dataset_inds[cluster_ind] = []
        if cluster_ind not in cluster_inds_dataset_ids:
          cluster_inds_dataset_ids[cluster_ind] = []
        cluster_inds_dataset_inds[cluster_ind].append(dataset_ind)
        cluster_inds_dataset_ids[cluster_ind].append(_id)

    for cluster_ind, dataset_ids in cluster_inds_dataset_ids.items():
      for _id in dataset_ids:
        self.ids_weight[_id] = 1.0 / len(cluster_inds_dataset_ids[cluster_ind])

    # re numbering cluster by size
    inds = []
    np.random.seed(seed)
    clusters = list(cluster_inds_dataset_inds.values())
    np.random.shuffle(clusters)
    for dataset_inds in sorted(clusters, key=len, reverse=True):
      inds.extend(dataset_inds)

    np.testing.assert_almost_equal(frac_train + frac_valid + frac_test, 1.)
    data_len = len(dataset)
    train_cutoff = int(frac_train * data_len)
    valid_cutoff = int((frac_train + frac_valid) * data_len)
    train_inds = inds[:train_cutoff]
    valid_inds = inds[train_cutoff:valid_cutoff]
    test_inds = inds[valid_cutoff:]
    return train_inds, valid_inds, test_inds


class PocketSplitter(Splitter):
  """
    Class for doing data splits based on clustering of protein sequence.
    Need uclust file from UCLUST

    O(N**2) algorithm
  """

  def __init__(self, clust_file, reweight=True, *args, **kwargs):
    self.clust_file = clust_file
    self.reweight = reweight
    self.ids_weight = {}
    super(PocketSplitter, self).__init__(*args, **kwargs)

  def train_valid_test_split(self,
                             dataset,
                             train_dir=None,
                             valid_dir=None,
                             test_dir=None,
                             frac_train=.8,
                             frac_valid=.1,
                             frac_test=.1,
                             seed=None,
                             log_every_n=1000,
                             verbose=True,
                             **kwargs):
    """
      Splits self into train/validation/test sets.

      Returns Dataset objects.
      """
    train_inds, valid_inds, test_inds = self.split(
        dataset,
        seed=seed,
        frac_train=frac_train,
        frac_test=frac_test,
        frac_valid=frac_valid,
        log_every_n=log_every_n,
        **kwargs)
    import tempfile
    if train_dir is None:
      train_dir = tempfile.mkdtemp()
    if valid_dir is None:
      valid_dir = tempfile.mkdtemp()
    if test_dir is None:
      test_dir = tempfile.mkdtemp()

    def reweight(dataset):
      print(self.ids_weight)

      def fn(x, y, w, ids):
        for i, _id in enumerate(ids):
          w[i][0] = np.float64(self.ids_weight[_id])
        return x, y, w, ids

      return dataset.transform(fn)

    if self.reweight:
      dataset = reweight(dataset)

    train_dataset = dataset.select(train_inds, train_dir)
    if frac_valid != 0:
      valid_dataset = dataset.select(valid_inds, valid_dir)
    else:
      valid_dataset = None
    test_dataset = dataset.select(test_inds, test_dir)

    return train_dataset, valid_dataset, test_dataset

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
    with open(self.clust_file) as f:
      import json
      clust_ids = json.load(f)
    for clust in clust_ids:
      for id_ in clust:
        self.ids_weight[id_] = 1.0 / len(clust)

    dataset_ids = dataset.ids
    weights = np.ones_like(dataset_ids)
    for i, id_ in enumerate(dataset_ids):
      if id_ in self.ids_weight:
        weights[i] = self.ids_weight[id_]
      else:
        self.ids_weight[id_] = 1.0

    # shuffle for not stable sort
    np.random.seed(seed)
    shuff_ids = np.random.permutation(len(weights))
    shuff_ws = weights[shuff_ids]
    # sort index by weight.
    inds = shuff_ids[np.argsort(shuff_ws)]

    np.testing.assert_almost_equal(frac_train + frac_valid + frac_test, 1.)
    data_len = len(dataset)
    train_cutoff = int(frac_train * data_len)
    valid_cutoff = int((frac_train + frac_valid) * data_len)
    train_inds = inds[:train_cutoff]
    valid_inds = inds[train_cutoff:valid_cutoff]
    test_inds = inds[valid_cutoff:]
    return train_inds, valid_inds, test_inds


def load_pdbbind(reload=True,
                 data_dir=None,
                 version="2015",
                 subset="core",
                 shard_size=4096,
                 load_binding_pocket=False,
                 featurizer="grid",
                 split="random",
                 split_complex=False,
                 split_seed=None,
                 clust_file=None,
                 reweight=True,
                 save_dir=None,
                 transform=False,
                 same_protein=False,
                 same_ligand=False,
                 save_timestamp=False):
  """Load raw PDBBind dataset by featurization and split.

  Parameters
  ----------
  reload: Bool, optional
    Reload saved featurized and splitted dataset or not.
  data_dir: Str, optional
    Specifies the directory storing the raw dataset.
  load_binding_pocket: Bool, optional
    Load binding pocket or full protein.
  subset: Str
    Specifies which subset of PDBBind, only "core" or "refined" for now.
  shard_size: Int, optinal
    Specifies size of shards when load general_PL subset considering its large scale.
    Default values are None for core/refined subset and 4096 for general_PL subset.
  featurizer: Str
    Either "grid" or "atomic" for grid and atomic featurizations.
  split: Str
    Either one of "random", "index", "fp", "mfp" and "seq" for random, index, ligand
    Fingerprints, butina clustering with Morgan Fingerprints of ligands, sequence 
    clustering of proteins splitting.
  split_seed: Int, optional
    Specifies the random seed for splitter.
  save_dir: Str, optional
    Specifies the directory to store the featurized and splitted dataset when
    reload is False. If reload is True, it will load saved dataset inside save_dir. 
  save_timestamp: Bool, optional
    Save featurized and splitted dataset with timestamp or not. Set it as True
    when running similar or same jobs simultaneously on multiple compute nodes.
  """

  pdbbind_tasks = ["-logKd/Ki"]

  deepchem_dir = deepchem.utils.get_data_dir()

  if data_dir == None:
    data_dir = deepchem_dir
  data_folder = os.path.join(data_dir, "pdbbind", "v" + version)

  if save_dir == None:
    save_dir = deepchem_dir
  if load_binding_pocket:
    feat_dir = os.path.join(save_dir, "feat-pdbbind", "v" + version,
                            "protein_pocket-%s-%s" % (subset, featurizer))
  else:
    feat_dir = os.path.join(save_dir, "feat-pdbbind", "v" + version,
                            "full_protein-%s-%s" % (subset, featurizer))

  if save_timestamp:
    feat_dir = "%s-%s-%s" % (feat_dir, time.strftime("%Y%m%d",
                                                     time.localtime()),
                             re.search(r"\.(.*)", str(time.time())).group(1))

  loaded = False
  if split is not None:
    if split_seed:
      split_dir = os.path.join(feat_dir, split + str(split_seed))
    else:
      split_dir = os.path.join(feat_dir, str(split))
    if transform:
      split_dir += '.trans'
    if reload:
      print("\nReloading splitted dataset from:\n%s\n" % split_dir)
      loaded, all_dataset, transformers = deepchem.utils.save.load_dataset_from_disk(
          split_dir)
      if loaded:
        return pdbbind_tasks, all_dataset, transformers
      else:
        print('Fail to reload splitted dataset.')

  if reload and loaded == False:
    print("Reloading featurized dataset:\n%s\n" % feat_dir)
    try:
      dataset = deepchem.data.DiskDataset(feat_dir)
      loaded = True
    except ValueError:
      print('Fail to reload featurized dataset.')

  if loaded == False:
    print('Start to featurize dataset form raw data ...')
    if os.path.exists(data_folder):
      logger.info("PDBBind full dataset already exists.")
    else:
      dataset_file = os.path.join(data_dir, "pdbbind_v2015.tar.gz")
      if not os.path.exists(dataset_file):
        logger.warning(
            "About to download PDBBind full dataset. Large file, 2GB")
        deepchem.utils.download_url(
            'http://deepchem.io.s3-website-us-west-1.amazonaws.com/datasets/' +
            "pdbbind_v2015.tar.gz",
            dest_dir=data_dir)

      print("Untarring full dataset...")
      deepchem.utils.untargz_file(
          dataset_file, dest_dir=os.path.join(data_dir, "pdbbind"))

    print("\nRaw dataset:\n%s" % data_folder)
    print("\nFeaturized dataset:\n%s" % feat_dir)

    if version == "2015":
      if subset == "core":
        index_labels_file = os.path.join(data_folder, "index",
                                         "INDEX_core_data.2013")
      elif subset == "refined":
        index_labels_file = os.path.join(data_folder, "index",
                                         "INDEX_refined_data.2015")
      elif subset == "general_PL":
        index_labels_file = os.path.join(data_folder, "index",
                                         "INDEX_general_PL_data.2015")
      else:
        raise ValueError(
            "%s subsets not supported for version %s" % (subset, version))

    elif version == "2018":
      if subset == "refined":
        index_labels_file = os.path.join(data_folder, "index",
                                         "INDEX_refined_data.2018")
      elif subset == "general_PL":
        index_labels_file = os.path.join(data_folder, "index",
                                         "INDEX_general_PL_data.2018")
      else:
        raise ValueError(
            "%s subsets not supported for version %s" % (subset, version))
    else:
      raise ValueError("version %s not supported" % (version))

    # Extract locations of data
    with open(index_labels_file, "r") as g:
      pdbs = [line[:4] for line in g.readlines() if line[0] != "#"]
    if load_binding_pocket:
      protein_files = [
          os.path.join(data_folder, pdb, "%s_pocket.pdb" % pdb) for pdb in pdbs
      ]
    else:
      protein_files = [
          os.path.join(data_folder, pdb, "%s_protein.pdb" % pdb) for pdb in pdbs
      ]
    if same_protein:
      protein_files = [protein_files[0] for i in protein_files]
    ligand_files = [
        os.path.join(data_folder, pdb, "%s_ligand.sdf" % pdb) for pdb in pdbs
    ]
    if same_ligand:
      ligand_files = [ligand_files[0] for i in ligand_files]
    # Extract labels
    with open(index_labels_file, "r") as g:
      labels = np.array([
          # Lines have format
          # PDB code, resolution, release year, -logKd/Ki, Kd/Ki, reference, ligand name
          # The base-10 logarithm, -log kd/pk
          float(line.split()[3]) for line in g.readlines() if line[0] != "#"
      ])

    # Featurize Data
    if featurizer == "grid":
      featurizer = rgf.RdkitGridFeaturizer(
          voxel_width=2.0,
          feature_types=[
              'ecfp', 'splif', 'hbond', 'salt_bridge', 'pi_stack', 'cation_pi',
              'charge'
          ],
          flatten=True)
    elif featurizer == "atomic" or featurizer == "atomic_conv":
      # Pulled from PDB files. For larger datasets with more PDBs, would use
      # max num atoms instead of exact.
      frag1_num_atoms = 70  # for ligand atoms
      if load_binding_pocket:
        frag2_num_atoms = 1000
        complex_num_atoms = 1070
      else:
        frag2_num_atoms = 24000  # for protein atoms
        complex_num_atoms = 24070  # in total
      max_num_neighbors = 4
      # Cutoff in angstroms
      neighbor_cutoff = 4
      if featurizer == "atomic":
        featurizer = ComplexNeighborListFragmentAtomicCoordinates(
            frag1_num_atoms=frag1_num_atoms,
            frag2_num_atoms=frag2_num_atoms,
            complex_num_atoms=complex_num_atoms,
            max_num_neighbors=max_num_neighbors,
            split_complex=split_complex,
            neighbor_cutoff=neighbor_cutoff)
      if featurizer == "atomic_conv":
        featurizer = AtomicConvFeaturizer(
            labels=labels,
            frag1_num_atoms=frag1_num_atoms,
            frag2_num_atoms=frag2_num_atoms,
            complex_num_atoms=complex_num_atoms,
            neighbor_cutoff=neighbor_cutoff,
            max_num_neighbors=max_num_neighbors,
            batch_size=64)
    else:
      raise ValueError("Featurizer %s not supported" % (featurizer))

    def get_shards(inputs, shard_size):
      if len(inputs) <= shard_size:
        yield inputs
      else:
        assert isinstance(shard_size, int) and 0 < shard_size <= len(inputs)
        print("About to start loading files.\n")
        for shard_ind in range(len(inputs) // shard_size + 1):
          if (shard_ind + 1) * shard_size < len(inputs):
            print("Loading shard %d of size %s." % (shard_ind + 1,
                                                    str(shard_size)))
            yield inputs[shard_ind * shard_size:(shard_ind + 1) * shard_size]
        else:
          print("\nLoading shard %d of size %s." %
                (shard_ind + 1, str(len(inputs) % shard_size)))
          yield inputs[shard_ind * shard_size:len(inputs)]

    def shard_generator(inputs, shard_size):
      for shard_num, shard in enumerate(get_shards(inputs, shard_size)):
        time1 = time.time()
        ligand_files, protein_files, labels, pdbs = zip(*shard)
        features, failures = featurizer.featurize_complexes(
            ligand_files, protein_files)
        labels = np.delete(labels, failures)
        labels = labels.reshape((len(labels), 1))
        weight = np.ones_like(labels)
        ids = np.delete(pdbs, failures)
        assert len(features) == len(labels) == len(weight) == len(ids)
        time2 = time.time()
        print("[%s] Featurizing shard %d took %0.3f s\n" % (time.strftime(
            "%Y-%m-%d %H:%M:%S", time.localtime()), shard_num, time2 - time1))
        yield features, labels, weight, ids

    print(
        "\n[%s] Featurizing and constructing dataset without failing featurization for"
        % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "\"%s\"\n" % data_folder)
    feat_t1 = time.time()
    zipped = list(zip(ligand_files, protein_files, labels, pdbs))
    dataset = deepchem.data.DiskDataset.create_dataset(
        shard_generator(zipped, shard_size),
        data_dir=feat_dir,
        tasks=pdbbind_tasks,
        verbose=True)
    feat_t2 = time.time()
    print("\n[%s] Featurization and construction finished, %0.3f s passed.\n" %
          (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
           feat_t2 - feat_t1))

  # Default: No transformations of data
  if transform:
    transformers = [
        deepchem.trans.NormalizationTransformer(
            transform_y=True, dataset=dataset)
    ]
  else:
    transformers = []
  for transformer in transformers:
    dataset = transformer.transform(dataset)

  # Split dataset
  print("\nSplit dataset...\n")
  if split == None:
    return pdbbind_tasks, (dataset, None, None), transformers

  # TODO(rbharath): This should be modified to contain a cluster split so
  # structures of the same protein aren't in both train/test
  uclust_file_dir = "/pubhome/cshen/docus/projects/can-ai-do/pdbbind/usearch_result"
  if subset == "core" and version == "2015":
    version_suffix = "2013"
  else:
    version_suffix = version
  uclust_file = os.path.join(uclust_file_dir,
                             "INDEX_%s_data.%s.uc" % (subset, version_suffix))
  if split == 'seq' and clust_file:
    uclust_file = clust_file
  splitters = {
      'index': deepchem.splits.IndexSplitter(),
      'random': deepchem.splits.RandomSplitter(),
      'fp': FingerprintSplitter4Pdbbind(data_folder),
      'mfp': ButinaSplitter4pdbbind(data_folder, reweight=reweight),
      'seq': SequenceSplitter(uclust_file, reweight=reweight),
      'scaffold': ScaffoldSplitter4pdbbind(data_folder, reweight=reweight),
      'pocket': PocketSplitter(clust_file, reweight=reweight),
  }
  splitter = splitters[split]
  train, valid, test = splitter.train_valid_test_split(dataset, seed=split_seed)

  all_dataset = (train, valid, test)
  print("\nSaving dataset to \"%s\" ..." % split_dir)
  deepchem.utils.save.save_dataset_to_disk(split_dir, train, valid, test,
                                           transformers)
  return pdbbind_tasks, all_dataset, transformers

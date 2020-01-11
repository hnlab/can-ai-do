"""
Script that trains Atomic Conv models on PDBbind dataset.
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

__author__ = "Bharath Ramsundar"
__copyright__ = "Copyright 2016, Stanford University"
__license__ = "MIT"

import json
import argparse
import numpy as np
import tensorflow as tf
import deepchem as dc

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-max_epoch", type=int, default=100)
parser.add_argument("-patience", type=int, default=3)
parser.add_argument("-version", default='2015')
parser.add_argument("-subset", default='core')
parser.add_argument("-component", default='binding')
parser.add_argument("-split", default='random')
parser.add_argument("-seed", type=int, default=111)
parser.add_argument("-clust_file")
parser.add_argument("-save_dir", default='/tmp')
parser.add_argument("-data_dir")
parser.add_argument("-reload", action='store_true')
parser.add_argument("-trans", action='store_true')
parser.add_argument("-feat_only", action='store_true')
parser.add_argument("-same_protein", action='store_true')
parser.add_argument("-same_ligand", action='store_true')
parser.add_argument("-timestamp", action='store_true')
parser.add_argument("-split_complex", action='store_true')
args = parser.parse_args()

# np seed for split only
np.random.seed(args.seed)
# tf seed not work, every training will different.
tf.set_random_seed(args.seed)

pdbbind_tasks, pdbbind_datasets, transformers = dc.molnet.load_pdbbind(
    reload=args.reload,
    featurizer="atomic",
    version=args.version,
    split=args.split,
    split_seed=args.seed,
    clust_file=args.clust_file,
    split_complex=args.split_complex,
    same_protein=args.same_protein,
    same_ligand=args.same_ligand,
    subset=args.subset,
    load_binding_pocket=True,
    data_dir=args.data_dir,
    save_dir=args.save_dir,
    save_timestamp=args.timestamp,
    transform=args.trans,
)

if args.feat_only:
  raise SystemExit(0)

train_dataset, valid_dataset, test_dataset = pdbbind_datasets

metrics = [
    dc.metrics.Metric(dc.metrics.pearson_r2_score),
    dc.metrics.Metric(dc.metrics.mean_absolute_error)
]

config = tf.ConfigProto()
config.gpu_options.allow_growth = True
batch_size = 16
frag1_num_atoms = 70  # for ligand atoms
# frag2_num_atoms = 24000  # for protein atoms
frag2_num_atoms = 1000  # for pocket atoms
complex_num_atoms = frag1_num_atoms + frag2_num_atoms
model = dc.models.AtomicConvModel(
    batch_size=batch_size,
    frag1_num_atoms=frag1_num_atoms,
    frag2_num_atoms=frag2_num_atoms,
    complex_num_atoms=complex_num_atoms,
    component=args.component,
    configproto=config,
)

train_y = train_dataset.y
valid_y = valid_dataset.y
test_y = test_dataset.y
if args.trans:
  for transformer in reversed(transformers):
    if transformer.transform_y:
      train_y = transformer.untransform(train_y)
      valid_y = transformer.untransform(valid_y)
      test_y = transformer.untransform(test_y)

# Fit trained model
print("Fitting model on train dataset")
patience = 0
best_r2 = 0
best_scores = None
train_evaluator = dc.utils.evaluate.Evaluator(model, train_dataset,
                                              transformers)
valid_evaluator = dc.utils.evaluate.Evaluator(model, valid_dataset,
                                              transformers)
test_evaluator = dc.utils.evaluate.Evaluator(model, test_dataset, transformers)


def copy_checkpoint(source, target='best_checkpoint'):
  import os
  from shutil import copyfile
  dirname = os.path.dirname(source)
  if '/' not in target:
    target = os.path.join(dirname, target)
  for item in os.listdir(dirname):
    item = os.path.join(dirname, item)
    name, ext = os.path.splitext(item)
    if name == source:
      copyfile(source + ext, target + ext)
  return target


best_checkpoint = None
for i in range(args.max_epoch):
  model.fit(train_dataset, nb_epoch=1)

  print("Evaluating model at {} epoch".format(i + 1))
  valid_scores = valid_evaluator.compute_model_performance(metrics)
  print("Validation scores")
  print(valid_scores)

  if valid_scores['pearson_r2_score'] < best_r2:
    patience += 1
    if patience > args.patience:
      break
  else:
    last_checkpoint = model.get_checkpoints()[-1]
    best_checkpoint = copy_checkpoint(last_checkpoint)
    patience = 0
    best_r2 = valid_scores['pearson_r2_score']
    print('### Better at epoch {}\n'.format(i + 1))

model.restore(checkpoint=best_checkpoint)
train_scores = train_evaluator.compute_model_performance(
    metrics, csv_out="train.csv")
valid_scores = valid_evaluator.compute_model_performance(
    metrics, csv_out="valid.csv")
test_scores = test_evaluator.compute_model_performance(
    metrics, csv_out="test.csv")

best_scores = {
    'train': train_scores,
    'valid': valid_scores,
    'test': test_scores
}
print('peformances of model best on validation dataset:')
print(json.dumps(best_scores, indent=2))

with open('best_scores.json', 'w') as f:
  data = vars(args)
  data['best_scores'] = best_scores
  json.dump(data, f, indent=2)

with open('splitted_ids.json', 'w') as f:
  data['splitted_ids'] = {
      'train': list(train_dataset.ids),
      'valid': list(valid_dataset.ids),
      'test': list(test_dataset.ids)
  }
  json.dump(data, f, indent=2)

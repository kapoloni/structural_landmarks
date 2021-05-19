#!/usr/bin/python3
"""
Save classifier with best params
"""

# Standard library imports
import os
import sys
import argparse
import numpy as np
import pandas as pd

# Third party imports
from sklearn import svm
from joblib import dump
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# Local application imports
from utils import create_folders


def parse_args(args):
    """!@brief
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description='Classification models')
    parser.add_argument('--dest_folder', help='Destination folder',
                        default='../experiment', type=str)
    parser.add_argument('--labels', help='Labels to perform the experiment',
                        default='cn_ad', type=str)
    parser.add_argument('--side', help='Hippocampal side',
                        default='L', type=str)
    parser.add_argument('--ft', help='Feature',
                        type=str)
    parser.add_argument('--fold', help='Fold number',
                        type=str)
    parser.add_argument('--c', help='C from the SVM',
                        type=float)
    return parser.parse_args(args)


def svm_rbf(X, y, output):

    clf = svm.SVC(kernel='rbf',
                  gamma='scale',
                  C=float(args.c),
                  class_weight='balanced',
                  probability=True)

    steps = [("scaler", StandardScaler()),
             ("classifier", clf)]

    pipe = Pipeline(steps=steps)
    pipe.fit(X, y)
    dump(pipe, output)


def get_descriptors(X):
    attributes = {}
    desc_size = 96
    # gm
    attributes['gm'] = X.iloc[:, :desc_size].astype(
        np.float64).reset_index(drop=True)
    # wm
    attributes['wm'] = X.iloc[:, desc_size:desc_size * 2].astype(
        np.float64).reset_index(drop=True)
    # csf
    attributes['csf'] = X.iloc[:, desc_size * 2:].astype(
        np.float64).reset_index(drop=True)

    attributes['tissues'] = pd.concat([attributes['gm'],
                                       attributes['wm'],
                                       attributes['csf']],
                                      axis=1).reset_index(drop=True)
    return attributes


def get_data(feature):
    input = os.path.join(final_desc_folder, "dataset_" + args.side + ".csv")
    reader = pd.read_csv(input, sep='\t', header=None)

    attributes = get_descriptors(reader.iloc[:, :-2])

    X = attributes[feature]
    y = reader.iloc[:, -1]

    print(np.unique(y, return_counts=True))

    if args.labels == "mci_ad":
        y = np.where(y == args.labels.split("_")[0].upper(), 1, 0)
    else:
        y = np.where(y == args.labels.split("_")[0].upper(), 0, 1)

    return X, y


def classifier(feature):
    X, y = get_data(feature)
    print(X.shape, y.shape, np.unique(y, return_counts=True))
    output = os.path.join(classifier_folder, "model_" + args.side + ".pkl")

    if not os.path.exists(output):
        svm_rbf(np.array(X), y, output)


if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    # Folder
    basename = os.path.join(args.dest_folder,
                            "hippocampus", args.fold)
    final_desc_folder = os.path.join(basename,
                                     "svm_descriptor_train",
                                     args.labels)
    classifier_folder = os.path.join(basename, "svm_models",
                                     args.ft, args.labels)
    create_folders(classifier_folder)

    classifier(args.ft)

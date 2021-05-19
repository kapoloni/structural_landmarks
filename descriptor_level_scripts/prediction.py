#!/usr/bin/python3
"""
Classifier predictions
"""

# Standard library imports
import os
import sys
import argparse
import numpy as np
import pandas as pd

# Third party imports
from joblib import load

# Local application imports
from utils import (create_folders,
                   get_directories_class,
                   get_position_descriptor_txt)


def parse_args(args):
    """!@brief
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description='Predictions')
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
    return parser.parse_args(args)


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


def count_brains_texture(in_, feature, c1, c2):
    counts = []
    for v in range(0, in_.shape[0]):
        if type(in_.iloc[v]) is str:
            imgs_ = {}
            _, desc = get_position_descriptor_txt(in_.iloc[v])

            attributes = get_descriptors(desc)

            X = attributes[feature]
            y = in_.iloc[v].split("/")[-2]
            if args.labels == "mci_ad":
                y = np.where(y == args.labels.split("_")[0].upper(), 1, 0)
            else:
                y = np.where(y == args.labels.split("_")[0].upper(), 0, 1)

            if X.shape[0] > 0:
                out_file = os.path.join(classifier_folder,
                                        "model_" + args.side + ".pkl")
                p0, p1, pp0, pp1, N = predict_decisionF(X,
                                                        out_file)

                subj = in_.iloc[v].split("/")[-1].split("_" + args.side)[0]
                imgs_['subj'], imgs_['0'], \
                    imgs_['1'], imgs_['class'] = subj, p0, p1, y
                imgs_['p0'], imgs_['p1'] = pp0, pp1
                imgs_['N'] = N
                counts.append(imgs_)

    return pd.DataFrame(counts)


def predict_decisionF(X, train_):

    clf = load(train_)

    result = pd.DataFrame(clf.predict_proba(X), columns=clf.classes_)
    result['pred'] = clf.predict(X)
    result['prob_pred'] = result.apply(lambda x: x[int(x['pred'])], axis=1)

    if args.labels == 'mci_ad':
        n0 = len(result[result['pred'] == 1])
        n1 = len(result[result['pred'] == 0])
        p0 = result[1].sum()
        p1 = result[0].sum()
    else:
        n0 = len(result[result['pred'] == 0])
        n1 = len(result[result['pred'] == 1])
        p0 = result[0].sum()
        p1 = result[1].sum()

    return n0, n1, p0, p1, len(result)


def save_csv_header(data, name):
    data.to_csv(name, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    basename = os.path.join(args.dest_folder,
                            "hippocampus", args.fold)
    classifier_folder = os.path.join(basename, "svm_models",
                                     args.ft, args.labels)
    final_desc_folder = os.path.join(basename,
                                     "svm_descriptor_train",
                                     args.labels)

    for step in ['validation', 'test']:
        descriptors_folder = os.path.join(basename,
                                          "tissues_descriptors_" + step,
                                          args.labels)
        pred_folder = os.path.join(basename,
                                   "predicted_descriptors_" + step,
                                   args.ft, args.labels)
        create_folders(pred_folder)

        images = os.path.join(args.dest_folder, "splits",
                              step + "_" + args.fold + ".csv")

        filenames = get_directories_class(images,
                                          args.labels.split("_")[0].upper(),
                                          args.labels.split("_")[1].upper())
        # names = filenames.str.split("/", expand=True)[1]

        input = descriptors_folder + "/" + filenames + "_" + \
                                     args.side + "_landmarks.txt"

        # if not os.path.exists(pred_folder + "/data"+_ss+".csv"):
        output = count_brains_texture(input, args.ft,
                                      args.labels.split("_")[0].upper(),
                                      args.labels.split("_")[1].upper())
        save_csv_header(output, os.path.join(pred_folder,
                                             "data_" + args.side + ".csv"))

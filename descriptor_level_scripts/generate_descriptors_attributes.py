#!/usr/bin/python3
"""
    Generate all descriptors attributes
"""
# Standard library imports
import argparse
import numpy as np
import pandas as pd
import sys
import os

# Local application imports
from utils import (get_directories_class,
                   save_csv,
                   create_folders,
                   get_position_descriptor_txt,
                   save_csv_desc,
                   get_position_descriptor)

# Third party imports
import multiprocessing


def parse_args(args):
    """!@brief
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description='Generate descriptors')
    parser.add_argument('--dest_folder', help='Destination folder',
                        default='../experiment', type=str)
    parser.add_argument('--labels', help='Labels to perform the experiment',
                        default='cn_ad', type=str)
    parser.add_argument('--side', help='Hippocampal side',
                        default='L', type=str)
    parser.add_argument('--fold', help='Fold number',
                        type=str)
    return parser.parse_args(args)


def match_structural(params_):
    input_land, output, structural_atlas = params_
    match_in = os.path.join("..", "bin", "bin",
                            "match-landmarks-maxd-prob")

    atlas = os.path.join(structural_atlas, "structural" + args.side)
    create_folders('/'.join(output.split("/")[:-1]))
    if not os.path.exists(output + "_" + args.side + ".csv"):
        os.system(match_in + " -io " + atlas + " " + input_land +
                  "_" + args.side + " -to 0.5 > " +
                  output + args.side + ".csv")


def multproc(params_, func, n_proc):
    p = multiprocessing.Pool(n_proc)
    p.map(func, params_)
    p.close()


def associate_attributes(params_):
    points, descriptors, outputs = params_

    create_folders('/'.join(outputs.split("/")[:-1]))

    fileP = pd.read_csv(points + ".csv", header=None)
    fileP = fileP.iloc[:, -4:].T.reset_index(drop=True).T

    posic, desc = get_position_descriptor(descriptors)
    land = pd.concat([posic, desc], axis=1).astype(str)
    land.columns = land.columns.astype(str)

    fileP = fileP.astype(int).astype(str)
    fileP.columns = fileP.columns.astype(str)
    final = land.merge(fileP, on=['0', '1', '2'])
    # adding first row
    final.loc[-1] = [np.nan for x in range(0, final.shape[1])]
    final.index = final.index + 1
    final = final.sort_index()
    final.iloc[0, :3] = [final.shape[0]-1, '0', '288']
    final = final.iloc[:, :-1]
    save_csv_desc(final, outputs+"_landmarks.txt")


def generate_final_dataset(input):
    data = []

    for i in range(len(input)):
        _, desc = get_position_descriptor_txt(
            input[i])
        desc = desc.T.reset_index(drop=True).T
        desc['subj'] = input[i].split("/")[-1].\
            split("_landmarks")[0]
        desc['label'] = input[i].split("/")[-2]

        data.append(desc)
    dt = pd.concat(data, sort=True)

    return dt


def match_atlas_structural(filenames):
    input_land = [os.path.join(sw_landmarks, file) for file in filenames]
    outputs = [os.path.join(match_struc_folder, file) for file in filenames]
    params = zip(input_land, outputs, [structural_atlas
                                       for _ in range(len(input_land))])
    multproc(params, match_structural, 35)
    print("Done -- Match structural")


def get_attributes():
    descriptors = [os.path.join(
        desc_landmarks, file + "_" + args.side + "_landmarks.txt")
                        for file in filenames]
    points = [os.path.join(match_struc_folder, file + "_" +
                           args.side) for file in filenames]
    outputs = [os.path.join(desc_text_folder, file + "_" +
                            args.side) for file in filenames]
    params = zip(points, descriptors, outputs)
    multproc(params, associate_attributes, 35)
    print("Done -- get_attributes")


def generate_train_dataset():
    inputs = [os.path.join(desc_text_folder, file + "_" +
                           args.side + "_landmarks.txt") for file in filenames]
    data = generate_final_dataset(inputs)
    save_csv(data, os.path.join(final_desc_folder,
                                "dataset_" + args.side + ".csv"))
    print("Done -- Final dataset")


# ********************************************************************
#
# ********************************************************************

if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    # Define folders
    sw_landmarks = os.path.join('..', 'images', 'hippocampus',
                                'descriptor', 'spiderweb')
    desc_landmarks = os.path.join('..', 'images', 'hippocampus',
                                  'descriptor', 'tissues')

    basename = os.path.join(args.dest_folder, "hippocampus", args.fold)
    structural_atlas = os.path.join(basename, "atlas", args.labels)

    for step in ['train', 'validation', 'test']:
        print("#"*5, step.capitalize(), "#"*5)

        match_struc_folder = os.path.join(basename,
                                          "match_struc_" + step,
                                          args.labels)
        desc_text_folder = os.path.join(basename,
                                        "tissues_descriptors_" + step,
                                        args.labels)

        create_folders(match_struc_folder)
        create_folders(desc_text_folder)

        images = os.path.join(args.dest_folder, "splits",
                              step + "_" + args.fold + ".csv")
        filenames = get_directories_class(images,
                                          args.labels.split("_")[0].upper(),
                                          args.labels.split("_")[1].upper())
        input_land = str(sw_landmarks) + "/" + filenames

        match_atlas_structural(filenames)

        get_attributes()

        if step == 'train':
            final_desc_folder = os.path.join(basename,
                                             "svm_descriptor_" + step,
                                             args.labels)
            create_folders(final_desc_folder)
            generate_train_dataset()

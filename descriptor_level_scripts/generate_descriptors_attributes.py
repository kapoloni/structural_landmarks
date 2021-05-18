#!/usr/bin/python3

# Standard library imports
import getopt
import sys
import os

# Local application imports
from utils import (save_csv,
                   create_folders)
from feature_extraction import FeatureExtraction


def showUsage():
    print('./structural_descriptors.py -n <result_folder>\
          -c <c1> -d <c2> -s <L|R> -f <fold num>')
    sys.exit(2)


def get_params(argv):
    try:
        opts, args = getopt.getopt(argv,
                                   "hn:c:d:s:f:t:",
                                   ["nfile=",
                                    "c1file=", "c2file=",
                                    "sfile=", "fnum="])
    except getopt.GetoptError:
        showUsage()
    for opt, arg in opts:
        if opt == '-h':
            showUsage()
        elif opt in ("-n", "--nfile"):
            result_folder_ = arg
        elif opt in ("-c", "--c1file"):
            c1_ = arg
        elif opt in ("-d", "--c2file"):
            c2_ = arg
        elif opt in ("-s", "--sfile"):
            side_ = arg
        elif opt in ("-f", "--fnum"):
            fold = arg
    return result_folder_, c1_, c2_, side_, fold


def match_atlas_structural(filenames):
    input_land = [os.path.join(match_landmarks, file) for file in filenames]
    outputs = [os.path.join(match_struc_folder, file) for file in filenames]
    atlas_newland_folder = [atlas_new_land_folder
                            for i in range(len(input_land))]
    params = zip(input_land, outputs, atlas_newland_folder)
    fe.multproc(params, fe.match_structural, 35)
    print("Done -- Match structural")


def get_attributes():
    descriptors = [os.path.join(
        desc_landmarks, file + fe._ss + "_landmarks.txt")
                        for file in filenames]
    points = [os.path.join(match_struc_folder, file +
                           fe._ss) for file in filenames]
    outputs = [os.path.join(desc_text_folder, file +
                            fe._ss) for file in filenames]
    params = zip(points, descriptors, outputs)
    fe.multproc(params, fe.associate_attributes, 35)
    print("Done -- get_attributes")


def generate_dataset():
    inputs = [os.path.join(desc_text_folder, file +
                           fe._ss + "_landmarks.txt") for file in filenames]
    data = fe.generate_final_dataset(inputs)
    save_csv(data, os.path.join(final_desc_folder,
                                "dataset" + fe._ss + ".csv"))
    print("Done -- Final dataset")


# ********************************************************************
#
# ********************************************************************

if __name__ == "__main__":

    if len(sys.argv) < 6:
        showUsage()
        exit()

    result_folder, class1, class2, side, fold = get_params(
        sys.argv[1:])
    c1c2 = class1.lower()+"_"+class2.lower()
    fe = FeatureExtraction(c1c2, side)

    # Define folders
    match_landmarks = os.path.join('..', 'images', 'hippocampus',
                                   'descriptor', 'spiderweb')
    desc_landmarks = os.path.join('..', 'images', 'hippocampus',
                                  'descriptor_and_tissues', 'log_sph')

    basename = os.path.join(result_folder, "hippocampus", fold)
    count_matching_folder = os.path.join(basename, "counts", c1c2)
    struc_th_folder = os.path.join(basename, "struc_th", c1c2)
    atlas_new_land_folder = os.path.join(basename, "atlas", c1c2)

    # For split2 train: i.e., training images
    match_struc_folder = os.path.join(basename, "match_struc", c1c2)
    desc_text_folder = os.path.join(basename, "desc", c1c2)
    final_desc_folder = os.path.join(basename, "final_desc", c1c2)
    create_folders(match_struc_folder)
    create_folders(desc_text_folder)
    create_folders(final_desc_folder)

    images = os.path.join(result_folder, "splits", "split_2",
                          "train_" + fold + ".csv")
    filenames = fe.get_directories_class(images, class1, class2)
    input_land = str(match_landmarks) + "/" + filenames

    match_atlas_structural(filenames)

    get_attributes()
    generate_dataset()

    # For split2 test: i.e., validation images
    match_struc_folder = os.path.join(basename, "match_struc_t2", c1c2)
    desc_text_folder = os.path.join(basename, "desc_t2", c1c2)
    create_folders(match_struc_folder)
    create_folders(desc_text_folder)

    images = os.path.join(result_folder, "splits", "split_2",
                          "test_" + fold + ".csv")
    filenames = fe.get_directories_class(images, class1, class2)
    match_atlas_structural(filenames)
    get_attributes()

    # For split split 1 test: i.e., test images
    match_struc_folder = os.path.join(basename, "match_struc_t1", c1c2)
    desc_text_folder = os.path.join(basename, "desc_t1", c1c2)
    create_folders(match_struc_folder)
    create_folders(desc_text_folder)

    images = os.path.join(result_folder, "splits", "split_1",
                          "test_" + fold + ".csv")
    filenames = fe.get_directories_class(images, class1, class2)
    match_atlas_structural(filenames)
    get_attributes()

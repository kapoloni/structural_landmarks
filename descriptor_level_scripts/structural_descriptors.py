#!/usr/bin/python3

# Standard library imports
import getopt
import sys
import os

# Local application imports
from utils import (concat_match_maxd,
                   save_csv,
                   associate_descriptor_atlas,
                   save_csv_desc,
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


def count_matches(filenames):
    in_ = str(match_folder) + "/" + filenames + fe._ss
    c1, c2 = concat_match_maxd(in_, class1, class2)

    maxd_c1 = fe.concat_maxd(c1)
    maxd_c2 = fe.concat_maxd(c2)

    print("Done -- count_points_matches")
    return maxd_c1, maxd_c2


def t_test(c1, c2, alpha):
    struc = fe.distribution(c1, c2, alpha)
    save_csv(struc.iloc[:, :3].astype(int),
             os.path.join(struc_th_folder,
                          "struc" + fe._ss + ".csv"))
    print("Done -- Points label")


def get_structural_atlas_descriptor():
    # Atlas associate_desc
    atlas = fe.atlas
    pts = os.path.join(struc_th_folder, "struc" + fe._ss + ".csv")
    desc = associate_descriptor_atlas(pts, atlas + "_landmarks.txt")
    save_csv_desc(desc, os.path.join(atlas_new_land_folder,
                  "structural" + fe._ss + "_landmarks.txt"))
    print("Done -- get_structural_atlas_descriptor")


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
    match_folder = os.path.join(result_folder, "hippocampus", "match")

    basename = os.path.join(result_folder, "hippocampus", fold)
    # count_matching_folder = os.path.join(basename, "counts", c1c2)
    struc_th_folder = os.path.join(basename, "struc_th", c1c2)
    atlas_new_land_folder = os.path.join(basename, "atlas", c1c2)

    # create_folders(count_matching_folder)
    create_folders(struc_th_folder)
    create_folders(atlas_new_land_folder)

    # Label atlas
    images = os.path.join(result_folder, "splits", "split_1",
                          "train_" + fold + ".csv")

    filenames = fe.get_directories_class(images, class1, class2)

    c1, c2 = count_matches(filenames)

    if c1c2 == 'mci_ad':
        t_test(c1, c2, 0.05)
    else:
        t_test(c1, c2, 0.01)

    get_structural_atlas_descriptor()

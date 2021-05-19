#!/usr/bin/python3
"""
    Match all descritors
"""
# Standard library imports
import os
import argparse
import sys

# Third party imports
import multiprocessing
import pandas as pd

# Local application imports
from utils import (create_folders)


def parse_args(args):
    """!@brief
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description='Match all descriptors')
    parser.add_argument('--dest_folder', help='Destination folder',
                        default='../experiment', type=str)

    return parser.parse_args(args)


def get_directories(csv_file):
    return pd.read_csv(csv_file, header=None)[0]


def match(params_):
    input_land, output, side = params_
    match_in = "../bin/bin/match-landmarks-maxd-prob"
    create_folders('/'.join(output.split("/")[:-1]))
    side_name = "right" if side == '_R' else "left"
    if not os.path.exists(output + side + ".csv"):
        os.system(match_in + " -io " + atlas.replace("side", side_name) + " " +
                  input_land + side + " -to 0.5 > " +
                  output + side + ".csv")


if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    images = "imgs.csv"

    landmarks = os.path.join("..", "images", "hippocampus",
                             "descriptor", "spiderweb")
    match_folder = os.path.join(args.dest_folder, "hippocampus", "match")

    create_folders(match_folder)

    atlas = os.path.join("..", "atlas", "hippocampus", "side", "output",
                         "IXI_20-80", "probabilistic_atlas_averaged")

    filenames = get_directories(images)

    input_land = [os.path.join(landmarks, file) for file in filenames]
    output_ = [os.path.join(match_folder, file) for file in filenames]

    for side in ["L", "R"]:
        sides = ["_" + side for i in range(len(input_land))]
        params_ = zip(input_land, output_, sides)
        p = multiprocessing.Pool(5)
        p.map(match, params_)
        p.close()
    print("Match -- done!")

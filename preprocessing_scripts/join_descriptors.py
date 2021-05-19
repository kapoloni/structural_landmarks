#!/usr/bin/python3
"""
    Join all tissues descriptors
"""
# Standard library imports
import sys
import argparse
import os
import pandas as pd
import numpy as np

# Third party imports
import multiprocessing


def parse_args(args):
    """!@brief
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(
        description='Join all tissues descriptors')
    parser.add_argument('--new_folder', help='New folder',
                        type=str)

    return parser.parse_args(args)


def create_folders(new_folder):
    if not os.path.exists(new_folder):
        os.makedirs(new_folder, exist_ok=True)


def get_directories(root_folder):
    list_in = []
    for class_ in os.listdir(root_folder):
        if os.path.isdir(os.path.join(root_folder, class_)):
            create_folders(os.path.join(args.new_folder, class_))
        for filename in os.listdir(os.path.join(root_folder, class_)):
            if "_L_landmarks.txt" in filename:
                filename = filename.split("_L_landmarks.txt")[0]
                list_in.append(os.path.join(class_, filename))
    return list_in


def get_position_descriptor(filename):
    reader = pd.read_csv(filename, header=None)
    data = reader[0].str.split(" ", expand=True)
    land = data.iloc[1:, :3]
    desc = data.iloc[1:, 3:]
    if np.shape(np.where(desc.iloc[:, -1] == ''))[1] == int(data.iloc[0, 0]):
        desc = desc.iloc[:, :-1]
    return land, desc


def join_descriptors(params_):
    inputfile = params_[0]
    for side in ["L", "R"]:
        print(inputfile)
        output = os.path.join(args.newfolder,
                              inputfile + "_" + side + "_landmarks.txt")
        inputf = DESCGM + inputfile + "_" + side + "_landmarks.txt"
        posic, desc = get_position_descriptor(inputf)
        inputf = pd.concat([posic, desc.iloc[:, 3:]], axis=1)
        inputf.columns = inputf.columns.astype(str)

        for name in [DESCWM, DESCCSF]:
            tmp_inputf = name + inputfile + "_" + side + "_landmarks.txt"
            tmp_posic, tmp_desc = get_position_descriptor(tmp_inputf)
            tmp_inputf = pd.concat([tmp_posic, tmp_desc.iloc[:, 3:]], axis=1)
            tmp_inputf.columns = tmp_inputf.columns.astype(str)

            inputf = inputf.merge(tmp_inputf, on=['0', '1', '2'])
            inputf.columns = [str(i) for i in range(len(inputf.columns))]
        inputf.to_csv(output, index=False, header=False)


if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    DESCGM = "../images/hippocampus/descriptor/gm/"
    DESCWM = "../images/hippocampus/descriptor/wm/"
    DESCCSF = "../images/hippocampus/descriptor/csf/"

    path_in = get_directories(DESCGM)
    p = multiprocessing.Pool(50)
    params = zip(path_in)
    p.map(join_descriptors, params)
    p.close()

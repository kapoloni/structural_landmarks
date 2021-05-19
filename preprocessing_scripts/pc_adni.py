#!/usr/bin/python3
"""
    Phase congruency on ADNI images
"""
# Standard library imports
import sys
import argparse
import os
import subprocess

# Third party imports
import multiprocessing


def parse_args(args):
    """!@brief
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description='Phase Conguency ADNI')
    parser.add_argument('--ref', help='Input images folder',
                        type=str)
    parser.add_argument('--new_folder', help='New folder',
                        type=str)

    return parser.parse_args(args)


def create_folders(new_folder):
    if not os.path.exists(new_folder):
        os.makedirs(new_folder, exist_ok=True)


def get_directories(root_folder):
    list_in = []
    list_out = []
    for class_ in os.listdir(root_folder):
        if os.path.isdir(os.path.join(root_folder, class_)):
            create_folders(os.path.join(args.new_folder, class_))
        for filename in os.listdir(os.path.join(root_folder, class_)):
            if "_mask.nii.gz" not in filename:
                filename = filename.split(".nii.gz")[0]
                list_in.append(os.path.join(
                    root_folder, class_, filename))
                list_out.append(os.path.join(class_, filename))
    return list_in, list_out


def pc(params_):
    inputfile, outputfile = params_
    if not os.path.exists(os.path.join(args.new_folder, outputfile +
                                       "_eigenvalues_0.nii.gz")):
        subprocess.call(["time", PHASE, "-i", inputfile + ".nii.gz",
                         "-m", inputfile + "_mask.nii.gz",
                         "-bof", BCO, "-o",
                         os.path.join(args.new_folder, outputfile)])


if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    BCO = os.path.join("..", "filters",
                       "logGab_#_256_256_256_04_06_03_033_209_055_120_0.bof")
    PHASE = os.path.join("..", "bin", "bin", "phase-congruency")

    path_in, path_out = get_directories(args.ref)
    p = multiprocessing.Pool(10)
    params = zip(path_in, path_out)
    p.map(pc, params)
    p.close()

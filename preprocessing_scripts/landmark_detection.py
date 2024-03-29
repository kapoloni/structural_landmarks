#!/usr/bin/python3
"""
    Describes a point position based on the pc map
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
    parser = argparse.ArgumentParser(description='Landmarks detection')
    parser.add_argument('--ref', help='Input images folder',
                        type=str)
    parser.add_argument('--pc', help='PC images folder',
                        type=str)
    parser.add_argument('--land', help='Landmarks folder',
                        type=str)
    parser.add_argument('--type', help='sw (spider web) or ls(log spherical)',
                        type=str)

    return parser.parse_args(args)


def create_folders(new_folder):
    if not os.path.exists(new_folder):
        os.makedirs(new_folder, exist_ok=True)


def get_directories(root_folder):
    list_in = []
    list_mask = []
    list_out = []
    for class_ in os.listdir(root_folder):
        if os.path.isdir(os.path.join(root_folder, class_)):
            create_folders(os.path.join(args.pc, class_))
            create_folders(os.path.join(args.land, class_))
        for filename in os.listdir(os.path.join(root_folder, class_)):
            if ".nii.gz" in filename and "mask" not in filename:
                filename = filename.split(".nii.gz")[0]

                list_in.append(os.path.join(root_folder, class_, filename))
                list_mask.append(os.path.join(MASK, class_, filename))
                list_out.append(os.path.join(class_, filename))
    return list_in, list_mask, list_out


def pc_land(params_):
    inputfile, maskfile, outputfile = params_
    if not os.path.exists(os.path.join(args.land,
                                       outputfile + "_L_landmarks.txt")):
        subprocess.call(["time", LAND, "-r", inputfile + ".nii.gz", "-io",
                         os.path.join(args.pc, outputfile), desc, "-st",
                         "0.01", "-na", na, '-nr', nr, g, rmax,
                         "-getp", "-mask", maskfile+"_L"+struc+".nii.gz",
                         "-o", os.path.join(args.land, outputfile + "_L")])

    if not os.path.exists(os.path.join(args.land,
                          outputfile + "_R_landmarks.txt")):
        subprocess.call(["time", LAND, "-r", inputfile + ".nii.gz",
                         "-io", os.path.join(args.pc, outputfile),
                         desc, "-st", "0.01", "-na", na, '-nr', nr, g, rmax,
                         "-getp", "-mask", maskfile+"_R"+struc+".nii.gz",
                         "-o", os.path.join(args.land, outputfile + "_R")])

# ********************************************************************
#
# ********************************************************************


if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    LAND = os.path.join("..", "bin", "bin", "bin/landmark-detector")
    MASK = os.path.join("databases", "data1_study1", "AD", "MRI",
                        "Katia", "ADNI", "experiment", "ADNI", "mesh", "r8")

    # Descriptor params spider web
    if args.type == 'sw':
        nr, na, rmax = ['5', '10', '32']
        desc = '-getl'
        g = '-gl'
    else:
        # Descriptor params log spherical
        nr, na, rmax = ['3', '8', '32']
        desc = '-getg'
        g = '-gr'

    struc = "H"
    path_in, path_mask, path_out = get_directories(args.ref)
    p = multiprocessing.Pool(30)
    params = zip(path_in, path_mask, path_out)
    p.map(pc_land, params)
    p.close()

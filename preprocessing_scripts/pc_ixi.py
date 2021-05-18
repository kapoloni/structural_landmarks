#!/usr/bin/python3
import sys
import getopt
import os
import subprocess
from pathlib import Path
import multiprocessing


BCO = "../atlas/filter/logGab_#_256_256_256_04_06_03_033_209_055_120_0.bof"
PHASE = "../bin/bin/phase-congruency"


# ********************************************************************
#
# ********************************************************************


def showUsage():
    print('./pc_ixi.py -r <root_folder> -n <new_folder>')
    sys.exit(2)
# ********************************************************************
#
# ********************************************************************


def create_folders(new_folder):
    if not os.path.exists(new_folder):
        os.makedirs(new_folder, exist_ok=True)

# ********************************************************************
#
# ********************************************************************


def get_directories(root_folder):
    list_in = []
    list_out = []
    create_folders(new_folder)
    for filename in os.listdir(root_folder):
        if "_mask.nii.gz" in filename:
            file_out = filename.split("_mask.nii.gz")[0]
            file_root = root_folder+"/"+filename.split("_mask.nii.gz")[0]
            list_in.append(file_root)
            list_out.append(file_out)
    return list_in, list_out

# ********************************************************************
#
# ********************************************************************


def pc(params_):
    inputfile, outputfile = params_
    if not Path(new_folder + "/" + outputfile +
                "_eigenvalues_2.nii.gz").exists():
        print("time", PHASE, "-i", inputfile + ".nii.gz",
              "-m", inputfile + "_mask.nii.gz", "-bof", BCO,
              "-o", new_folder + "/" + outputfile)
        subprocess.call(["time", PHASE, "-i", inputfile + ".nii.gz",
                         "-m", inputfile + "_mask.nii.gz",
                         "-bof", BCO, "-o", new_folder + "/" + outputfile])
# ********************************************************************
#
# ********************************************************************


if __name__ == "__main__":
    if len(sys.argv) < 5:
        showUsage()
        exit()
    argv = sys.argv[1:]

    root_file = ''

    try:
        opts, args = getopt.getopt(argv, "hr:n:", ["rfile=", "nfile="])
    except getopt.GetoptError:
        showUsage()
    for opt, arg in opts:
        if opt == '-h':
            showUsage()
        elif opt in ("-r", "--rfile"):
            root_folder = arg
        elif opt in ("-n", "--npcfile"):
            new_folder = arg

    path_in, path_out = get_directories(root_folder)
    p = multiprocessing.Pool(5)
    params = zip(path_in, path_out)
    p.map(pc, params)
    p.close()
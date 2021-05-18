#!/usr/bin/python3
import sys
import getopt
import os
import subprocess
from pathlib import Path
import multiprocessing


LAND = "../bin/bin/landmark-detector-tissues"
MASK = "/databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/mesh/r8/"
TISSUE = "/databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/Fast/"


# ********************************************************************
#
# ********************************************************************
def showUsage():
    print('./land_adni.py -r <root_folder> -p <new_folder>' +
          '-l <new_folder> -t <csf | gm | wm>')
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
    list_mask = []
    list_out = []
    for class_ in os.listdir(root_folder):
        if os.path.isdir(root_folder + "/" + class_):
            create_folders(new_folder_pc + "/" + class_)
            create_folders(new_folder_l + "/" + class_)
        for filename in os.listdir(root_folder + "/" + class_):
            if ".nii.gz" in filename:
                if "mask" not in filename:
                    filename = filename.split(".nii.gz")[0]
                    file_out = class_ + "/" + filename
                    file_root = root_folder+"/"+class_+"/" + filename
                    file_mask = class_ + "/" + filename
                    list_in.append(file_root)
                    list_mask.append(file_mask)
                    list_out.append(file_out)
    return list_in, list_mask, list_out


# ********************************************************************
#
# ********************************************************************


def pc_land(params_):
    st = "0.01"
    inputfile, filename, outputfile = params_
    if not Path(new_folder_l + "/" + outputfile + "_L_landmarks.txt").exists():
        # print(TISSUE + filename)
        print("time", LAND, "-r", inputfile + ".nii.gz", "-io",
              new_folder_pc + "/" + outputfile, desc, "-st", st,
              "-na", na, '-nr', nr, g, rmax, "-getp",
              "-m", TISSUE + filename + "_pve_" + pve + ".nii.gz", "-mask",
              MASK + filename + "_L"+struc+".nii.gz", "-o",
              new_folder_l + "/" + outputfile + "_L")
        subprocess.call(["time", LAND, "-r", inputfile + ".nii.gz", "-io",
                         new_folder_pc + "/" + outputfile, desc, "-st",
                         st, "-na", na, '-nr', nr, g, rmax,
                         "-getp", "-m", TISSUE + filename +
                         "_pve_" + pve + ".nii.gz",
                         "-mask", MASK + filename + "_L"+struc+".nii.gz",
                         "-o", new_folder_l + "/" + outputfile + "_L"])
    if not Path(new_folder_l + "/" + outputfile + "_R_landmarks.txt").exists():
        # print(TISSUE + filename)
        print("time", LAND, "-r", inputfile + ".nii.gz", "-io",
              new_folder_pc + "/" + outputfile, desc, "-st", st,
              "-na", na, '-nr', nr, g, rmax, "-getp",
              "-m", TISSUE + filename + "_pve_" + pve + ".nii.gz", "-mask",
              MASK + filename + "_R"+struc+".nii.gz", "-o",
              new_folder_l + "/" + outputfile + "_R")
        subprocess.call(["time", LAND, "-r", inputfile + ".nii.gz",
                         "-io", new_folder_pc + "/" + outputfile,
                         desc, "-st", st, "-na", na, '-nr', nr, g, rmax,
                         "-getp",  "-m", TISSUE + filename + "_pve_" +
                         pve + ".nii.gz",
                         "-mask", MASK + filename + "_R"+struc+".nii.gz",
                         "-o", new_folder_l + "/" + outputfile + "_R"])

# ********************************************************************
#
# ********************************************************************


if __name__ == "__main__":
    if len(sys.argv) < 7:
        showUsage()
        exit()
    argv = sys.argv[1:]

    root_file = ''

    try:
        opts, args = getopt.getopt(argv, "hr:p:l:t:",
                                         ["rfile=",
                                          "npcfile=",
                                          "nlfile=",
                                          "type="])
    except getopt.GetoptError:
        showUsage()
    # print(opts)
    for opt, arg in opts:
        if opt == '-h':
            showUsage()
        elif opt in ("-r", "--rfile"):
            root_folder = arg
        elif opt in ("-p", "--npcfile"):
            new_folder_pc = arg
        elif opt in ("-l", "--nlfile"):
            new_folder_l = arg
        elif opt in ("-t", "--type"):
            type = arg

    if type == 'cfs':
        pve = "0"  # csf
    elif type == 'gm':
        pve = "1"  # gm
    else:
        pve = "2"  # wm

    # Log Spherical
    nr, na, rmax = ['3', '8', '32']
    desc = '-getg'
    g = '-gr'
    struc = "H"

    path_in, path_mask, path_out = get_directories(root_folder)
    # print(path_in, path_out)
    p = multiprocessing.Pool(30)
    params = zip(path_in, path_mask, path_out)
    p.map(pc_land, params)
    p.close()

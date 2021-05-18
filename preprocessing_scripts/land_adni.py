#!/usr/bin/python3
import sys
import getopt
import os
import subprocess
from pathlib import Path
import multiprocessing


LAND = "../bin/bin/landmark-detector"
MASK = "/databases/data1_study1/AD/MRI/Katia/ADNI2/experiment/ADNI/mesh/r8/"
# ********************************************************************
#
# ********************************************************************


def showUsage():
    print('./land_adni.py -r <root_folder> -p <new_folder> -l ' +
          '<new_folder> -t <spider-web (s) | log_spherical (l)>')
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
                    file_mask = MASK + "/" + class_ + "/" + filename
                    list_in.append(file_root)
                    list_mask.append(file_mask)
                    list_out.append(file_out)
    return list_in, list_mask, list_out


# ********************************************************************
#
# ********************************************************************

def pc_land(params_):
    inputfile, maskfile, outputfile = params_
    if not Path(new_folder_l + "/" + outputfile + "_R_landmarks.txt").exists():
        print("time", LAND, "-r", inputfile + ".nii.gz", "-io",
              new_folder_pc + "/" + outputfile, desc, "-st", "0.01",
              "-na", na, '-nr', nr, g, rmax, "-getp", "-mask",
              maskfile+"_L"+struc+".nii.gz", "-o",
              new_folder_l + "/" + outputfile + "_L")
        subprocess.call(["time", LAND, "-r", inputfile + ".nii.gz", "-io",
                         new_folder_pc + "/" + outputfile, desc, "-st",
                         "0.01", "-na", na, '-nr', nr, g, rmax,
                         "-getp", "-mask", maskfile+"_L"+struc+".nii.gz",
                         "-o", new_folder_l + "/" + outputfile + "_L"])
        print("time", LAND, "-r", inputfile + ".nii.gz", "-io",
              new_folder_pc + "/" + outputfile, desc, "-st", "0.01",
              "-na", na, '-nr', nr, g, rmax, "-getp", "-mask",
              maskfile+"_R"+struc+".nii.gz", "-o",
              new_folder_l + "/" + outputfile + "_R")
        subprocess.call(["time", LAND, "-r", inputfile + ".nii.gz",
                         "-io", new_folder_pc + "/" + outputfile,
                         desc, "-st", "0.01", "-na", na, '-nr', nr, g, rmax,
                         "-getp", "-mask", maskfile+"_R"+struc+".nii.gz",
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
                                          "type"])
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

    # Descriptor params spider web
    if type == 's':
        nr, na, rmax = ['5', '10', '32']
        desc = '-getl'
        g = '-gl'
    else:
        # Descriptor params log spherical
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

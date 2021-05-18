#!/usr/bin/python3
import sys
import getopt
import os
import pandas as pd
import numpy as np
import multiprocessing

# ********************************************************************
#
# ********************************************************************


def showUsage():
    print('./join_descriptors.py -n <new_folder>')
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
    for class_ in os.listdir(root_folder):
        if os.path.isdir(root_folder + "/" + class_):
            create_folders(NEWFOLDER + "/" + class_)
        for filename in os.listdir(root_folder + "/" + class_):
            if "_L_landmarks.txt" in filename:
                filename = filename.split("_L_landmarks.txt")[0]
                input_img = class_ + "/" + filename
                list_in.append(input_img)
    return list_in


# ********************************************************************
#
# ********************************************************************

def get_position_descriptor(filename):
    reader = pd.read_csv(filename, header=None)
    data = reader[0].str.split(" ", expand=True)
    land = data.iloc[1:, :3]
    desc = data.iloc[1:, 3:]
    if np.shape(np.where(desc.iloc[:, -1] == ''))[1] == int(data.iloc[0, 0]):
        desc = desc.iloc[:, :-1]
    return land, desc

# ********************************************************************
#
# ********************************************************************


def join_descriptors(params_):
    inputfile = params_[0]
    for side in ["L", "R"]:
        print(inputfile)
        output = NEWFOLDER + "/" + inputfile + "_" + side + "_landmarks.txt"
        inputf = DESCSHAPE + inputfile + "_" + side + "_landmarks.txt"
        posic, desc = get_position_descriptor(inputf)
        inputf = pd.concat([posic, desc.iloc[:, 3:]], axis=1)
        inputf.columns = inputf.columns.astype(str)

        for name in [DESCGM, DESCWM, DESCCSF]:
            tmp_inputf = name + inputfile + "_" + side + "_landmarks.txt"
            tmp_posic, tmp_desc = get_position_descriptor(tmp_inputf)
            tmp_inputf = pd.concat([tmp_posic, tmp_desc.iloc[:, 3:]], axis=1)
            tmp_inputf.columns = tmp_inputf.columns.astype(str)

            inputf = inputf.merge(tmp_inputf, on=['0', '1', '2'])
            inputf.columns = [str(i) for i in range(len(inputf.columns))]
        inputf.to_csv(output, index=False, header=False)

# ********************************************************************
#
# ********************************************************************


if __name__ == "__main__":

    STRUCT = "hippo"

    DESCGM = "../../images/" + STRUCT + "/landgm/32/"
    DESCWM = "../../images/" + STRUCT + "/landwm/32/"
    DESCCSF = "../../images/" + STRUCT + "/landcsf/32/"
    DESCSHAPE = "../../images/" + STRUCT + "/landshape/32/"

    if len(sys.argv) < 2:
        showUsage()
        exit()
    argv = sys.argv[1:]

    root_file = ''

    try:
        opts, args = getopt.getopt(argv, "hn:", ["nfile="])
    except getopt.GetoptError:
        showUsage()
    # print(opts)
    for opt, arg in opts:
        if opt == '-h':
            showUsage()
        elif opt in ("-n", "--nfile"):
            NEWFOLDER = arg

    path_in = get_directories(DESCSHAPE)
    p = multiprocessing.Pool(50)
    params = zip(path_in)
    p.map(join_descriptors, params)
    p.close()

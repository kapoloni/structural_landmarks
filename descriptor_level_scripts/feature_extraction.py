"""
Feature extraction functions and paths
"""

# Standard library imports
from __future__ import absolute_import
from __future__ import print_function
import os
import numpy as np
import pandas as pd

# Third party imports
from pathlib import Path
import multiprocessing

# Local application imports
from utils import (create_folders,
                   get_position_descriptor,
                   save_csv_desc)


class FeatureExtraction():

    def __init__(self, c1c2, side):

        atlas = "../atlas/hippocampus/" + \
            "side/output/IXI_20-80/probabilistic_atlas_averaged"

        self.atlas = atlas.replace("side", "left" if side == "L" else "right")

        self.side = side

        self.labels = c1c2

    def match_structural(self, params_):
        input_land, output, atlas_new_land_folder = params_
        match_in = "../bin/bin/match-landmarks-maxd-prob"

        atlas_ = os.path.join(atlas_new_land_folder, "structural" + self._ss)
        create_folders('/'.join(output.split("/")[:-1]))
        if not Path(output + "_" + self.side + ".csv").exists():
            os.system(match_in + " -io " + atlas_ + " " + input_land +
                      "_" + self.side +" -to 0.5 > " + output + self._ss + ".csv")

    def multip(self, input_land, path_out_, func, n_proc):
        params_ = zip(input_land, path_out_)
        p = multiprocessing.Pool(n_proc)
        p.map(func, params_)
        p.close()

    def multproc(self, params_, func, n_proc):
        p = multiprocessing.Pool(n_proc)
        p.map(func, params_)
        p.close()

    def associate_attributes(self, params_):
        points, descriptors, outputs = params_

        create_folders('/'.join(outputs.split("/")[:-1]))

        fileP = pd.read_csv(points + ".csv", header=None)
        fileP = fileP.iloc[:, -4:].T.reset_index(drop=True).T

        posic, desc = get_position_descriptor(descriptors)
        land = pd.concat([posic, desc], axis=1).astype(str)
        land.columns = land.columns.astype(str)

        fileP = fileP.astype(int).astype(str)
        fileP.columns = fileP.columns.astype(str)
        final = land.merge(fileP, on=['0', '1', '2'])
        # adding first row
        final.loc[-1] = [np.nan for x in range(0, final.shape[1])]
        final.index = final.index + 1
        final = final.sort_index()
        final.iloc[0, :3] = [final.shape[0]-1, '0', '288']
        final = final.iloc[:, :-1]
        save_csv_desc(final, outputs+"_landmarks.txt")

    # def associate_descriptor_database(self, params_):
    #     points_, descriptors, outputs = params_
    #     points = points_ + ".csv"
    #     create_folders('/'.join(outputs.split("/")[:-1]))
    #     fileP = pd.read_csv(points, header=None)
    #     fileP = fileP.iloc[:, -4:].T.reset_index(drop=True).T
    #     posic, desc = get_position_descriptor(descriptors,
    #                                           first=True)
    #     land = pd.concat([posic, desc], axis=1)
    #     land = land.astype(str)
    #     land.columns = land.columns.astype(str)
    #     fileP = fileP.astype(int).astype(str)
    #     fileP.columns = fileP.columns.astype(str)
    #     final = land.merge(fileP, on=['0', '1', '2'])
    #     # adding first row
    #     final.loc[-1] = [np.nan for x in range(0, final.shape[1])]
    #     final.index = final.index + 1
    #     final = final.sort_index()
    #     final.iloc[0, :3] = [final.shape[0]-1, '0', '128']
    #     save_csv_desc(final, outputs+"_landmarks.txt")

    # def calculate_texture_descritor_glcm(self, params_):
    #     images, descriptors_, desc2_, outputs, outputs2 = params_
    #     glcm = "../../bin/bin/glcm_landmarks"
    #     descriptors = descriptors_ + "_landmarks.txt"
    #     desc2 = desc2_ + "_landmarks.txt"
    #     haralic_folder = outputs+"_text.txt"
    #     create_folders('/'.join(outputs.split("/")[:-1]))
    #     create_folders('/'.join(outputs2.split("/")[:-1]))
    #     if not Path(haralic_folder).exists():
    #         print(glcm + " -i " + images + " -p " +
    #               descriptors + " -ps " + str(self.patch_size) +
    #               " -nr 1 -s " + str(self.roi_size) +
    #               " > " + haralic_folder)

    #     posic_sc, desc_sc = get_position_descriptor(desc2)

    #     shape = pd.concat([posic_sc, desc_sc], axis=1)
    #     posic_hc, desc_hc = get_position_descriptor(haralic_folder)

    #     haralic = pd.concat([posic_hc, desc_hc], axis=1)
    #     final = shape.merge(haralic, on=[0, 1, 2])
    #     # adding first row
    #     final.loc[-1] = [np.nan for x in range(0, final.shape[1])]
    #     final.index = final.index + 1
    #     final = final.sort_index()
    #     final.iloc[0, :3] = [final.shape[0]-1, '0', '136']

    #     save_csv_desc(final, outputs2+"_text_desc.txt")

    # def calculate_texture_descritor_rlm(self, params_):
        # images, descriptors_, desc2_, outputs, outputs2 = params_
        # rlm = "../../bin/bin/rlm_landmarks"
        # descriptors = descriptors_ + "_landmarks.txt"
        # desc2 = desc2_ + "_text_desc.txt"
        # rlm_folder = outputs+"_text.txt"
        # create_folders('/'.join(outputs.split("/")[:-1]))
        # create_folders('/'.join(outputs2.split("/")[:-1]))
        # if not Path(rlm_folder).exists():
        #     os.system(rlm + " -i " + images + " -p " +
        #               descriptors + " -ps " + str(self.patch_size) +
        #               " -nr 1 -s " + str(self.roi_size) +
        #               " > " + rlm_folder)
        # posic_sc, desc_sc = get_position_descriptor(desc2)
        # shape = pd.concat([posic_sc, desc_sc], axis=1)
        # posic_hc, desc_hc = get_position_descriptor(rlm_folder)
        # haralic = pd.concat([posic_hc, desc_hc], axis=1)
        # final = shape.merge(haralic, on=[0, 1, 2])
        # # adding first row
        # final.loc[-1] = [np.nan for x in range(0, final.shape[1])]
        # final.index = final.index + 1
        # final = final.sort_index()
        # final.iloc[0, :3] = [final.shape[0]-1, '0', '146']

        # save_csv_desc(final, outputs2+"_text_desc.txt")

    # def get_directories_class(self, csv_file, c1, c2):
    #     reader = pd.read_csv(csv_file, sep=",", header=None)
    #     reader.columns = ['subj', 'label']
    #     c1 = reader[reader.label == c1]
    #     c2 = reader[reader.label == c2]
    #     return pd.concat([c1.subj, c2.subj]).reset_index(drop=True)

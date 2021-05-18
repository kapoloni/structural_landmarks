#!/usr/bin/python3
# Standard library imports
import os
import numpy as np
import pandas as pd

# Third party imports
from scipy import stats
from pathlib import Path
import multiprocessing

# Local application imports
from utils import (create_folders,
                   get_position_descriptor,
                   save_csv_desc,
                   get_position_descriptor_txt)


class FeatureExtraction():

    def __init__(self, c1c2, side):

        atlas_, template_, mask_, \
            self.roi_size = self.get_atributes()

        self.atlas, self.template, \
            self.mask = self.define_variables_names(
                        atlas_, template_, mask_, side)

        self.side = side

        self.class_ = c1c2

        self.patch_size = 32

        self.ss, self._ss = self.get_sides_strings(side)

        self.th0, self.th1, self.th_ = 0.22, 0.4, 0.1

    def define_variables_names(self, atlas_, template_, mask_, side):
        side_name = "left" if side == "L" else "right"
        template_ = template_.replace("Side", side_name.capitalize())
        atlas_ = atlas_.replace("side", side_name)
        mask_ = mask_.replace("Side", side_name.capitalize())
        return atlas_, template_, mask_

    def get_atributes(self):
        template_ = "../template/NAC_T1_RAI.nii.gz"

        atlas_ = "../atlas/hippocampus/" + \
                 "side/output/IXI_20-80/probabilistic_atlas_averaged"
        mask_ = "../template/hippocampus/Side_256.nii.gz"
        roi_size_ = 64

        return atlas_, template_, mask_, roi_size_

    def get_sides_strings(self, side):
        if side == "L":
            ss = "_L_"
        elif side == "R":
            ss = "_R_"
        elif side == "":
            ss = ""

        if side == "L" or side == "R":
            _ss = "_"+ss.replace("_", "")
        else:
            _ss = ""
        return str(ss), str(_ss)

    def match_structural(self, params_):
        input_land, output, atlas_new_land_folder = params_
        match_in = "../bin/bin/match-landmarks-maxd-prob"

        atlas_ = os.path.join(atlas_new_land_folder, "structural" + self._ss)
        create_folders('/'.join(output.split("/")[:-1]))
        if not Path(output + self._ss + ".csv").exists():
            os.system(match_in + " -io " + atlas_ + " " + input_land +
                      self._ss+" -to 0.5 > " + output + self._ss + ".csv")

    def concat_maxd(self, dt):
        point_posic, counts = np.unique(dt[[0, 1, 2]],
                                        return_counts=True, axis=0)
        max_ = np.max(counts)+3
        dist = np.full((len(point_posic), max_+3), -1, dtype=float)
        for i, point in enumerate(point_posic):
            v = dt[[0, 1, 2]] == point
            v = dt.iloc[np.where(v.all(axis=1))[0]]
            maxd = v.iloc[:, -1].transpose()
            dist[i, :3] = point
            dist[i, 3:len(maxd)+3] = maxd

        return pd.DataFrame(dist)

    def t_test(self, row, alpha):
        c1 = row[row.index.str.contains('_c1')]
        c2 = row[row.index.str.contains('_c2')]

        c1 = c1[c1 > 0]
        c2 = c2[c2 > 0]
        t_score = stats.ttest_ind_from_stats(mean1=c2.mean(),
                                             std1=c2.std(ddof=1),
                                             nobs1=len(c2),
                                             mean2=c1.mean(),
                                             std2=c1.std(ddof=1),
                                             nobs2=len(c1),
                                             equal_var=True)
        stat, p_value = t_score
        p_value /= 2
        if p_value < alpha and stat > 0:
            return stat
        else:
            return np.nan

    def distribution(self, c1, c2, alpha):
        c1 = c1.reset_index(drop=True)
        c2 = c2.reset_index(drop=True)
        c1.columns = [str(col) + '_c1' if col not in [0, 1, 2] else col 
                      for col in c1.columns]
        c2.columns = [str(col) + '_c2' if col not in [0, 1, 2] else col
                      for col in c2.columns]

        dt_ = c1.merge(c2, on=[0, 1, 2])
        dt = dt_.iloc[:, 3:]
        print("before", len(dt))
        dt['keep'] = dt.apply(self.t_test, axis=1, alpha=alpha)
        pts = pd.concat([dt_.iloc[:, :3], dt['keep']], axis=1)
        pts = pts.replace([np.inf, -np.inf], np.nan)
        pts = pts.dropna().reset_index(drop=True)
        print("after", len(pts))
        return pts

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
        final.iloc[0, :3] = [final.shape[0]-1, '0', '146']
        final = final.iloc[:, :-1]
        save_csv_desc(final, outputs+"_landmarks.txt")

    def associate_descriptor_database(self, params_):
        points_, descriptors, outputs = params_
        points = points_ + ".csv"
        create_folders('/'.join(outputs.split("/")[:-1]))
        fileP = pd.read_csv(points, header=None)
        fileP = fileP.iloc[:, -4:].T.reset_index(drop=True).T
        posic, desc = get_position_descriptor(descriptors,
                                              first=True)
        land = pd.concat([posic, desc], axis=1)
        land = land.astype(str)
        land.columns = land.columns.astype(str)
        fileP = fileP.astype(int).astype(str)
        fileP.columns = fileP.columns.astype(str)
        final = land.merge(fileP, on=['0', '1', '2'])
        # adding first row
        final.loc[-1] = [np.nan for x in range(0, final.shape[1])]
        final.index = final.index + 1
        final = final.sort_index()
        final.iloc[0, :3] = [final.shape[0]-1, '0', '128']
        save_csv_desc(final, outputs+"_landmarks.txt")

    def calculate_texture_descritor_glcm(self, params_):
        images, descriptors_, desc2_, outputs, outputs2 = params_
        glcm = "../../bin/bin/glcm_landmarks"
        descriptors = descriptors_ + "_landmarks.txt"
        desc2 = desc2_ + "_landmarks.txt"
        haralic_folder = outputs+"_text.txt"
        create_folders('/'.join(outputs.split("/")[:-1]))
        create_folders('/'.join(outputs2.split("/")[:-1]))
        if not Path(haralic_folder).exists():
            print(glcm + " -i " + images + " -p " +
                  descriptors + " -ps " + str(self.patch_size) +
                  " -nr 1 -s " + str(self.roi_size) +
                  " > " + haralic_folder)

        posic_sc, desc_sc = get_position_descriptor(desc2)

        shape = pd.concat([posic_sc, desc_sc], axis=1)
        posic_hc, desc_hc = get_position_descriptor(haralic_folder)

        haralic = pd.concat([posic_hc, desc_hc], axis=1)
        final = shape.merge(haralic, on=[0, 1, 2])
        # adding first row
        final.loc[-1] = [np.nan for x in range(0, final.shape[1])]
        final.index = final.index + 1
        final = final.sort_index()
        final.iloc[0, :3] = [final.shape[0]-1, '0', '136']

        save_csv_desc(final, outputs2+"_text_desc.txt")

    def calculate_texture_descritor_rlm(self, params_):
        images, descriptors_, desc2_, outputs, outputs2 = params_
        rlm = "../../bin/bin/rlm_landmarks"
        descriptors = descriptors_ + "_landmarks.txt"
        desc2 = desc2_ + "_text_desc.txt"
        rlm_folder = outputs+"_text.txt"
        create_folders('/'.join(outputs.split("/")[:-1]))
        create_folders('/'.join(outputs2.split("/")[:-1]))
        if not Path(rlm_folder).exists():
            os.system(rlm + " -i " + images + " -p " +
                      descriptors + " -ps " + str(self.patch_size) +
                      " -nr 1 -s " + str(self.roi_size) +
                      " > " + rlm_folder)
        posic_sc, desc_sc = get_position_descriptor(desc2)
        shape = pd.concat([posic_sc, desc_sc], axis=1)
        posic_hc, desc_hc = get_position_descriptor(rlm_folder)
        haralic = pd.concat([posic_hc, desc_hc], axis=1)
        final = shape.merge(haralic, on=[0, 1, 2])
        # adding first row
        final.loc[-1] = [np.nan for x in range(0, final.shape[1])]
        final.index = final.index + 1
        final = final.sort_index()
        final.iloc[0, :3] = [final.shape[0]-1, '0', '146']

        save_csv_desc(final, outputs2+"_text_desc.txt")

    def get_directories_class(self, csv_file, c1, c2):
        reader = pd.read_csv(csv_file, sep=",", header=None)
        reader.columns = ['subj', 'label']
        c1 = reader[reader.label == c1]
        c2 = reader[reader.label == c2]
        return pd.concat([c1.subj, c2.subj]).reset_index(drop=True)

    def generate_final_dataset(self, input):
        data = []

        for i in range(len(input)):
            _, desc = get_position_descriptor_txt(
                input[i])
            desc = desc.T.reset_index(drop=True).T
            desc['subj'] = input[i].split("/")[-1].\
                split("_landmarks")[0]
            desc['label'] = input[i].split("/")[-2]

            data.append(desc)
        dt = pd.concat(data, sort=True)

        return dt

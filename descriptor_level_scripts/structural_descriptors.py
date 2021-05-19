#!/usr/bin/python3
"""
    Generate structural points
"""
# Standard library imports
from __future__ import absolute_import
from __future__ import print_function
import argparse
import numpy as np
import pandas as pd
import sys
import os

# Third party imports
from scipy import stats

# Local application imports
from utils import (save_csv,
                   get_directories_class,
                   get_position_descriptor_txt,
                   save_csv_desc,
                   create_folders)


def parse_args(args):
    """!@brief
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description='Generate descriptors')
    parser.add_argument('--dest_folder', help='Destination folder',
                        default='../experiment', type=str)
    parser.add_argument('--labels', help='Labels to perform the experiment',
                        default='cn_ad', type=str)
    parser.add_argument('--side', help='Hippocampal side',
                        default='L', type=str)
    parser.add_argument('--fold', help='Fold number',
                        type=str)
    return parser.parse_args(args)


def get_maxd(input):
    reader = pd.read_csv(input+".csv", header=None)
    out = pd.concat([reader.iloc[:, :3], reader.iloc[:, -1]], axis=1)
    return out


def concat_match_maxd(input, c1, c2):
    c1_ = pd.DataFrame()
    c2_ = pd.DataFrame()
    for v in range(input.shape[0]):
        if type(input.iloc[v]) is str:
            if not os.stat(input.iloc[v]+".csv").st_size == 0:
                if "/"+c1+"/" in input.iloc[v]:
                    c1_ = pd.concat([c1_, get_maxd(input.iloc[v])])
                if "/"+c2+"/" in input.iloc[v]:
                    c2_ = pd.concat([c2_, get_maxd(input.iloc[v])])

    return c1_, c2_


def concat_maxd(dt):
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


def count_matches(filenames):
    in_ = str(match_folder) + "/" + filenames + "_" + args.side
    c1, c2 = concat_match_maxd(in_, args.labels.split("_")[0].upper(),
                               args.labels.split("_")[1].upper())
    maxd_c1 = concat_maxd(c1)
    maxd_c2 = concat_maxd(c2)
    print("Done -- count_points_matches")
    return maxd_c1, maxd_c2


def t_test(row, alpha):
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


def distribution(c1, c2, alpha):
    c1 = c1.reset_index(drop=True)
    c2 = c2.reset_index(drop=True)
    c1.columns = [str(col) + '_c1' if col not in [0, 1, 2] else col
                  for col in c1.columns]
    c2.columns = [str(col) + '_c2' if col not in [0, 1, 2] else col
                  for col in c2.columns]

    dt_ = c1.merge(c2, on=[0, 1, 2])
    dt = dt_.iloc[:, 3:]
    print("before", len(dt))
    dt['keep'] = dt.apply(t_test, axis=1, alpha=alpha)
    pts = pd.concat([dt_.iloc[:, :3], dt['keep']], axis=1)
    pts = pts.replace([np.inf, -np.inf], np.nan)
    pts = pts.dropna().reset_index(drop=True)
    print("after", len(pts))
    return pts


def get_structural_points(c1, c2, alpha):
    struc = distribution(c1, c2, alpha)
    save_csv(struc.iloc[:, :3].astype(int),
             os.path.join(struc_th_folder,
                          "struc_" + args.side + ".csv"))
    print("Done -- Points label")


def associate_descriptor_atlas(fileP, fileD):
    fileP = pd.read_csv(fileP, sep="\t", header=None)
    posic, desc = get_position_descriptor_txt(fileD)
    land = pd.concat([posic, desc], axis=1)
    land.iloc[:, :3] = land.iloc[:, :3].astype(np.int64)
    land.columns = land.columns.astype(str)
    land = land.astype(str)
    fileP.columns = fileP.columns.astype(str)
    fileP = fileP.astype(str)
    final = land.merge(fileP, on=['0', '1', '2'])
    # adding first row
    final.loc[-1] = [np.nan for x in range(0, final.shape[1])]
    final.index = final.index + 1
    final = final.sort_index()
    final.iloc[0, :3] = [final.shape[0]-1, '250', '0']
    # final.iloc[0, :3] = [final.shape[0]-1, '0', '200']
    return final


def get_structural_atlas_descriptor():
    side_name = "right" if args.side == 'R' else "left"
    atlas = os.path.join("..", "atlas", "hippocampus", side_name,
                         "output", "IXI_20-80",
                         "probabilistic_atlas_averaged")
    pts = os.path.join(struc_th_folder, "struc_" + args.side + ".csv")
    desc = associate_descriptor_atlas(pts, atlas + "_landmarks.txt")
    save_csv_desc(desc, os.path.join(atlas_new_land_folder,
                  "structural_" + args.side + "_landmarks.txt"))
    print("Done -- get_structural_atlas_descriptor")


if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    # Define folders
    match_folder = os.path.join(args.dest_folder, "hippocampus", "match")

    basename = os.path.join(args.dest_folder, "hippocampus", args.fold)
    struc_th_folder = os.path.join(
        basename, "structural_positions", args.labels)
    atlas_new_land_folder = os.path.join(basename, "atlas", args.labels)

    create_folders(struc_th_folder)
    create_folders(atlas_new_land_folder)

    images = os.path.join(args.dest_folder, "splits",
                          "train_" + args.fold + ".csv")

    filenames = get_directories_class(images,
                                      args.labels.split("_")[0].upper(),
                                      args.labels.split("_")[1].upper())

    c1, c2 = count_matches(filenames)

    if args.labels == 'mci_ad':
        get_structural_points(c1, c2, 0.05)
    else:
        get_structural_points(c1, c2, 0.01)

    get_structural_atlas_descriptor()

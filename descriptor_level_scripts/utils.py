"""
Utility functions
"""

# Standard library imports
from __future__ import absolute_import
from __future__ import print_function
import os

# Third party imports
import numpy as np
import pandas as pd


def save_csv(data, name):
    data.to_csv(name, sep='\t', encoding='utf-8', index=False, header=False)


def save_csv_desc(data, name):
    data.to_csv(name, sep=' ', encoding='utf-8', index=False, header=False)


def get_directories_class(csv_file, c1, c2):
    reader = pd.read_csv(csv_file, sep=",", header=None)
    reader.columns = ['subj', 'label']
    c1 = reader[reader.label == c1]
    c2 = reader[reader.label == c2]
    return pd.concat([c1.subj, c2.subj]).reset_index(drop=True)


def get_position_descriptor_txt(filename, first=False):

    reader = pd.read_csv(filename, header=None)
    data = reader[0].str.split(" ", expand=True)

    land = data.iloc[1:, :3]
    if first:
        desc = data.iloc[1:, 6:]
    else:
        desc = data.iloc[1:, 3:]

    # se a ultima coluna for vazia
    if np.shape(np.where(desc.iloc[:, -1] == ''))[1] == int(data.iloc[0, 0]):
        desc = desc.iloc[:, :-1]

    return land, desc


def get_position_descriptor(filename, first=False):
    # print("read", filename)
    data = pd.read_csv(filename, header=None)

    land = data.iloc[1:, :3]
    if first:
        desc = data.iloc[1:, 6:]
    else:
        desc = data.iloc[1:, 3:]

    return land, desc


def create_folders(name):
    if not os.path.exists(name):
        os.makedirs(name, exist_ok=True)

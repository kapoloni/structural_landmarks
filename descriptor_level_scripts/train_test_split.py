#!/usr/bin/python3

# Standard library imports
import os
import sys
import argparse

# Third party imports
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold


def parse_args(args):
    """!@brief
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description='Split database')
    parser.add_argument('--dest_folder', help='Destination folder',
                        default='../experiment', type=str)

    return parser.parse_args(args)


def split_and_save_train_test(dt, result_folder):
    dt = get_groups(dt)

    skf = StratifiedKFold(n_splits=10)
    skf2 = StratifiedKFold(n_splits=2)
    for i, (train_idx, test_idx) in enumerate(skf.split(dt.subject,
                                                        dt.age_class)):
        train = dt.iloc[train_idx].reset_index(drop=True)
        test = dt.iloc[test_idx].reset_index(drop=True)

        print("Fold", str(i))
        print("Split 1")
        print("Train", train.subject.shape, np.unique(train.label,
                                                      return_counts=True))
        print("Test", test.subject.shape, np.unique(test.label,
                                                    return_counts=True))

        test1 = pd.concat([test.subject, test.label], axis=1)
        test1.to_csv(os.path.join(result_folder, "test_" + str(i) +
                     ".csv"), index=False, header=False)

        # Split2
        dt2 = get_groups(train)

        idx = [(train_idx, test_idx) for train_idx, test_idx in
               skf2.split(dt2.subject, dt2.age_class)][0]

        train = dt2.iloc[idx[0]]
        test = dt2.iloc[idx[1]]

        print('Split2')
        print("Train", train.subject.shape, np.unique(train.label,
                                                      return_counts=True))
        print("Test", test.subject.shape, np.unique(test.label,
                                                    return_counts=True))

        train2 = pd.concat([train.subject, train.label], axis=1)
        test2 = pd.concat([test.subject, test.label], axis=1)

        train2.to_csv(os.path.join(result_folder, "train_" +
                      str(i) + ".csv"), index=False, header=False)
        test2.to_csv(os.path.join(result_folder, "validation_" +
                     str(i) + ".csv"), index=False, header=False)


def get_directories(root_folder, split_):
    list_in = []
    for class_ in os.listdir(root_folder):
        if os.path.isdir(os.path.join(root_folder, class_)):
            for filename in os.listdir(os.path.join(root_folder, class_)):
                if split_ in filename:
                    file_root = os.path.join(root_folder, class_,
                                             filename.split(split_)[0])
                    list_in.append(file_root)
    return list_in


def organize_data_class(data_):
    dt = pd.DataFrame(data_)
    dt = dt[0].str.split("/", expand=True)
    dt = dt[[5, 6]]
    dt.columns = ['label', 'subject']
    dt['subject'] = dt.apply(lambda x: '/'.join(x)[:-1], axis=1).astype(str)

    adni = pd.read_csv("imgs_data.csv")
    adni = adni[['Subject ID', 'Age', 'Sex']]
    adni.columns = ['subj', 'age', 'sex']

    dt['subj'] = dt['subject'].apply(
        lambda x: '_'.join(x.split('/')[-1].split("_")[:3]))
    dt = adni.merge(dt, left_on='subj', right_on='subj')
    dt = dt.reset_index(drop=True)
    return dt


def get_groups(df):
    delta = df.age.max() - df.age.min()
    discrete = pd.cut(df.age,
                      bins=int(delta/4),
                      right=False
                      ).reset_index(drop=True).astype(str)
    group = [(un, i) for i, un in enumerate(discrete.unique())]
    df['age_class'] = discrete.apply(lambda x: str(dict(group)[x]))
    df.age_class = df[['age_class', 'label']].apply(
        lambda x: x[1] + "_" + x[0], axis=1)

    return df


def create_folders(name):
    if not os.path.exists(name):
        os.makedirs(name, exist_ok=True)


if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    result_folder = os.path.join(args.dest_folder, "splits")

    create_folders(result_folder)

    database = os.path.join("..", "images", "hippocampus",
                            "descriptor", "spiderweb")

    input_ = get_directories(database, "L_landmarks.txt")

    dt = organize_data_class(input_)

    split_and_save_train_test(dt,
                              result_folder)

#!/usr/bin/python3
# Standard library imports
import numpy as np
import pandas as pd
import os

"""
    Código para organizar os atributos de treino e teste extraídos
    das predições dos descritores e gerar um arquivo de entrada
    para o classificador no nível da imagem.
"""


def fill_data(dt, idx):
    dt_frac = dt.iloc[idx]
    dt0 = dt_frac[['0', 'p0']]
    dt1 = dt_frac[['1', 'p1']]
    for col0, col1 in zip(dt0.columns, dt1.columns):
        dt.loc[idx, col0] = np.nan
        dt.loc[idx, col1] = np.nan
    dt = dt.interpolate(method='linear', limit_direction='both')
    return dt


def correct_wrongs(dt, c1c2):
    if c1c2 == 'mci_ad':
        c0 = 1
        c1 = 0
    else:
        c0 = 0
        c1 = 1

    zero_wrong = dt[(dt['0'] == 0.0) & (dt['class'] == c0)].index
    one_wrong = dt[(dt['1'] == 0.0) & (dt['class'] == c1)].index

    dt = fill_data(dt, one_wrong)
    dt = fill_data(dt, zero_wrong)

    return dt


def read_input(left, right, ft, train=True):

    if ft != features[0]:
        left = left.replace(features[0],
                            ft)
        right = right.replace(features[0],
                              ft)

    dataL = pd.read_csv(left, sep='\t')
    dataR = pd.read_csv(right, sep='\t')
    if train:
        dataL = correct_wrongs(dataL, c1c2)
        dataR = correct_wrongs(dataR, c1c2)

    data = dataL.merge(dataR, on="subj", suffixes=["R", "L"])
    left_cols = rename_columns(cols, 'L')[:-1]
    right_cols = rename_columns(cols, 'R')[1:]
    left_cols.extend(right_cols)
    selected_cols = left_cols
    renamed_cols = [name.replace('R', 'D')
                        .replace('L', 'E')
                        .replace('classD', 'class')
                    for name in selected_cols]

    data = data[selected_cols]
    data.columns = renamed_cols

    return data


def get_adni_data():
    adni = pd.read_csv("imgs_data.csv")
    adni = adni[['Subject ID', 'Age', 'Sex']]
    adni.columns = ['subj', 'age', 'sex']
    adni = adni[['subj', 'age']]
    return adni


def get_left(dt):
    return dt.iloc[:, dt.columns.str.contains("E")]


def get_right(dt):
    return dt.iloc[:, dt.columns.str.contains("D")]


def get_X_Y(train, test):
    X_train = train.iloc[:, ~train.columns.str.contains('class')]
    y_train = train.iloc[:, train.columns.str.contains(
        'class')].values.ravel()

    X_test = test.iloc[:, ~test.columns.str.contains('class')]
    y_test = test.iloc[:, test.columns.str.contains(
        'class')].values.ravel()

    return X_train, y_train, X_test, y_test


def rename_columns(col, prefix):
    new_name = []
    for col_name in col:
        if 'subj' not in col_name:
            new_name.append(col_name + prefix)
        else:
            new_name.append(col_name)
    return new_name


def concat_features(left, right):
    train = read_input(left, right, features[0], train=True)
    test = read_input(left.replace("validation", "test"),
                      right.replace("validation", "test"),
                      features[0], train=False)
    X_train, y_train, X_test, y_test = get_X_Y(train, test)
    X_train.columns = rename_columns(X_train.columns, "_" + features[0])
    X_test.columns = rename_columns(X_test.columns, "_" + features[0])

    # Concat all features in one data
    for ft in features[1:]:
        train = read_input(left, right, ft, train=True)
        test = read_input(left.replace("validation", "test"),
                          right.replace("validation", "test"),
                          ft, train=False)
        tmp_X_train, _, tmp_X_test, _ = get_X_Y(train, test)
        tmp_X_train.columns = rename_columns(tmp_X_train.columns, "_" + ft)
        tmp_X_test.columns = rename_columns(tmp_X_test.columns, "_" + ft)
        X_train = X_train.merge(tmp_X_train, on=['subj'])
        X_test = X_test.merge(tmp_X_test, on=['subj'])

    X_train.columns = rename_columns(X_train.columns, "_hippocampus")
    X_test.columns = rename_columns(X_test.columns, "_hippocampus")

    train = pd.concat([X_train,
                       pd.DataFrame(y_train, columns=['class'])
                       ], axis=1)
    test = pd.concat([X_test,
                      pd.DataFrame(y_test, columns=['class'])
                      ], axis=1)

    return train, test


if __name__ == "__main__":

    path = os.path.join("predicted_descriptors_validation", "gm")
    cols = ['subj', '0', '1', 'p0', 'p1', 'N', 'class']
    features = ['gm', 'wm', 'csf', 'tissues']
    output = 'image_train_test_data'

    for fold in range(10):
        for c1c2 in ['cn_ad', 'cn_mci', 'mci_ad']:
            fold = str(fold)
            left = os.path.join("../", fold, path, c1c2, "data_L.csv")
            right = os.path.join("../", fold, path, c1c2, "data_R.csv")
            out_prefix = fold + "_" + c1c2

            train, test = concat_features(left, right)

            print(train.head(), test.shape)

            train.to_csv(os.path.join(output,
                                      "fts_train_" + out_prefix + ".csv"),
                         index=False)

            test.to_csv(os.path.join(output,
                                     "fts_test_" + out_prefix + ".csv"),
                        index=False)

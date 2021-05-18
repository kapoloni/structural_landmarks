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
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler
from joblib import load
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score


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


def save_csv(data, name):
    data.to_csv(name, sep='\t', encoding='utf-8', index=False, header=False)


def save_csv_header(data, name):
    data.to_csv(name, sep='\t', encoding='utf-8', index=False)


def save_csv_desc(data, name):
    data.to_csv(name, sep=' ', encoding='utf-8', index=False, header=False)

def get_directories_class(csv_file, c1, c2):
    reader = pd.read_csv(csv_file, sep=",", header=None)
    reader.columns = ['subj', 'label']
    c1 = reader[reader.label == c1]
    c2 = reader[reader.label == c2]
    return pd.concat([c1.subj, c2.subj]).reset_index(drop=True)

def create_folder_and_delete(name):
    if os.path.exists(name):
        os.system("rm -rf " + name)
        os.makedirs(name, exist_ok=True)
    if not os.path.exists(name):
        os.makedirs(name, exist_ok=True)


def get_directories(root_folder, new_folder, split_):
    list_in, list_out = [], []
    for class_ in os.listdir(root_folder):
        if os.path.isdir(root_folder + "/" + class_):
            for filename in os.listdir(root_folder + "/" + class_):
                if split_ in filename:
                    file_root = root_folder+"/"+class_+"/"+filename.split(split_)[0]
                    file_out = new_folder+"/"+class_+"/"+filename.split(split_)[0]
                    list_in.append(file_root)
                    list_out.append(file_out)
    return list_in, list_out


def get_directories_fold(csv_file, root_path):
    reader = pd.read_csv(csv_file, sep='\t', header=None)
    return reader


def get_directories_fold_class(csv_file, root_path, c1, c2):
    reader = pd.read_csv(csv_file, sep="\t", header=None)
    dt = []
    for v in range(0, reader.shape[0]):
        dt.append(reader.iloc[v, np.where(np.logical_or(
                                reader.iloc[v].str.contains(c1+"/"),
                                reader.iloc[v].str.contains(c2+"/")))[0]
                                ].T.reset_index(drop=True).T)
    return pd.DataFrame(dt)


def create_folders(name):
    if not os.path.exists(name):
        os.makedirs(name, exist_ok=True)


def delete_document(name):
    if os.path.exists(name):
        os.remove(name)


def split_data(X, y):
    sss = StratifiedShuffleSplit(n_splits=10, test_size=0.1,
                                    random_state=0)  # 90 train 10 test_
    train_, test_ = [], []
    for train_index, test_index in sss.split(X, y):
        x_train, x_test = X.iloc[train_index], X.iloc[test_index]
        train_.append(x_train.values)
        test_.append(x_test.values)
    return test_, train_


def organize_data_3_class(data_):
    dt = pd.DataFrame(data_)
    dt = dt[0].str.split("/", expand=True)
    dt[5] = dt[[3, 4]].apply(lambda x: '/'.join(x), axis=1)
    return dt[5], dt[3]


def organize_data_(data_, c1, c2):
    dt = pd.DataFrame(data_)
    dt = dt[0].str.split("/", expand=True)
    dt = dt.iloc[np.where(dt[3].isin([c1, c2]))].reset_index()
    dt[5] = dt[[3, 4]].apply(lambda x: '/'.join(x), axis=1)
    return dt[5], dt[3]


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


def concat_match(input, c1, c2):
    c1_ = pd.DataFrame()
    c2_ = pd.DataFrame()
    for v in range(input.shape[0]):
        if type(input.iloc[v]) is str:
            if not os.stat(input.iloc[v]+".csv").st_size == 0:
                if "/"+c1+"/" in input.iloc[v]:
                    c1_ = pd.concat([c1_, count_match(input.iloc[v])])
                if "/"+c2+"/" in input.iloc[v]:
                    c2_ = pd.concat([c2_, count_match(input.iloc[v])])
    print(c1_)
    uniq, count = np.unique(c1_[[0, 1, 2]], return_counts=True, axis=0)
    c1_ = pd.concat([pd.DataFrame(uniq, columns=[0, 1, 2]),
                     pd.DataFrame(count, columns=[3])], axis=1)

    uniq, count = np.unique(c2_[[0, 1, 2]], return_counts=True, axis=0)
    c2_ = pd.concat([pd.DataFrame(uniq, columns=[0, 1, 2]),
                    pd.DataFrame(count, columns=[3])], axis=1)

    return c1_, c2_

def structural_conc(c1, c2, th):
    # c1[3] = c1[3].div(c1[3].max())
    # c2[3] = c2[3].div(c2[3].max())

    structural = c1.merge(c2, on=[0, 1, 2], suffixes=('_c1', '_c2'))
    structural['diff'] = structural['3_c2'] - structural['3_c1']
    structural['diff2'] = structural['diff']
    structural['diff'] = structural['diff'].div(structural['diff'].max())
    # print(structural['diff'], len(structural['diff']))
    structural.to_csv("diff"+str(th)+".csv")
    structural = structural.drop(structural.index[np.where(
                                    structural['diff'] < np.float(th))])
    return structural


def structural_conc_th(c1, c2, th):
    c1[3] = c1[3].div(c1[3].max())
    c2[3] = c2[3].div(c2[3].max())

    c1_c2 = c1.merge(c2, on=[0, 1, 2], suffixes=('_c1', '_c2'))
    c1_c2['diff'] = c1_c2['3_c2'] - c1_c2['3_c1']
    structural = c1_c2

    return structural


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


def generate_dataset(input, struc):
    data = []
    for v in range(0, input.shape[0]):
        if type(input[v]) is str:
            posic, desc = get_position_descriptor_txt(input[v])
            desc = desc.T.reset_index(drop=True).T
            if struc is True:
                desc['label'] = input[v].split("/")[-2]
            else:
                desc['label'] = "CC"
            data.append(desc)
            # print(input[v], input[v].split("/")[-2], data)
    dt = pd.concat(data, sort=True)
    return dt


def generate_dataset_cluster(input, atlas, struc):
    data = []
    for v in range(0, input.shape[0]):
        if type(input[v]) is str:
            file_atlas = pd.read_csv(atlas[v], header=None)
            file_atlas = file_atlas.iloc[:, :3].T.reset_index(drop=True).T
            file_atlas = file_atlas.iloc[:, :3].reset_index(drop=True)
            posic, desc = get_position_descriptor(input[v])
            desc = desc.T.reset_index(drop=True).T
            desc = desc.reset_index(drop=True)
            land = pd.concat([file_atlas, desc], axis=1)
            if struc is True:
                land['label'] = input[v].split("/")[-2]
            else:
                land['label'] = "CC"

            data.append(land)
    dt = pd.concat(data, sort=True)
    return dt


def predict_(X, y, train_, c1, c2):
    X = StandardScaler().fit_transform(X.astype(np.float64))
    clf = load(train_)
    result = clf.predict(X)
    class_, counts = np.unique(result, return_counts=True)
    ad = np.where(c2 == class_)
    cn = np.where(c1 == class_)
    ad = 0 if np.shape(ad)[1] == 0 else int(counts[ad[0]])
    cn = 0 if np.shape(cn)[1] == 0 else int(counts[cn[0]])
    return ad, cn


def predict_probability(X, y, train_, c1, c2):
    X = StandardScaler().fit_transform(X.astype(np.float64))
    clf = load(train_)
    result = pd.DataFrame(clf.predict_proba(X), columns=clf.classes_)
    result['true_class'] = y
    result['pred'] = ''

    for i in range(result.shape[0]):
        if result[c2].iloc[i] > result[c1].iloc[i]:
            if result[c2].iloc[i] > 0.6:
                result.loc[i, 'pred'] = c2
        else:
            if result[c1].iloc[i] > 0.6:
                result.loc[i, 'pred'] = c1

    result = result[result['pred'] != '']

    class_, counts = np.unique(result['pred'], return_counts=True)

    ad = np.where(c2 == class_)
    cn = np.where(c1 == class_)
    ad = 0 if np.shape(ad)[1] == 0 else int(counts[ad[0]])
    cn = 0 if np.shape(cn)[1] == 0 else int(counts[cn[0]])
    return ad, cn, np.sum(result[c1]), np.sum(result[c2])


def count_brains(desc, classifier, n_fold, c1, c2):
    left, right = [], []
    input = desc + str("_R_"+str(n_fold)+"_landmarks.txt")
    in_ = input.iloc[n_fold]
    for v in range(0, len(in_)):
        if type(in_[v]) is str:
            imgs_ = {}
            reader = pd.read_csv(in_[v], sep=' ', header=None)
            X = reader.iloc[1:, 3:-1]
            y = in_[v].split("/")[-2]
            a, b = predict_(X, y,
                                    classifier+"/rf_R_"+str(n_fold)+".pkl", c1, c2)
            subj = in_[v].split("/")[-1].split("_R")[0]
            imgs_['subj'], imgs_['a'], imgs_['b'], imgs_['class'] = subj, a, b, y
            right.append(imgs_)
    dataR = pd.DataFrame(right)

    input = desc + str("_L_"+str(n_fold)+"_landmarks.txt")
    in_ = input.iloc[n_fold]
    for v in range(0, len(in_)):
        if type(in_[v]) is str:
            imgs_ = {}
            reader = pd.read_csv(in_[v], sep=' ', header=None)
            X = reader.iloc[1:, 3:-1]
            y = in_[v].split("/")[-2]
            a, b = predict_(X, y,
                                    classifier+"/rf_L_"+str(n_fold)+".pkl", c1, c2)
            subj = in_[v].split("/")[-1].split("_L")[0]
            imgs_['subj'], imgs_['a'], imgs_['b'], imgs_['class'] = subj, a, b, y
            left.append(imgs_)
    dataL = pd.DataFrame(left)

    data = dataL.merge(dataR, on=['subj', 'class'], suffixes=('_L', '_R'))
    data['suma'] = data['a_L'] + data['a_R']
    data['sumb'] = data['b_L'] + data['b_R']
    return data


def pre_proc(data):
    na = data.isna().sum(axis=1)
    na = na.sort_values(ascending=False)
    data = data.iloc[na[na < 110].index.values]
    data = data.fillna(0)
    return data


def count_brains_texture(desc, classifier, n_fold, c1, c2, feature):
    left, right = [], []
    input = desc + str("_R_"+str(n_fold)+"_texture_desc.txt")
    in_ = input.iloc[n_fold]
    for v in range(0, len(in_)):
        if type(in_[v]) is str:
            imgs_ = {}
            reader = pd.read_csv(in_[v], sep=' ', header=None)
            # se a ultima coluna for nan
            if np.shape(np.where(
                        reader.iloc[1:, -1].isna()))[1] == int(reader.iloc[0, 0]):
                reader = reader.iloc[:, :-1]
            reader = pre_proc(reader)
            if feature == 'offset':
                X = reader.iloc[:, 3:]
            if feature == 'glcm':
                X = reader.iloc[:, 3:reader.shape[1]-104]
            if feature == 'shape':
                X = reader.iloc[:, 3:reader.shape[1]-104-8]
            y = in_[v].split("/")[-2]
            a, b, p_a, p_b = predict_probability(X, y,
                                                        classifier+"/rf_R_" +
                                                        str(n_fold) + ".pkl", c1, c2)
            subj = in_[v].split("/")[-1].split("_R")[0]
            imgs_['subj'], imgs_['a'], imgs_['b'], imgs_['class'] = subj, a, b, y
            imgs_['pa'], imgs_['pb'] = p_a, p_b
            right.append(imgs_)
    dataR = pd.DataFrame(right)

    input = desc + str("_L_"+str(n_fold)+"_texture_desc.txt")
    in_ = input.iloc[n_fold]
    for v in range(0, len(in_)):
        if type(in_[v]) is str:
            imgs_ = {}
            reader = pd.read_csv(in_[v], sep=' ', header=None)
            # se a ultima coluna for nan
            if np.shape(np.where(
                        reader.iloc[1:, -1].isna()))[1] == int(reader.iloc[0, 0]):
                reader = reader.iloc[:, :-1]
            reader = pre_proc(reader)
            if feature == 'offset':
                X = reader.iloc[:, 3:]
            if feature == 'glcm':
                X = reader.iloc[:, 3:reader.shape[1]-104]
            if feature == 'shape':
                X = reader.iloc[:, 3:reader.shape[1]-104-8]
            y = in_[v].split("/")[-2]
            a, b, p_a, p_b = predict_probability(X, y,
                                                        classifier+"/rf_L_" +
                                                        str(n_fold) + ".pkl", c1, c2)
            subj = in_[v].split("/")[-1].split("_L")[0]
            imgs_['subj'], imgs_['a'], imgs_['b'], imgs_['class'] = subj, a, b, y
            imgs_['pa'], imgs_['pb'] = p_a, p_b
            left.append(imgs_)
    dataL = pd.DataFrame(left)

    data = dataL.merge(dataR, on=['subj', 'class'], suffixes=('_L', '_R'))

    data['suma'] = data['a_L'] + data['a_R']
    data['sumb'] = data['b_L'] + data['b_R']
    data['sum_pb'] = data['pa_L'] + data['pb_L']
    data['sum_pa'] = data['pa_R'] + data['pb_R']
    return data


def transform_pca(X, input):
    pca_ = load(input)
    return pca_.transform(X)


def count_brains_pca(desc, classifier, n_fold, c1, c2):
    left, right = [], []
    input = desc + str("_R_"+str(n_fold)+"_landmarks.txt")
    in_ = input.iloc[n_fold]
    for v in range(0, len(in_)):
        if type(in_[v]) is str:
            imgs_ = {}
            reader = pd.read_csv(in_[v], sep=' ', header=None)
            X = reader.iloc[1:, 3:-1]
            y = in_[v].split("/")[-2]
            X = transform_pca(X, classifier+"/rf_R_"+str(n_fold)+"_pca.pkl")
            a, b = predict_(X, y,
                                    classifier+"/rf_R_"+str(n_fold)+".pkl", c1, c2)
            subj = in_[v].split("/")[-1].split("_R")[0]
            imgs_['subj'], imgs_['a'], imgs_['b'], imgs_['class'] = subj, a, b, y
            right.append(imgs_)
    dataR = pd.DataFrame(right)

    input = desc + str("_L_"+str(n_fold)+"_landmarks.txt")
    in_ = input.iloc[n_fold]
    for v in range(0, len(in_)):
        if type(in_[v]) is str:
            imgs_ = {}
            reader = pd.read_csv(in_[v], sep=' ', header=None)
            X = reader.iloc[1:, 3:-1]
            y = in_[v].split("/")[-2]
            X = transform_pca(X, classifier+"/rf_L_"+str(n_fold)+"_pca.pkl")
            a, b = predict_(X, y,classifier+"/rf_L_"+str(n_fold)+".pkl", c1, c2)
            subj = in_[v].split("/")[-1].split("_L")[0]
            imgs_['subj'], imgs_['a'], imgs_['b'], imgs_['class'] = subj, a, b, y
            left.append(imgs_)
    dataL = pd.DataFrame(left)

    data = dataL.merge(dataR, on=['subj', 'class'], suffixes=('_L', '_R'))
    data['suma'] = data['a_L'] + data['a_R']
    data['sumb'] = data['b_L'] + data['b_R']
    return data


def t_parameter(pred_folder, c1, c2):
    final = []
    th = 0
    while th < 5.05:
        acc, sens, spec = [], [], []
        dt = {}
        for v in range(0, 10):
            reader = pd.read_csv(pred_folder+"/data_"+str(v)+".csv", sep='\t')
            data = reader.iloc[:, np.where(reader.columns != 'class')[0]]
            label = reader['class']
            data['div'] = data['sumb'] / data['suma']
            data['pred'] = 'AD'
            data.loc[np.where(data['div'] > th)[0], 'pred'] = c1
            data.loc[np.where(data['div'] <= th)[0], 'pred'] = c2
            tn, fp, fn, tp = confusion_matrix(label, data['pred']).ravel()
            sens.append(tp / (tp+fn))
            spec.append(tn / (tn+fp))
            acc.append(accuracy_score(label, data['pred']))

        dt["th"] = th
        dt["Accuracy"] = np.mean(acc)
        dt["Sensitivity"] = np.mean(sens)
        dt["Specifity"] = np.mean(spec)
        final.append(dt)
        th = th + 0.05
    data = pd.DataFrame(final)
    return data

#!/usr/bin/python3
"""
    Save all classification models  (pickles) and metrics
"""
# Standard library imports
import numpy as np
import os
import pandas as pd
import pickle
import argparse
import sys

# Third party imports
from sklearn import svm
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import balanced_accuracy_score


def parse_args(args):
    """!@brief
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description='Classification Image level')
    parser.add_argument('--dest_folder', help='Destination folder',
                        default='../experiment_cmpb', type=str)
    return parser.parse_args(args)


def create_folders(name):
    if not os.path.exists(name):
        os.makedirs(name, exist_ok=True)


class ColumnExtractor(object):
    """
        Class to perform column extraction/selection
        in sklearn pipeline.
    """
    def __init__(self, str_):
        self.str_ = str_

    def transform(self, X):
        return X.iloc[:, X.columns.str.contains(self.str_)]

    def fit(self, X, y=None):
        return self


def metrics(y_pred, y_test, y_score, tauc):
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    bacc = balanced_accuracy_score(y_test, y_pred)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    auc_ = roc_auc_score(y_test, y_score)
    sensitivity = tp / (tp+fn)
    specificity = tn / (tn+fp)
    results = {'train': tauc,
               'f1': f1,
               'acc': accuracy,
               'bacc': bacc,
               'sens': sensitivity,
               'spe': specificity,
               'auc': auc_}
    return results


def get_columns_names():
    return ['gm', 'wm', 'csf', 'tissues']


def coarse_grid(ft):
    skf = StratifiedKFold(n_splits=5)
    pipe = Pipeline([("sel", ColumnExtractor(str_=ft)),
                     ("scale", StandardScaler()),
                     ('clf', svm.SVC(gamma='scale',
                                     kernel='poly',
                                     class_weight='balanced',
                                     probability=True))])

    param_grid = {"clf__C": [2**i for i in np.arange(-5, 10, 0.25)],
                  "clf__degree": [1]}

    clf = GridSearchCV(pipe,
                       param_grid=param_grid,
                       cv=skf,  n_jobs=-1,
                       scoring='roc_auc')
    return clf


def finer_grid(clf, ft):
    skf = StratifiedKFold(n_splits=5)
    pipe = Pipeline([("sel", ColumnExtractor(str_=ft)),
                     ("scale", StandardScaler()),
                     ('clf', svm.SVC(gamma='scale',
                                     kernel='poly',
                                     class_weight='balanced',
                                     probability=True))])

    c_exp = np.log2(np.float(clf.best_params_['clf__C']))
    c_par1 = [2**i for i in np.arange(c_exp-2, c_exp+2.1, 0.125)]
    param_grid2 = {}
    param_grid2['clf__C'] = c_par1
    param_grid2['clf__degree'] = [clf.best_params_['clf__degree']]

    clf = GridSearchCV(pipe, param_grid=param_grid2,
                       cv=skf, n_jobs=-1,
                       scoring='roc_auc')
    return clf


def get_X_Y(data):
    X = data.iloc[:, ~(data.columns.str.contains('class') |
                       data.columns.str.contains('subj'))]
    y = data.iloc[:, data.columns.str.contains('class')].values.ravel()
    return X, y


def fit_all_classifiers(df_train, df_test):
    X_train, y_train = get_X_Y(df_train)
    X_test, y_test = get_X_Y(df_test)
    classifiers = []
    for ft in get_columns_names():
        classifier = {}
        # Coarse grid classifier
        classifier['fts'] = ft
        classifier['clf0'] = coarse_grid(ft)
        classifier['clf0'].fit(X_train, y_train)
        # Finer grid classifier
        classifier['clf'] = finer_grid(classifier['clf0'], ft)
        classifier['clf'].fit(X_train, y_train)
        # Predictions
        y_pred = classifier['clf'].predict(X_test)
        y_score = classifier['clf'].predict_proba(X_test)[:, 1]
        scores = metrics(y_pred, y_test, y_score,
                         classifier['clf'].best_score_)
        classifier['scores'] = scores
        classifiers.append(classifier)

    return classifiers


def one_fold(folder, exp, fold, result_folder):
    name_train = 'fts_train.csv'
    name_test = 'fts_test.csv'
    df_train = pd.read_csv(os.path.join(folder, name_train))
    df_test = pd.read_csv(os.path.join(folder, name_test))

    result = fit_all_classifiers(df_train, df_test)
    outfile = open('{}/clf.pickle'.
                   format(result_folder), 'wb')

    pickle.dump(result, outfile)
    outfile.close()

    return result


def all_folds(exps, folds):

    list_results = []
    for exp in exps:
        print(exp)
        for i in folds:
            input_ = os.path.join(input_folder.replace("fold", str(i)), exp)
            result_ = os.path.join(result_folder.replace("fold", str(i)), exp)
            create_folders(result_)

            result = one_fold(input_, exp, i, result_)
            list_results.append({'exp': exp, 'res': result})

    return list_results


if __name__ == "__main__":

    # Parse arguments
    args = sys.argv[1:]
    args = parse_args(args)
    print(args)

    input_folder = os.path.join(args.dest_folder, "hippocampus",
                                "fold", "image_train_test_data")
    result_folder = os.path.join(args.dest_folder, "hippocampus",
                                 "fold", 'image_results')

    exps = ['mci_ad']  # 'cn_ad', 'cn_mci',
    fts = ['gm', 'wm', 'csf', 'tissues']
    folds = range(10)

    all_folds(exps, folds)

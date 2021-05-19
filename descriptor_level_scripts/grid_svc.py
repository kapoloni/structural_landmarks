#!/usr/bin/python3
"""
SVM grid search classifier for the descriptors level
"""

# Standard library imports
import pandas as pd
import numpy as np
import os
import getopt
import sys

# Third party imports
import multiprocessing
from sklearn import svm
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import make_scorer
from sklearn.model_selection import GroupKFold
from joblib import dump


def showUsage():
    print('./grid_svc.py -n <result_folder>\
          -c -d')
    sys.exit(2)


def get_params(argv):
    try:
        opts, args = getopt.getopt(argv,
                                   "hn:c:d:f:t:", ["nfile=",
                                                   "c1file=",
                                                   "c2file=",
                                                   "fnum="])
    except getopt.GetoptError:
        showUsage()
    for opt, arg in opts:
        if opt == '-h':
            showUsage()
        elif opt in ("-n", "--nfile"):
            result_folder_ = arg
        elif opt in ("-c", "--c1file"):
            c1_ = arg
        elif opt in ("-d", "--c2file"):
            c2_ = arg
        elif opt in ("-f", "--fnum"):
            fold = arg
    return result_folder_, c1_, c2_, fold


def multip(params_, func, n_proc):
    p = multiprocessing.Pool(n_proc)
    p.map(func, params_)
    p.close()


def save_csv(data, name):
    data.to_csv(name, sep='\t', encoding='utf-8', index=False)


def create_folders(name):
    if not os.path.exists(name):
        os.makedirs(name, exist_ok=True)


def getSens(estimator, x, y):
    yPred = estimator.predict(x)
    _, _, fn, tp = confusion_matrix(y, yPred).ravel()
    sensitivity = tp / (tp+fn)
    return (sensitivity)


def getSpec(estimator, x, y):
    yPred = estimator.predict(x)
    tn, fp, _, _ = confusion_matrix(y, yPred).ravel()
    specificity = tn / (tn+fp)
    return (specificity)


def get_descriptors(X):
    print(X.shape)

    attributes = {}
    desc_size = 96
    # gm
    attributes['gm'] = X.iloc[:, :desc_size].astype(
        np.float64).reset_index(drop=True)
    # wm
    attributes['wm'] = X.iloc[:, desc_size:desc_size * 2].astype(
        np.float64).reset_index(drop=True)
    # csf
    attributes['csf'] = X.iloc[:, desc_size * 2:].astype(
        np.float64).reset_index(drop=True)

    attributes['tissues'] = pd.concat([attributes['gm'],
                                       attributes['wm'],
                                       attributes['csf']],
                                      axis=1).reset_index(drop=True)
    return attributes


def classifier_config():
    # Define classifier configs
    config = {"classifier": svm.SVC(kernel='rbf',
                                    gamma='scale',
                                    class_weight='balanced',
                                    probability=True),
              "params": {'classifier__C':
                         [2**i for i in np.arange(-5, 5, .5)]}
              }
    return config


def my_scorer_groups(y_true, y_pred, groups):

    y_pred = pd.Series(y_pred)
    y_pred.index = y_true.index
    groups = groups.iloc[y_true.index]
    dt = pd.concat([y_true, y_pred, groups], axis=1)
    dt.columns = ['true', 'pred', 'group']

    f1s = []
    for values in dt.groupby('group'):
        value = np.array(values)
        y_true, y_pred = value[1]['true'], value[1]['pred']
        if y_true.unique() == 0:
            y_true = pd.Series(np.where(y_true == 0, 1, 0))
            y_pred = pd.Series(np.where(y_pred == 0, 1, 0))
        f1 = f1_score(y_true, y_pred, average='weighted')
        f1s.append(f1)
    return np.mean(f1s)


def my_custom_loss_func(y_true, y_pred):
    return 0


def grid_search(X, y, njobs, output_, groups):
    my_score = make_scorer(my_scorer_groups, groups=groups)
    # my_score = make_scorer(my_custom_loss_func, groups=groups)
    # print(my_score)
    cv = GroupKFold(n_splits=5)
    print("##"*5, 'Grid Search', "##"*5)
    print(X.shape)
    X = pd.DataFrame(X)
    y = pd.Series(y)
    groups = pd.Series(groups)

    config = classifier_config()
    classifier = config['classifier']
    param_grid = config['params']

    steps = [("scaler", StandardScaler()),
             ("classifier", classifier)]
    pipe = Pipeline(steps=steps)

    # Coarse grid
    clf = GridSearchCV(estimator=pipe,
                       param_grid=param_grid,
                       cv=cv, n_jobs=njobs,
                       scoring=my_score)
    clf.fit(X, y, groups)

    # Finer grid
    c_exp = np.log2(np.float(clf.best_params_['classifier__C']))
    c_par = [2**i for i in np.arange(c_exp-2, c_exp+2.1, .25)]
    param_grid2 = param_grid
    param_grid2['classifier__C'] = c_par

    clf = GridSearchCV(estimator=pipe,
                       param_grid=param_grid2,
                       cv=cv, n_jobs=njobs,
                       scoring=my_score)
    clf.fit(X, y, groups=groups)

    dump(clf, output_)
    print(clf.best_params_, clf.best_score_)
    return clf.best_params_, clf.best_score_


def svm_grid(input, side):
    classifier_folder = os.path.join(result_folder, "hippocampus",
                                     fold, "classifiers")
    reader = pd.read_csv(input, sep='\t', header=None)
    y = reader.iloc[:, -1]
    groups = reader.iloc[:, -2]

    if c1c2 == "mci_ad":
        y = np.where(y == class1, 1, 0)
    else:
        y = np.where(y == class1, 0, 1)

    print("\nSide: ", side, "\nClass: ", np.unique(y, return_counts=True))
    attributes = get_descriptors(reader.iloc[:, :-2])

    results = []
    for tissue, X in attributes.items():
        result_dict = {}
        print(tissue, X.shape)

        output_folder = os.path.join(classifier_folder + tissue, c1c2)
        create_folders(output_folder)
        output = os.path.join(output_folder, "model_" + side + ".pkl")

        best, auc = grid_search(X, y, -1, output, groups)
        result_dict['best'], result_dict['auc'] = best, auc
        result_dict['feature'] = tissue
        results.append(result_dict)

    print(pd.DataFrame(results), "\n\n")
    print(results)
    save_csv(pd.DataFrame(results),
             output_folder+"output"+side+".csv")


# ********************************************************************
#
# ********************************************************************


if __name__ == "__main__":

    if len(sys.argv) < 5:
        showUsage()
        exit()

    result_folder, class1, class2, fold = get_params(sys.argv[1:])

    c1c2 = class1.lower()+"_"+class2.lower()

    basename = os.path.join(result_folder, "hippocampus", fold)

    descriptors = os.path.join(basename, "svm_descriptor_train", c1c2)
    output_folder = os.path.join(basename, "svm_results", c1c2)
    create_folders(output_folder)

    svm_grid(descriptors + "/dataset_L.csv", "L")
    svm_grid(descriptors + "/dataset_R.csv", "R")

#!/usr/bin/python3
import pandas as pd
import numpy as np
import multiprocessing
import os
import getopt
import sys
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
          -r <region (hippocampus, ventricles, brain)> \
          -c -d')
    sys.exit(2)


def get_params(argv):
    try:
        opts, args = getopt.getopt(argv,
                                   "hn:r:c:d:f:t:", ["nfile=", "rfile=",
                                                     "c1file=", "c2file=",
                                                     "fnum="])
    except getopt.GetoptError:
        showUsage()
    for opt, arg in opts:
        if opt == '-h':
            showUsage()
        elif opt in ("-n", "--nfile"):
            result_folder_ = arg
        elif opt in ("-r", "--rfile"):
            region_ = arg
        elif opt in ("-c", "--c1file"):
            c1_ = arg
        elif opt in ("-d", "--c2file"):
            c2_ = arg
        elif opt in ("-f", "--fnum"):
            fold = arg
    return result_folder_, region_, c1_, c2_, fold

# ********************************************************************
#
# ********************************************************************


def multip(params_, func, n_proc):
    p = multiprocessing.Pool(n_proc)
    p.map(func, params_)
    p.close()

# ********************************************************************
#
# *****************************************************************pipeline***


def save_csv(data, name):
    data.to_csv(name, sep='\t', encoding='utf-8', index=False)

# ********************************************************************
#
# ********************************************************************


def create_folders(name):
    if not os.path.exists(name):
        os.makedirs(name, exist_ok=True)

# ********************************************************************
#
# ********************************************************************


def getSens(estimator, x, y):
    yPred = estimator.predict(x)
    tn, fp, fn, tp = confusion_matrix(y, yPred).ravel()
    sensitivity = tp / (tp+fn)
    return (sensitivity)


def getSpec(estimator, x, y):
    yPred = estimator.predict(x)
    tn, fp, fn, tp = confusion_matrix(y, yPred).ravel()
    specificity = tn / (tn+fp)
    return (specificity)


def get_descriptors(X):
    # gm
    print(X.shape)
    desc_size = 96
    shape = X.iloc[:, :desc_size].astype(np.float64)
    # gm
    gm = X.iloc[:, desc_size:desc_size*2].astype(np.float64)
    # wm
    wm = X.iloc[:, desc_size*2:desc_size*3].astype(np.float64)
    # csf
    csf = X.iloc[:, desc_size*3:desc_size*4].astype(np.float64)
    # GLCM
    glcm = X.iloc[:, desc_size*4:-10].astype(np.float64)
    # RLM
    rlm = X.iloc[:, -8:].astype(np.float64)

    t1 = pd.concat([shape, glcm, rlm], axis=1)
    gm = pd.concat([gm], axis=1)
    wm = pd.concat([wm], axis=1)
    csf = pd.concat([csf], axis=1)

    return t1, gm, wm, csf

# ********************************************************************
#
# ********************************************************************


def combinations(t1, gm, wm, csf):
    combs = {}
    combs['t1'] = t1
    combs['gm'] = gm
    combs['wm'] = wm
    combs['csf'] = csf
    combs['tissues'] = pd.concat([gm, wm, csf], axis=1)
    combs['all'] = pd.concat([t1, gm, wm, csf], axis=1)
    return combs


# ********************************************************************
#
# ********************************************************************


def svm_grid(input, side):
    classifier_folder = result_folder+"/"+region + \
        "/" + fold + "/classifiers/"
    reader = pd.read_csv(input, sep='\t', header=None)
    y = reader.iloc[:, -1]
    groups = reader.iloc[:, -2]
    # print('y', y)
    # print('g', groups)
    if c1c2 == "mci_ad":
        y = np.where(y == class1, 1, 0)
    else:
        y = np.where(y == class1, 0, 1)

    print("\nSide: ", side, "\nClass: ", np.unique(y, return_counts=True))
    t1, gm, wm, csf = get_descriptors(reader.iloc[:, :-2])
    combs = combinations(t1, gm, wm, csf)

    results = []
    for k, v in combs.items():
        # if k != 'all':
        #     continue
        X = v
        result_dict = {}
        print(k, X.shape)

        create_folders(classifier_folder + k + "/" + c1c2)
        output = classifier_folder + k + "/" + c1c2 + \
            "/model_" + side + ".pkl"

        best, auc = grid_search(X, y, -1, output, groups)
        result_dict['best'], result_dict['auc'] = best, auc
        result_dict['feature'] = k
        results.append(result_dict)

    print(pd.DataFrame(results), "\n\n")
    print(results)
    save_csv(pd.DataFrame(results),
             output_folder+"output"+side+".csv")

# ********************************************************************
#
# ********************************************************************


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


def grid_search(X_, y_, njobs, output_, groups):
    my_score = make_scorer(my_scorer_groups, groups=groups)
    selected_classifier = "SVM-RBF"
    # cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=47)
    cv = GroupKFold(n_splits=5)
    print("##"*5, selected_classifier, "##"*5)
    print(X_.shape)
    X_ = pd.DataFrame(X_)
    y_ = pd.Series(y_)
    groups = pd.Series(groups)
    # print(groups)
    classifier = classifiers[selected_classifier]
    param_grid = parameters[selected_classifier]

    steps = [("scaler", StandardScaler()),
             ("classifier", classifier)]
    pipe = Pipeline(steps=steps)

    # Coarse grid
    clf = GridSearchCV(estimator=pipe,
                       param_grid=param_grid,
                       cv=cv, n_jobs=njobs,
                       scoring=my_score)
    clf.fit(X_, y_, groups=groups)
    # print(clf.best_params_, clf.best_score_)

    # Finer grid
    c_exp = np.log2(np.float(clf.best_params_['classifier__C']))
    c_par = [2**i for i in np.arange(c_exp-2, c_exp+2.1, .25)]
    param_grid2 = param_grid
    param_grid2['classifier__C'] = c_par

    clf = GridSearchCV(estimator=pipe,
                       param_grid=param_grid2,
                       cv=cv, n_jobs=njobs,
                       scoring=my_score)
    clf.fit(X_, y_, groups=groups)
    # print(clf.best_params_, clf.best_score_)

    dump(clf, output_)

    # scoring = {'acc': 'accuracy',
    #            'f1': 'f1',
    #            'sens': getSens,
    #            'spec': getSpec,
    #            'auc': 'roc_auc',
    #            'ap': 'average_precision'}
    # scores = cross_validate(clf, X_, y_, groups=groups, scoring=scoring,
    #                         cv=cv, return_train_score=True)

    # print(scores.keys())
    # print(scores['test_acc'])

    print(clf.best_params_, clf.best_score_)
    return clf.best_params_, clf.best_score_

# ********************************************************************
#
# ********************************************************************


if __name__ == "__main__":

    if len(sys.argv) < 5:
        showUsage()
        exit()

    result_folder, region, class1, class2, fold = get_params(sys.argv[1:])

    c1c2 = class1.lower()+"_"+class2.lower()

    # Create list of tuples with classifier label and classifier object
    classifiers, parameters = {}, {}
    classifiers.update({"SVM-RBF": svm.SVC(kernel='rbf',
                                           gamma='scale',
                                           class_weight='balanced',
                                           probability=True)})
    # classifiers.update({"SVM-poly": svm.SVC(kernel='poly',
    #                                         gamma='scale',
    #                                         class_weight='balanced',
    #                                         probability=True)})

    c_par = [2**i for i in np.arange(-5, 5, .5)]

    parameters.update({"SVM-poly": {
                       "classifier__C": c_par,
                       "classifier__class_weight": ['balanced'],
                       "classifier__degree": [1, 2, 3]}})

    parameters.update({"SVM-RBF": {"classifier__C": c_par}})

    # Variables
    descriptors = result_folder+"/"+region+"/" + \
        fold+"/final_desc/"+c1c2+"/"
    output_folder = result_folder+"/"+region+"/" + \
        fold+"/th_results/"+c1c2+"/"
    create_folders(output_folder)

    svm_grid(descriptors + "/dataset_L.csv", "L")
    svm_grid(descriptors + "/dataset_R.csv", "R")

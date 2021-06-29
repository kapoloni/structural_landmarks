#!/usr/bin/python3
"""
ROC curve plot
"""

# Standard library imports
import pandas as pd
import numpy as np
import pickle
import os

# Third party imports
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns


class ColumnExtractor(object):
    def __init__(self, str_):
        self.str_ = str_

    def transform(self, X):
        return X.iloc[:, X.columns.str.contains(self.str_)]

    def fit(self, X, y=None):
        return self


def get_results(exps):
    results = []
    for exp in exps:
        for fold in folds:
            path = os.path.join(result_folder.replace(
                "fold", str(fold)), exp, "clf.pickle")

            infile = open(path, 'rb')
            new_dict = pickle.load(infile)
            infile.close()
            for j in range(len(new_dict)):
                res = new_dict[j]['scores']
                res['exp'] = exp
                res['fold'] = fold
                res['ft'] = new_dict[j]['fts']
                res['clf'] = new_dict[j]['clf']

                results.append(res)
    return pd.DataFrame(results)


def get_roc_probs(df, exp, ft):
    auc, tprs, roc_values = [], [], []
    df = df[df['exp'] == exp]
    mean_fpr = np.linspace(0, 1, 100)
    for fold in range(10):
        # Best result
        df_fold = df[(df.fold == fold) & (df.ft == ft)]
        name_train = 'fts_train.csv'
        name_test = 'fts_test.csv'
        folder = os.path.join(input_folder.replace("fold", str(fold)), exp)
        df_train = pd.read_csv(os.path.join(folder, name_train))
        df_test = pd.read_csv(os.path.join(folder, name_test))

        X_train = df_train.iloc[:, 1:-1]
        y_train = df_train.iloc[:, -1]

        X_test = df_test.iloc[:, 1:-1]
        y_test = df_test.iloc[:, -1]

        c = df_fold['clf'].values[0].best_params_['clf__C']
        d = df_fold['clf'].values[0].best_params_['clf__degree']

        pipe = Pipeline([("sel", ColumnExtractor(str_='tissues')),
                        ("scale", StandardScaler()),
                        ('clf', svm.SVC(gamma='scale',
                                        kernel='poly',
                                        class_weight='balanced',
                                        probability=True,
                                        C=c, degree=d))])

        pipe.fit(X_train, y_train)

        y_score = pipe.predict_proba(X_test)[:, 1]
        auc_ = roc_auc_score(y_test, y_score)
        fpr, tpr, thresholds = roc_curve(y_test, y_score, pos_label=1)
        tprs.append(np.interp(mean_fpr, fpr, tpr))
        auc.append(auc_)
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    roc_values.append([np.linspace(0, 1, 100), mean_tpr, auc, tprs])
    return roc_values


def plot_roc(roc_values, save=False, pt=False):
    sns.set(context='paper', style='white',
            font='sans-serif', font_scale=3,
            rc={'figure.figsize': (12, 10)})

    mean_fpr, mean_tpr, aucs, tprs = roc_values['cn_mci'][0]

    if pt:
        name1 = 'CNxCCL'
        name2 = 'CCLxDA'
        name3 = 'CNxDA'
    else:
        name1 = 'CNxMCI'
        name2 = 'MCIxAD'
        name3 = 'CNxAD'
    plt.plot(mean_fpr, mean_tpr,
             label=r'%s = %0.2f' % (name1, np.mean(aucs)),
             lw=2)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, alpha=.1,)

    mean_fpr, mean_tpr, aucs, tprs = roc_values['mci_ad'][0]
    plt.plot(mean_fpr, mean_tpr,
             label=r'%s = %0.2f' % (name2, np.mean(aucs)),
             lw=2)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, alpha=.1,)

    mean_fpr, mean_tpr, aucs, tprs = roc_values['cn_ad'][0]

    plt.plot(mean_fpr, mean_tpr,
             label=r'%s = %0.2f' % (name3, np.mean(aucs)),
             lw=2)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, alpha=.1,)

    plt.xlim([0., 1.0])
    plt.ylim([0., 1.0])
    if pt:
        plt.xlabel('1-especificidade')
        plt.ylabel('Sensibilidade')
        outname = "plots/roc_final_pt.pdf"
    else:
        plt.xlabel('1-specificity')
        plt.ylabel('Sensitivity')
        outname = "plots/roc_final.pdf"
    plt.legend(loc="lower right", fontsize='x-small')
    if save:
        plt.tight_layout()
        plt.savefig(outname, bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    input_folder = os.path.join("..", "experiment_cmpb", "hippocampus",
                                "fold", "image_train_test_data")
    result_folder = os.path.join("..", "experiment_cmpb", "hippocampus",
                                 "fold", 'image_results')
    exps = ['cn_ad', 'cn_mci', 'mci_ad']
    fts = ['gm', 'wm', 'csf', 'tissues']
    folds = range(10)

    results = get_results(exps)

    roc_values = {}
    for exp in exps:
        roc_values[exp] = get_roc_probs(results, exp, ft='tissues')
    print(np.mean(roc_values['mci_ad'][0][3]))
    plot_roc(roc_values, save=True, pt=True)
    plot_roc(roc_values, save=True, pt=False)

    outfile = open('results/roc_cnad.pickle', 'wb')
    pickle.dump(roc_values['cn_ad'][0], outfile)
    outfile.close()

    outfile = open('results/roc_mciad.pickle', 'wb')
    pickle.dump(roc_values['mci_ad'][0], outfile)
    outfile.close()

    outfile = open('results/roc_cnmci.pickle', 'wb')
    pickle.dump(roc_values['cn_mci'][0], outfile)
    outfile.close()

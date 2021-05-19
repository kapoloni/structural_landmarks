#!/usr/bin/python3
"""
Print image results
"""

# Standard library imports
import pandas as pd
import pickle
import os


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

                results.append(res)
    return pd.DataFrame(results)


def print_mean_results(dt, exp, round=3):
    cols = ['auc', 'acc', 'sens', 'spe', 'f1', 'ft']
    dt = dt[dt['exp'] == exp][cols]
    means = (dt.groupby(['ft']).mean()*100).round(round)
    stds = (dt.groupby(['ft']).std(ddof=0)*100).round(round)
    metrics = means.columns.values
    fts = means.index.values
    results = means.copy()
    for metric in metrics:
        for ft in fts:
            if metric in ['auc', 'f1']:
                means.loc[ft, metric] = (means.loc[ft, metric]/100
                                         ).round(round)
                stds.loc[ft, metric] = (stds.loc[ft, metric]/100).round(round)
            results.loc[ft, metric] = str(means.loc[ft, metric]) + \
                " +- " + str(stds.loc[ft, metric])
    return results


def get_latex_table(dt, exp, round=3):

    fts_name = {'gm': 'GM',
                'wm': 'WM',
                'csf': 'CSF',
                'tissues': 'GM\\&WM\\&CSF'}

    cols = ['auc', 'acc', 'sens', 'spe', 'f1', 'ft']
    dt = dt[dt['exp'] == exp][cols]
    means = (dt.groupby(['ft']).mean()*100).round(round)
    stds = (dt.groupby(['ft']).std(ddof=0)*100).round(round)
    metrics = means.columns.values
    fts = means.index.values

    for ft in fts:
        text = "\\textbf{"+fts_name[ft]+"}"
        for metric in metrics:
            if metric in ['auc', 'f1']:
                means.loc[ft, metric] = (
                    means.loc[ft, metric]/100).round(round)
                stds.loc[ft, metric] = (stds.loc[ft, metric]/100).round(round)

            text += ' & $' + str(means.loc[ft, metric]).replace(".", ',') + \
                    "\\pm" + str(stds.loc[ft, metric]).replace(".", ',')+"$ "
        text += "\\tabularnewline"
        print(text)
        print("\\midrule")


if __name__ == "__main__":

    result_folder = os.path.join("..", "experiment_cmpb", "hippocampus",
                                 "fold", 'image_results')
    exps = ['cn_ad', 'cn_mci', 'mci_ad']
    fts = ['gm', 'wm', 'csf', 'tissues']
    folds = range(10)

    results = get_results(exps)
    for exp in exps:
        print(exp)
        print(print_mean_results(results, exp, round=2))

    # Latex print
    for exp in exps:
        print(exp)
        print(get_latex_table(results, exp, round=2))

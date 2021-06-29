#!/usr/bin/python3
"""
Get best and worst subjects results
"""

# Standard library imports
import os
import pandas as pd


if __name__ == "__main__":

    input_folder = os.path.join('..', 'experiment_cmpb', "hippocampus",
                                "fold", "image_results")

    exps = ['cn_ad', 'cn_mci', 'mci_ad']
    # Folds belonging to the graphs
    folds = {'cn_ad': [6, 7],
             'cn_mci': [2, 3],
             'mci_ad': [5, 6]}
    for exp in exps:
        print("#"*5, exp, "#"*5)
        for i in range(folds[exp][0], folds[exp][1]):
            input_ = os.path.join(input_folder.replace("fold", str(i)), exp)
            preds = pd.read_csv('{}/predictions.csv'.format(input_))
            preds = preds[['Subject', 'True', 'Pred_tissues',
                           'Score_tissues']]
            # preds['pred_sc'] = preds['Score_tissues'].apply(lambda x: 0 if x < 0.5 else 1)
            preds['Corr_score'] = preds['Score_tissues'].apply(lambda x: 1-x if x < 0.5 else x)
            # preds['Abs_df'] = preds['Score_tissues'].apply(lambda x: abs(x))
            for label in [0, 1]:
                label_frac = preds[(preds['True'] == label) &
                                   (preds['Pred_tissues'] == label)]
                best = label_frac[(label_frac['True'] ==
                                   label_frac.Pred_tissues)
                                  ]['Corr_score'].argmax()
                # print("Best", label, label_frac.iloc[best])
            
                print("Best", label, label_frac.iloc[best][[
                      'Subject', 'Corr_score']].values)
                worst = label_frac[(label_frac['True'] ==
                                    label_frac.Pred_tissues)
                                   ]['Corr_score'].argmin()
                # print("Worst", label, label_frac.iloc[worst])
                print("Worst", label, label_frac.iloc[worst][[
                      'Subject', 'Corr_score']].values)

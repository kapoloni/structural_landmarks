#!/usr/bin/python3
import os
import pandas as pd


def read_params(c1c2, side, fold, ft):
    folder = os.path.join(result_folder, 'hippocampus',
                          fold, "grid_search",
                          'tissues', c1c2,
                          "output" + side + ".csv")
    reader = pd.read_csv(folder, sep='\t')
    idx = reader[reader['feature'] == ft].index
    best_params_ = reader.loc[idx, 'best'].values[0]
    c = dict(eval(best_params_))['classifier__C']

    return str(c)


if __name__ == "__main__":

    result_folder = os.path.join("..", "experiment")

    features = ['gm', 'wm', 'csf', 'tissues']

    # Split database
    if not os.path.exists(os.path.join(result_folder, 'splits')):
        os.system("./train_test_split.py -n " + result_folder)

    # Matches
    os.system("./match_all_landmarks.py -n " + result_folder)

    for exp in [['CN', 'AD']]:  # , ['CN', 'MCI'], ['MCI', 'AD']
        for fold in range(0, 10):
            fold = str(fold)
            c1c2 = exp[0].lower() + "_" + exp[1].lower()
            # for side in ["L", "R"]:
            #     print(side, fold)
            #     # Get structural points/descriptors
            #     os.system("./structural_descriptors.py --dest_folder " + result_folder +
            #               " --labels " + c1c2 + " --side " + side + " --fold " + fold)

            #     # Get attributes points/descriptors
            #     os.system("./generate_descriptors_attributes.py --dest_folder " +
            #               result_folder + " --labels " + c1c2 +
            #               " --side " + side + " --fold " + fold)

            # # Grid search
            # os.system("./grid_svc.py --dest_folder " + result_folder +
            #           " --labels " + c1c2 + " --fold " + fold)
            # exit()
            for ft in features:
                for side in ["R", "L"]:
                    C = read_params(c1c2, side, fold, ft)
                    # Classifiers
                    os.system("./fit_classifier.py --dest_folder " +
                              result_folder + " --labels " + c1c2 +
                              " --side " + side + " --fold " + fold +
                              " --ft " + ft + " --c " + C)
            #         # Predictions
            #         os.system("./predict_descriptors.py -n " +
            #                     result_folder + " -r " + region +
            #                     " -d " + exp[0] + " -e " + exp[1] +
            #                     " -s " + side + " -f " + ft + " -t " +
            #                     fold)

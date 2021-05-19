#!/usr/bin/python3
"""
Execute all pipeline for the image level
"""

# Standard library imports
import os


if __name__ == "__main__":

    result_folder = os.path.join("..", "experiment_cmpb")

    # Organize the extracted training and test attributes
    # of the descriptions and generate an input file
    # for the classifier at the image level.
    # os.system("./get_image_train_test_data.py " +
    #           "--dest_folder " + result_folder)

    # Save all classification models (pickles) and metrics
    os.system("./classification_concatenated_features.py " +
              "--dest_folder " + result_folder)

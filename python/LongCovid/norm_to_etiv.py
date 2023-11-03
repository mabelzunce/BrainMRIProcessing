import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

def normalize_to_etiv(value, etiv):
    # takes a value and a parameter and return normalized value (in%)
    return value/etiv * 100



if __name__ == '__main__':

    path_parcellation_volumes  = "/home/sol/COVID/CSV_subjects/parcellation.csv"
    path_segmentation_volumes  = "/home/sol/COVID/CSV_subjects/segmentation.csv"
    path_brainvol = "/home/sol/COVID/CSV_subjects/brain_volumes.csv"

    # Read csv
    parcellation_volumes = pd.read_csv(path_parcellation_volumes)
    segmentation_volumes = pd.read_csv(path_segmentation_volumes)
    brainvol_data = pd.read_csv(path_brainvol)

    # Initialize new CSV
    parcellation_volumes_to_etiv = pd.DataFrame(parcellation_volumes)
    segmentation_volumes_to_etiv = pd.DataFrame(segmentation_volumes)

    # # Rename: TO DO
    # brainvol_data.columns[0] = "Subject"
    # print(brainvol_data)

    # Parcellation columns to apply normalization to etiv
    parcellation_volumes_columns = list(parcellation_volumes.columns[1:-3])
    segmentation_volumes_columns = list(segmentation_volumes.columns[1:-2])


    for index, row in parcellation_volumes_to_etiv.iterrows():
        # row[0] == subject

        # Extract Estimated Total Intercranial Volume
        etiv = brainvol_data.loc[brainvol_data["subject"] == row[0]]["Estimated Total Intracranial Volume"].values[0]

        etiv = float(etiv)

        # Normalize to eTIV for selected column
        for column in parcellation_volumes_columns:
            parcellation_volumes_to_etiv.at[index, column] = normalize_to_etiv(float(row[column]), etiv)

    for index, row in segmentation_volumes_to_etiv.iterrows():
        # row[0] == subject

        # Extract Estimated Total Intercranial Volume
        etiv = brainvol_data.loc[brainvol_data["subject"] == row[0]]["Estimated Total Intracranial Volume"].values[0]
        etiv = float(etiv)

        # Normalize to eTIV for selected column
        for column in segmentation_volumes_columns:
            segmentation_volumes_to_etiv.at[index, column] = normalize_to_etiv(row[column], etiv)


    # Write CSV
    path_parcellation_to_etiv = "/home/sol/COVID/CSV_subjects/parcellation_etiv.csv"
    path_segmentation_to_etiv = "/home/sol/COVID/CSV_subjects/segmentation_etiv.csv"

    parcellation_volumes_to_etiv.to_csv(path_parcellation_to_etiv, index=False)
    segmentation_volumes_to_etiv.to_csv(path_segmentation_to_etiv, index=False)
    #
    # # Read Brainvol to extract ETIV value
    # etiv = brainvol_data.loc[brainvol_data["subject"] == subject]["group"].values[0]
    # etiv = CSV_covid_group.loc[CSV_covid_group["subject"] == subject]["group"].values[0]

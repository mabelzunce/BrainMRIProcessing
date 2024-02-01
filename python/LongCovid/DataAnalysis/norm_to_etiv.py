import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

def normalize_to_etiv(value, etiv):
    # takes a value and a parameter and return normalized value (in%)
    return value/etiv * 100

def parcellation_to_etiv(df_brainvol, df_parcellation):

    '''Normalize df parcellation to etiv values'''
    df_parcellation_to_etiv = pd.DataFrame(df_parcellation)

    # Parcellation columns to apply normalization to etiv
    parcellation_volumes_columns = list(df_parcellation_to_etiv.columns[1:-3])

    # Iterate over the rows and normalize values
    for index, row in df_parcellation_to_etiv.iterrows():

        # Extract Estimated Total Intercranial Volume
        subject = row['subject']
        etiv = df_brainvol.loc[df_brainvol["subject"] == subject]["Estimated Total Intracranial Volume"].values[0]
        etiv = float(etiv)

        # Normalize to eTIV for selected column
        for column in parcellation_volumes_columns:
            df_parcellation_to_etiv.at[index, column] = normalize_to_etiv(float(row[column]), etiv)

    return df_parcellation_to_etiv

def segmentation_to_etiv(df_brainvol, df_segmentation):

    '''Normalize df segmentation to etiv values'''
    df_segmentation_to_etiv = pd.DataFrame(df_segmentation)

    # Segmentation columns to apply normalization to etiv
    segmentation_volumes_columns = list(df_segmentation_to_etiv.columns[1:-2])

    # Iterate over the rows and normalize values
    for index, row in df_segmentation_to_etiv.iterrows():

        # Extract Estimated Total Intercranial Volume
        subject = row['subject']

        etiv = df_brainvol.loc[df_brainvol["subject"] == subject]["Estimated Total Intracranial Volume"].values[0]
        etiv = float(etiv)

        # Normalize to eTIV for selected column
        for column in segmentation_volumes_columns:
            df_segmentation_to_etiv.at[index, column] = normalize_to_etiv(row[column], etiv)

    return df_segmentation_to_etiv


if __name__ == '__main__':

    path_parcellation_volumes  = "/home/sol/COVID/CSV_subjects/parcellation.csv"
    path_parcellation_DKT_volumes = "/home/sol/COVID/CSV_subjects/parcellation_DKT.csv"
    path_segmentation_volumes  = "/home/sol/COVID/CSV_subjects/segmentation.csv"
    path_brainvol = "/home/sol/COVID/CSV_subjects/brain_volumes.csv"

    # New CSVs
    path_parcellation_to_etiv = "/home/sol/COVID/CSV_subjects/parcellation_etiv.csv"
    path_parcellation_DKT_to_etiv = "/home/sol/COVID/CSV_subjects/parcellation_DKT_etiv.csv"
    path_segmentation_to_etiv = "/home/sol/COVID/CSV_subjects/segmentation_etiv.csv"

    # Read CSV
    parcellation_volumes = pd.read_csv(path_parcellation_volumes)
    parcellation_DKT_volumes = pd.read_csv(path_parcellation_DKT_volumes)
    segmentation_volumes = pd.read_csv(path_segmentation_volumes)
    brainvol_data = pd.read_csv(path_brainvol)

    parcellation_volumes_to_etiv = parcellation_to_etiv(brainvol_data, parcellation_volumes)
    parcellation_DKT_volumes_to_etiv = parcellation_to_etiv(brainvol_data, parcellation_DKT_volumes)
    segmentation_volumes_to_etiv = segmentation_to_etiv(brainvol_data, segmentation_volumes)


    # Write CSV
    parcellation_volumes_to_etiv.to_csv(path_parcellation_to_etiv, index=False)
    parcellation_DKT_volumes_to_etiv.to_csv(path_parcellation_DKT_to_etiv, index=False)
    segmentation_volumes_to_etiv.to_csv(path_segmentation_to_etiv, index=False)


    #
    # # Read Brainvol to extract ETIV value
    # etiv = brainvol_data.loc[brainvol_data["subject"] == subject]["group"].values[0]
    # etiv = CSV_covid_group.loc[CSV_covid_group["subject"] == subject]["group"].values[0]

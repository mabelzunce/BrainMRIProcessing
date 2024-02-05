import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

def convert_mm3_to_cm3(row):
    # Function to convert a row's values in mm^3 to cm³
    return row / 1000


def convert_mm3_to_cm3_brainvol(df_brainvol):

    # List of column names to convert from mm3 to cm3
    columns_brainvol_to_apply = [col for col in df_brainvol.columns if col not
                                 in ['Grey Matter Normalization (notFS)', 'Data', 'Mean Thickness']]

    df_brainvol[columns_brainvol_to_apply] = df_brainvol[columns_brainvol_to_apply].apply(convert_mm3_to_cm3)

    df_brainvol["Units"] = "cm3"
    return df_brainvol

def convert_mm3_to_cm3_parcellation(df_parcellation):

    # List of column names to convert from mm3 to cm3
    columns_parcellation_to_apply = [col for col in df_parcellation.columns if
                                     col not in ['Hemisphere', 'Data']]

    df_parcellation[columns_parcellation_to_apply] = df_parcellation[columns_parcellation_to_apply].apply(convert_mm3_to_cm3)

    df_parcellation["Units"] = "cm3"
    return df_parcellation


def convert_mm3_to_cm3_segmentation(df_segmentation):

    # List of column names to convert from mm3 to cm3
    columns_segmentation_to_apply = [col for col in df_segmentation.columns if
                                     col not in ['Data']]

    df_segmentation[columns_segmentation_to_apply] = df_segmentation[columns_segmentation_to_apply].apply(convert_mm3_to_cm3)

    df_segmentation["Units"] = "cm3"
    return df_segmentation


def new_subject_to_df(df, csv_subject, type_df='brainvol', name_subject="Unknown", data=None, hemisphere="UK", thickness = False):
    """Append a row with the information of the subject to the dataframe"""

    df_subject = pd.read_csv(csv_subject)  # read csv with data of volumes

    if type_df == 'brainvol':
        new_df_subject = pd.DataFrame(data=df_subject[["Region", "Volume"]]).set_index("Region").T
        new_df_subject.rename(index={'Volume': name_subject}, inplace=True)

    elif type_df == 'parcellation':
        if not thickness:
            new_df_subject = pd.DataFrame(data=df_subject[["StructName", "GrayVol"]]).set_index("StructName").T
            new_df_subject.rename(index={'GrayVol': name_subject}, inplace=True)
            new_df_subject["Hemisphere"] = hemisphere
        else:
            new_df_subject = pd.DataFrame(data=df_subject[["StructName", "ThickAvg"]]).set_index("StructName").T
            new_df_subject.rename(index={'ThickAvg': name_subject}, inplace=True)
            new_df_subject["Hemisphere"] = hemisphere


    elif type_df == 'segmentation':
        new_df_subject = pd.DataFrame(data=df_subject[["StructName", "Volume_mm3"]]).set_index("StructName").T
        new_df_subject.rename(index={'Volume_mm3': name_subject}, inplace=True)

    new_df_subject.columns.name = None

    if data:
        new_df_subject["Data"] = data
    return pd.concat([df, new_df_subject])


def create_df(path_subjects, df_brainvol, df_parcellation, df_segmentation, df_parcellation_thick, df_parcellation_DKT
              , df_parcellation_thick_DKT, name_subjects):

    '''Create df for brainvol.csv, parcellation.csv, segmentation.csv'''

    subjects = [f for f in sorted(os.listdir(path_subjects)) if
                os.path.isdir(os.path.join(path_subjects, f))]  # list of subjects

    existent_subjects = list(df_parcellation.index)

    subjects = [subject for subject in subjects if subject not in existent_subjects]

    # read the csv of every subject and append to the dataframe
    for subject in subjects:
        file_brainvol = f'{subject}/{subject}_brainvol.csv'  # brainvol
        file_rh_aparca2009 = f'{subject}/{subject}_rh_aparca2009s.csv'  # parcellation rh
        file_lh_aparca2009 = f'{subject}/{subject}_lh_aparca2009s.csv'  # parcellation lh
        file_rh_aparcDKT= f'{subject}/{subject}_rh_aparc.csv'  # parcellation rh
        file_lh_aparDKT= f'{subject}/{subject}_lh_aparc.csv'  # parcellation lh
        file_aseg = f'{subject}/{subject}_aseg.csv'  # segmentation

        new_csv_brainvol = os.path.join(path_subjects, file_brainvol)  # path of csv brainvol file
        new_csv_rh_aparca2009 = os.path.join(path_subjects, file_rh_aparca2009)  # path of csv rh parcellation file
        new_csv_lh_aparca2009 = os.path.join(path_subjects, file_lh_aparca2009)  # path of csv lh parcellation
        new_csv_rh_aparcDKT = os.path.join(path_subjects, file_rh_aparcDKT)  # path of csv rh parcellation file
        new_csv_lh_aparcDKT = os.path.join(path_subjects, file_lh_aparDKT)  # path of csv lh parcellation file
        new_csv_aseg = os.path.join(path_subjects, file_aseg)  # path of segmentation file

        df_brainvol = new_subject_to_df(df_brainvol, new_csv_brainvol, name_subject=subject,
                                        data=name_subjects)


        df_parcellation = new_subject_to_df(df_parcellation, new_csv_rh_aparca2009, type_df="parcellation",
                                            name_subject=subject, hemisphere="rh",
                                            data=name_subjects)

        df_parcellation_thick = new_subject_to_df(df_parcellation_thick, new_csv_rh_aparca2009, type_df="parcellation",
                                            name_subject=subject, hemisphere="rh",
                                            data=name_subjects, thickness=True)

        df_parcellation = new_subject_to_df(df_parcellation, new_csv_lh_aparca2009, type_df="parcellation",
                                            name_subject=subject, hemisphere="lh",
                                            data=name_subjects)

        df_parcellation_thick = new_subject_to_df(df_parcellation_thick, new_csv_lh_aparca2009, type_df="parcellation",
                                            name_subject=subject, hemisphere="lh",
                                            data=name_subjects, thickness=True)
        
        # DKT
        df_parcellation_DKT = new_subject_to_df(df_parcellation_DKT, new_csv_rh_aparcDKT, type_df="parcellation",
                                            name_subject=subject, hemisphere="rh",
                                            data=name_subjects)

        df_parcellation_thick_DKT = new_subject_to_df(df_parcellation_thick_DKT, new_csv_rh_aparcDKT, type_df="parcellation",
                                            name_subject=subject, hemisphere="rh",
                                            data=name_subjects, thickness=True)

        df_parcellation_DKT = new_subject_to_df(df_parcellation_DKT, new_csv_lh_aparcDKT, type_df="parcellation",
                                            name_subject=subject, hemisphere="lh",
                                            data=name_subjects)

        df_parcellation_thick_DKT = new_subject_to_df(df_parcellation_thick_DKT, new_csv_lh_aparcDKT, type_df="parcellation",
                                            name_subject=subject, hemisphere="lh",
                                            data=name_subjects, thickness=True)

        df_segmentation = new_subject_to_df(df_segmentation, new_csv_aseg, type_df="segmentation",
                                            name_subject=subject,
                                            data=name_subjects)
    return df_brainvol, df_parcellation, df_segmentation, df_parcellation_thick, df_parcellation_DKT, df_parcellation_thick_DKT

if __name__ == '__main__':

    path_COVID_sanos = "/home/sol/COVID/CSV_subjects/Control"
    path_COVID_prolongado = "/home/sol/COVID/CSV_subjects/Covid Prolongado"

    # Create DF
    df_brainvol = pd.DataFrame()
    df_parcellation = pd.DataFrame()
    df_segmentation = pd.DataFrame()
    df_parcellation_thick = pd.DataFrame()
    df_parcellation_DKT = pd.DataFrame()
    df_parcellation_thick_DKT = pd.DataFrame()

    df_brainvol, df_parcellation, df_segmentation, df_parcellation_thick, df_parcellation_DKT, df_parcellation_thick_DKT = create_df(path_COVID_sanos, df_brainvol, df_parcellation,
                                                              df_segmentation,  df_parcellation_thick, df_parcellation_DKT, df_parcellation_thick_DKT ,"Grupo Control")

    df_brainvol, df_parcellation, df_segmentation,  df_parcellation_thick,  df_parcellation_DKT, df_parcellation_thick_DKT= create_df(path_COVID_prolongado, df_brainvol, df_parcellation,
                                                              df_segmentation, df_parcellation_thick,df_parcellation_DKT, df_parcellation_thick_DKT, "COVID Prolongado")

    df_brainvol = convert_mm3_to_cm3_brainvol(df_brainvol)
    df_parcellation = convert_mm3_to_cm3_parcellation(df_parcellation)
    df_parcellation_DKT =convert_mm3_to_cm3_parcellation(df_parcellation_DKT)

    df_segmentation = convert_mm3_to_cm3_segmentation(df_segmentation)


    # CSV brain volumes general
    path_brain_volumes = "/home/sol/COVID/CSV_subjects/brain_volumes.csv"
    df_brainvol.index.rename('subject', inplace=True)
    df_brainvol.to_csv(path_brain_volumes, decimal='.')

    # CSV parcellation
    path_parcellation = "/home/sol/COVID/CSV_subjects/parcellation.csv"
    df_parcellation.index.rename('subject', inplace=True)
    df_parcellation.to_csv(path_parcellation, decimal='.')

    path_parcellation_thick = "/home/sol/COVID/CSV_subjects/parcellation_thick.csv"
    df_parcellation_thick.index.rename('subject', inplace=True)
    df_parcellation_thick.to_csv(path_parcellation_thick, decimal='.')

    path_parcellation_DKT  = "/home/sol/COVID/CSV_subjects/parcellation_DKT.csv"
    df_parcellation_DKT.index.rename('subject', inplace=True)
    df_parcellation_DKT.to_csv(path_parcellation_DKT, decimal='.')

    path_parcellation_thick_DKT = "/home/sol/COVID/CSV_subjects/parcellation_thick_DKT.csv"
    df_parcellation_thick_DKT.index.rename('subject', inplace=True)
    df_parcellation_thick_DKT.to_csv(path_parcellation_thick_DKT, decimal='.')

    # CSV segmentation
    path_segmentation = "/home/sol/COVID/CSV_subjects/segmentation.csv"
    df_segmentation.index.rename('subject', inplace=True)
    df_segmentation.to_csv(path_segmentation, decimal='.')




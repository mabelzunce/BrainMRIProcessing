import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os


# Function to convert a row's values in mm^3 to cm³
def convert_mm3_to_cm3(row):
    return row / 1000

def new_subject_to_df(df, csv_subject, type_df='brainvol', name_subject="Unknown", data=None, hemisphere="UK"):
    """Append a row with the information of the subject to the dataframe"""

    df_subject = pd.read_csv(csv_subject)  # read csv with data of volumes

    if type_df == 'brainvol':
        new_df_subject = pd.DataFrame(data=df_subject[["Region", "Volume"]]).set_index("Region").T
        new_df_subject.rename(index={'Volume': name_subject}, inplace=True)

    elif type_df == 'parcellation':
        new_df_subject = pd.DataFrame(data=df_subject[["StructName", "GrayVol"]]).set_index("StructName").T
        new_df_subject.rename(index={'GrayVol': name_subject}, inplace=True)
        new_df_subject["Hemisphere"] = hemisphere

    elif type_df == 'segmentation':
        new_df_subject = pd.DataFrame(data=df_subject[["StructName", "Volume_mm3"]]).set_index("StructName").T
        new_df_subject.rename(index={'Volume_mm3': name_subject}, inplace=True)


    new_df_subject.columns.name = None

    if data:
        new_df_subject["Data"] = data
    return pd.concat([df, new_df_subject])


def create_df(path_subjects, df_brainvol, df_parcellation, df_segmentation, name_subjects):
    '''Create df for brainvol.csv, parcellation.csv, segmentation.csv'''

    subjects = [f for f in sorted(os.listdir(path_subjects)) if
                os.path.isdir(os.path.join(path_subjects, f))]  # list of subjects

    # read the csv of every subject and append to the dataframe
    for subject in subjects:
        file_brainvol = f'{subject}/{subject}_brainvol.csv'  # brainvol
        file_rh_aparca2009 = f'{subject}/{subject}_rh_aparca2009s.csv'  # parcellation rh
        file_lh_aparca2009 = f'{subject}/{subject}_lh_aparca2009s.csv'  # parcellation lh
        file_aseg = f'{subject}/{subject}_aseg.csv'  # segmentation

        new_csv_brainvol = os.path.join(path_subjects, file_brainvol)  # path of csv brainvol file
        new_csv_rh_aparca2009 = os.path.join(path_subjects, file_rh_aparca2009)  # path of csv rh parcellation file
        new_csv_lh_aparca2009 = os.path.join(path_subjects, file_lh_aparca2009)  # path of csv lh parcellation file
        new_csv_aseg = os.path.join(path_subjects, file_aseg)  # path of segmentation file

        df_brainvol = new_subject_to_df(df_brainvol, new_csv_brainvol, name_subject=subject,
                                        data=name_subjects)


        df_parcellation = new_subject_to_df(df_parcellation, new_csv_rh_aparca2009, type_df="parcellation",
                                            name_subject=subject, hemisphere="rh",
                                            data=name_subjects)

        df_parcellation = new_subject_to_df(df_parcellation, new_csv_lh_aparca2009, type_df="parcellation",
                                            name_subject=subject, hemisphere="lh",
                                            data=name_subjects)
        df_segmentation = new_subject_to_df(df_segmentation, new_csv_aseg, type_df="segmentation",
                                            name_subject=subject,
                                            data=name_subjects)
    return df_brainvol, df_parcellation, df_segmentation

if __name__ == '__main__':

    path_COVID_sanos = "/home/sol/COVID/CSV_subjects/Control"
    path_COVID_prolongado = "/home/sol/COVID/CSV_subjects/Covid Prolongado"


    # Create DF
    df_brainvol = pd.DataFrame()
    df_parcellation = pd.DataFrame()
    df_segmentation = pd.DataFrame()

    df_brainvol, df_parcellation, df_segmentation = create_df(path_COVID_sanos, df_brainvol, df_parcellation,
                                                              df_segmentation, "Grupo Control")

    df_brainvol, df_parcellation, df_segmentation = create_df(path_COVID_prolongado, df_brainvol, df_parcellation,
                                                              df_segmentation, "COVID Prolongado")


    # List of column names to convert from mm3 to cm3
    columns_brainvol_to_apply = [col for col in df_brainvol.columns if col not
                                 in ['Grey Matter Normalization (notFS)', 'Data']]
    columns_parcellation_to_apply = [col for col in df_parcellation.columns if
                                     col not in ['Hemisphere', 'Data']]
    columns_segmentation_to_apply = [col for col in df_segmentation.columns if
                                     col not in ['Data']]

    # convert mm3 to cm3
    df_brainvol[columns_brainvol_to_apply] = df_brainvol[columns_brainvol_to_apply].apply(convert_mm3_to_cm3)
    df_parcellation[columns_parcellation_to_apply] = df_parcellation[columns_parcellation_to_apply].apply(convert_mm3_to_cm3)
    df_segmentation[columns_segmentation_to_apply] = df_segmentation[columns_segmentation_to_apply].apply(convert_mm3_to_cm3)

    # df units
    df_brainvol["Units"] = "cm3"
    df_parcellation["Units"] = "cm3"
    df_segmentation["Units"] = "cm3"

    # CSV brain volumes general
    path_brain_volumes = "/home/sol/COVID/CSV_subjects/brain_volumes.csv"
    df_brainvol.to_csv(path_brain_volumes)

    # CSV parcellation
    path_parcellation = "/home/sol/COVID/CSV_subjects/parcellation.csv"
    df_parcellation.to_csv(path_parcellation)

    # CSV segmentation
    path_segmentation = "/home/sol/COVID/CSV_subjects/segmentation.csv"
    df_segmentation.to_csv(path_segmentation)


    # regiones que voy a graficar
    regions = ['Brain Segmentation Volume', "Total cortical gray matter volume", "Grey Matter Normalization (notFS)",
               ]
    x_label = ["Volumen\ntotal cerebro", "Volumen\nmateria gris", "Materia gris\nnormalizada",
               "Volumen\nventrículos"]

    # size
    n=13


    # initialize figure with 3 subplots in a row
    fig, ax = plt.subplots(1, 3, figsize=(8, 6), gridspec_kw={'width_ratios': [1.5, 1.5, 1.5]})
    fig.subplots_adjust(hspace=0.4, wspace=0.7)  # Adjust the spacing between subplots

    # Boxplot covid
    for i, region in enumerate(regions):
        sns.boxplot(data=df_brainvol, x="Data", y=region, palette="Blues", ax=ax[i], width=0.4)

        # title
        ax[i].set_title(x_label[i])

        ax[i].title.set_size(n)

        if i == 2:
            ax[i].set_ylabel("Porcentaje (%)", fontsize=n,  labelpad=0)
        else:
            ax[i].set_ylabel(r"Volumen (cm$^3$)", fontsize=n)

        ax[i].set_xticks([y for y in range(2)],
                       labels=["Grupo\nControl", "COVID\nprolongado"], fontsize=10)

        # # xlabel
        ax[i].set_xlabel("")

    # Guardar la figura
    plt.savefig('/home/sol/COVID/Graficos/OnlyCOVID/boxplot_covid.png')



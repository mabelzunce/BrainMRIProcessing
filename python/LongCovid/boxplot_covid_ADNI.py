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


def make_boxplot(df, regions, ax, fontsize=15, widths=0.7, plot_title='Brain Volumes'):
    """
    Make a boxplot of different regions in the dataframe"""

    data = [df[region] for region in regions]
    x_label = [region for region in regions]

    bp = ax.boxplot(data, widths=widths, patch_artist=True)

    # titulo del boxplot y de los ejes
    ax.set(
        title=plot_title,
        xlabel="Regions",
        ylabel="Volume (mm3)"
    )

    ax.set_xticklabels(x_label, rotation=0, fontsize=fontsize)
    # plt.subplots_adjust(bottom=0.226)

    # fill with colors
    # colors = ['pink', 'pink', 'bisque', 'aquamarine']

    # color of boxes
    # for patch, color in zip(bp['boxes'], colors):
    #     patch.set_facecolor(color)

    # # color of median
    # for median in bp['medians']:
    #     median.set_color('red')

    return bp


if __name__ == '__main__':

    path_COVID_sanos = "/home/sol/COVID/CSV_subjects/SANOS"
    path_COVID_prolongado = "/home/sol/COVID/CSV_subjects/PROLONGADO"
    path_ADNI_AD = "/home/sol/COVID/CSV_subjects/AD"
    path_ADNI_CN = "/home/sol/COVID/CSV_subjects/CN"

    # create df
    df_brainvol = pd.DataFrame()
    df_parcellation = pd.DataFrame()
    df_segmentation = pd.DataFrame()

    df_brainvol, df_parcellation, df_segmentation = create_df(path_COVID_sanos, df_brainvol, df_parcellation,
                                                              df_segmentation, "Grupo Control")

    df_brainvol, df_parcellation, df_segmentation = create_df(path_COVID_prolongado, df_brainvol, df_parcellation,
                                                              df_segmentation, "COVID Prolongado")

    df_brainvol, df_parcellation, df_segmentation = create_df(path_ADNI_AD, df_brainvol, df_parcellation,
                                                              df_segmentation, "ADNI Alzheimer")

    df_brainvol, df_parcellation, df_segmentation = create_df(path_ADNI_CN, df_brainvol, df_parcellation,
                                                              df_segmentation, "ADNI Control")
    print(df_parcellation)

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

    # regiones que voy a graficar
    regions = ['Brain Segmentation Volume', "Total cortical gray matter volume", "Grey Matter Normalization (notFS)",
               'Volume of ventricles and choroid plexus']
    x_label = ["Volumen\ntotal cerebro", "Volumen\nmateria gris", "Materia gris\nnormalizada",
               "Volumen\nventrículos"]
    # size
    n=13


    # initialize figure with 3 subplots in a row
    fig, ax = plt.subplots(1, 4, figsize=(16, 4), gridspec_kw={'width_ratios': [1.5, 1.5, 1.5, 1.5]})
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

        ax[i].set_xticks([y for y in range(4)],
                       labels=["Grupo\nControl", "Long\nCOVID", "ADNI\nAD", "ADNI\nCN"], fontsize=10)

        # # xlabel
        ax[i].set_xlabel("")

    # Guardar la figura
    plt.savefig('/home/sol/COVID/Graficos/boxplot_covid.png')



    regions_segmentation = ['Left-Cerebellum-Cortex', 'Left-Thalamus', 'Left-Caudate', 'Left-Putamen', 'Left-Pallidum',
                            'Right-Cerebellum-Cortex', 'Right-Thalamus', 'Right-Caudate', 'Right-Putamen', 'Right-Pallidum', 'Right-Hippocampus']

    titles_segmentation = ['Cerebelo izquierdo', 'Tálamo izquierdo', 'Núcleo caudado izquierdo', 'Putamen izquierdo', 'Pallidum izquierdo',
              'Cerebelo derecho', 'Tálamo derecho', 'Núcleo caudado derecho', 'Putamen derecho', 'Pallidum derecho', 'Hipocampo derecho']
   #
   # # fig2, ax2 = plt.subplots(2, 5, figsize=(14, 10))
   #
   #  # add padding between the subplots
   #  # plt.subplots_adjust(wspace=0.5)
   #
    for i, region in enumerate(regions_segmentation):
        fig2, ax2 = plt.subplots(figsize=(3,5))

        sns.boxplot(data=df_segmentation, x="Data", y=region, palette="Blues", ax=ax2, width=0.4)

        # title
        ax2.set_title(titles_segmentation[i])

        ax2.set_ylabel(r"Volumen (cm$^3$)", fontsize=n)

        ax2.set_xticks([y for y in range(4)],
                      labels=["Grupo\nControl", "COVID", "ADNI\nAD", "ADNI\nCN"], fontsize=10)
       # ax2.set_xticks(["COVID\nControl", "COVID\nProlongado", "ADNI\nAD", "ADNI\nCN"], fontsize=10)

        # # xlabel
        ax2.set_xlabel("")

        if i == 2:
            plt.subplots_adjust(left=0.214)

        # Adjust layout
        plt.tight_layout()

        # Save fig
        plt.savefig(f'/home/sol/COVID/Graficos/{titles_segmentation[i]}.png')




regions_parcellation= ['G_pariet_inf-Angular', 'G_pariet_inf-Supramar']

titles_parcellation = ['Giro angular', 'G_parietal_inferior-Supramarginal']


for i, region in enumerate(regions_parcellation):
 fig2, ax2 = plt.subplots(figsize=(2,4))

 sns.boxplot(data=df_parcellation, x="Data", y=region, palette="Blues", ax=ax2, width=0.4)

 # title
 ax2.set_title(titles_parcellation[i])

 ax2.set_ylabel(r"Volumen (cm$^3$)", fontsize=n)

 ax2.set_xticks([y for y in range(4)],
               labels=["Grupo \nControl", "COVID", "ADNI\nAD", "ADNI\nCN"], fontsize=10)

 # # xlabel
 ax2.set_xlabel("")

 if i == 2:
     plt.subplots_adjust(left=0.214)

 # Adjust layout
 plt.tight_layout()

 # Save fig
 plt.savefig(f'/home/sol/COVID/Graficos/{titles_parcellation[i]}.png')





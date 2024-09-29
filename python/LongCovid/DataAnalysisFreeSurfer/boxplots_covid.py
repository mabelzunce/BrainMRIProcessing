import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
import os
matplotlib.use('Agg')  # Use the Agg backend (non-interactive)



def boxplot_general(df_plot, regions_plot, x_label_plot, y_label_plot, size=13):

    # initialize figure with 3 subplots in a row
    fig, ax = plt.subplots(1, 4, figsize=(10, 4), gridspec_kw={'width_ratios': [1.5, 1.5, 1.5, 1.5]})
    fig.subplots_adjust(hspace=0.4, wspace=0.7)  # Adjust the spacing between subplots

    # Boxplot covid
    for i, region in enumerate(regions_plot):
        sns.boxplot(data=df_plot, x="Data", y=region, palette="Blues", ax=ax[i], width=0.4)

        # title
        ax[i].set_title(x_label_plot[i])
        ax[i].title.set_size(size)

        # y_label
        ax[i].set_ylabel(y_label_plot[i], fontsize=size, labelpad=0)
        ax[i].set_xticks([y for y in range(2)],
                         labels=["Control", "Long\nCOVID"], fontsize=10)

        # # xlabel
        ax[i].set_xlabel("")
    return fig, ax

def do_boxplot_segmentation(regions, dataset,  output_dir, titles=None, normalized=True):
    # Boxplot Segmentation (all)

    if not titles:
        titles = regions


    for i, region in enumerate(regions):
        # Data
        control_group = dataset.loc[dataset["Data"] == "Grupo Control", region]
        covid_group = dataset.loc[dataset["Data"] == "COVID Prolongado", region]

        data = [control_group, covid_group]

        fig, ax = plt.subplots(figsize=(3, 5))

        # Create a boxplot
        plt.boxplot(data, labels=["Control", "Long COVID"])

        # title
        ax.set_title(titles[i])

        if normalized:
            ax.set_ylabel("Percentage (%)", fontsize=10)
        else:
            ax.set_ylabel("Volumen(cm3)", fontsize=10)


        # Adjust layout
        plt.tight_layout()

        output_plot = os.path.join(output_dir, f'{titles[i]}.png')

        # Save fig
        plt.savefig(output_plot)


def do_boxplot_parcellation(regions, dataset,  output_dir, titles=None, normalized=True, hemisphere='lh'):

    if not titles:
        titles = [f'{region}_{hemisphere}' for region in regions]

    for i, region in enumerate(regions):
        control_group = dataset.loc[
            (dataset["Data"] == "Grupo Control") & (dataset["Hemisphere"] == hemisphere), region]
        covid_group = dataset.loc[
            (dataset["Data"] == "COVID Prolongado") & (dataset["Hemisphere"] == hemisphere), region]

        data = [control_group, covid_group]

        fig0, ax0 = plt.subplots(figsize=(3, 5))

        # Create a boxplot
        plt.boxplot(data, labels=["Control", "Long\nCOVID"])

        # title
        ax0.set_title(titles[i])

        if normalized:
            ax0.set_ylabel("Percentage (%)", fontsize=10)
        else:
            ax0.set_ylabel(r"Volume (cm$^3$)", fontsize=10)

        # Adjust layout
        plt.tight_layout()

        output_fig = os.path.join(output_dir, f'{titles[i]}.png')

        # # Save fig
        plt.savefig(output_fig)



if __name__ == '__main__':

    # Read CSV
    # general volumes
    df_brain_volumes = pd.read_csv("/home/sol/COVID/CSV_subjects/brain_volumes.csv")

    # parcellation
    df_parcellation_to_etiv = pd.read_csv("/home/sol/COVID/CSV_subjects/parcellation_etiv.csv")
    df_parcellation = pd.read_csv("/home/sol/COVID/CSV_subjects/parcellation.csv")

    # segmentation
    df_segmentation_to_etiv = pd.read_csv("/home/sol/COVID/CSV_subjects/segmentation_etiv.csv")
    df_segmentation = pd.read_csv("/home/sol/COVID/CSV_subjects/segmentation.csv")


    # Cerebellum
    right_cerebellum = df_segmentation.loc[:, 'Right-Cerebellum-Cortex']
    left_cerebellum = df_segmentation.loc[:, 'Left-Cerebellum-Cortex']
    total_cerebellum = right_cerebellum + left_cerebellum
    total_cerebellum = pd.Series(total_cerebellum, name="Cerebellum")

    df_plot = pd.concat([df_brain_volumes, total_cerebellum],axis=1)



    # regiones que voy a graficar
    regions = ['Estimated Total Intracranial Volume', 'Total gray matter volume',
               "Mean Thickness", "Cerebellum"]
    x_label = ['Total Intracranial\nVolume', 'Total gray matter volume',
               "Mean Cortical\nThickness", "Cerebellum"]
    y_label = [r"Volume (cm$^3$)", "Percentage (%)", "Thickness (mm)", r"Volume (cm$^3$)"]


    fig1, ax1 = boxplot_general(df_plot, regions, x_label, y_label, size=13)
    plt.savefig('/home/sol/COVID/Graficos/1_2024_Boxplot_COVID.png')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


def brainvol_new_subject_to_df(df, path_subject, name_subject="Unknown"):
    """Append a row with the information of the subject to the dataframe"""

    df_brainvol = pd.read_csv(path_subject)
    new_df_brainvol = pd.DataFrame(data=df_brainvol[["Region", "Volume"]]).set_index("Region").T
    new_df_brainvol.rename(index={'Volume': name_subject}, inplace=True)

    return pd.concat([df, new_df_brainvol])


def make_boxplot_brain_volumes(regions, ax, fontsize=12, widths=0.15, plot_title='Brain Volumes'):
    """Make a boxplot of different regions in brainvol.csv"""

    data = [df[region] for region in regions]
    x_label = [region for region in regions]

    bp = ax.boxplot(data, widths=widths, patch_artist=True)

    ax.set(
        title=plot_title,
        xlabel="Regions",
        ylabel="Volume (mm3)"
    )

    ax.set_xticklabels(x_label, fontsize=fontsize)

    # fill with colors
    colors = ['pink', 'lightblue','bisque', 'aquamarine']

    # color of boxes
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    # color of median
    for median in bp['medians']:
        median.set_color('red')

    return bp


if __name__ == '__main__':

    path = r"C:\Users\ecyt\Documents\Sol\CovidMRI\stats\subjects"
    subjects = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]  # list of subjects

    df = pd.DataFrame()  # empty dataframe

    for subject in subjects:
        file = f'{subject}\{subject}_brainvol.csv'
        new_path = os.path.join(path, file)  # path of csv file
        df = brainvol_new_subject_to_df(df, new_path, name_subject=subject)

    print(df)

    regions = list(df.columns)

    # boxplots
    colors = []
    fig1, ax1 = plt.subplots()
    bp1 = make_boxplot_brain_volumes(regions[0:4], ax1)

    fig2, axs2 = plt.subplots(1, 2)
    bp2 = make_boxplot_brain_volumes(regions[5:7], axs2[0], 10, widths=0.20, plot_title="Gray matter")  # LH/RH gray matter

    bp3 = make_boxplot_brain_volumes(regions[9:11], axs2[1], 10, widths=0.20, plot_title="White matter")  # LH/RH white matter

    fig3, ax3 = plt.subplots()
    bp4 = make_boxplot_brain_volumes([regions[8], regions[11]], ax3, widths=0.08, plot_title="Total gray and white matter")  # total gray/ total white

    # statistics
    metrics = df.describe()
    print(metrics)
    metrics.to_csv("metrics.csv")

    plt.show()

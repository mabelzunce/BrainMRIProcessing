import os.path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def intensity_regions_bar_chart(variables, value_variables, errors, specific_legend=None, general_legend =None, legend_position=(1.14, 1.14)
                                ,colors=['tab:blue', 'tab:orange', 'tab:green', 'tab:red'], title=None):
    """
    Creates a bar chart using the brain regions in the x-axis and the signal intensity in the y-axis.

    Parameters
    ----------
    variable : list of variables
    intensity1: list of values for the first variable
    intensity2: list of values for the second variable
    error1: list of error values for the first variable (optional)
    error2: list of error values for the second variable (optional)
    title: title of the chart
    """

    y_pos = np.arange(len(variables))

    legend_names = [f'{i}_{name}' for i, name in enumerate(variables)]

    # Plot bars for each variable
    for i, variable in enumerate(value_variables):
        error = errors[i] if errors is not None else None
        color = colors[i]
        if specific_legend:
            label = specific_legend[i]
        else:
            label=None
        plt.bar(y_pos + (i - len(value_variables) / 2) * 0.2, variable, label=label, width=0.18, yerr=error,
                color=color, capsize=5)

    if general_legend is not None:
        plt.xticks(y_pos, [f'{name}' for name in general_legend], size=14)  # X-axis label
    else:
        plt.xticks(y_pos, [str(i + 1) for i in range(len(variables))])



    plt.xlabel("", size=16)
    plt.ylabel("Value", size=16)

    plt.legend(bbox_to_anchor=legend_position, loc="upper right", title=general_legend)

    if title is not None:
        plt.title(title, size=16)
    else:
        plt.title("Bar Chart Variable", size=16)
    plt.legend(bbox_to_anchor=legend_position)
    plt.subplots_adjust(right=0.80)


def extract_values_errors_ukbiobank(df_ukbiobank, variables):

    controls_variable_1, controls_error_variable_1 = [], []
    cases_variable_1, cases_error_variable_1 = [], []

    controls_variable_2, controls_error_variable_2 = [], []
    cases_variable_2, cases_error_variable_2 = [], []

    # Brain Volume bar-chart
    for variable in variables:
        # First time
        controls_variable_1.append(float(df_ukbiobank.loc[df_ukbiobank["variable"] == variable]['controls_mean_1']))
        controls_error_variable_1.append(
            float(df_ukbiobank.loc[df_ukbiobank["variable"] == variable]['controls_std_1']))

        cases_variable_1.append(float(df_ukbiobank.loc[df_ukbiobank["variable"] == variable]['cases_mean_1']))
        cases_error_variable_1.append(float(df_ukbiobank.loc[df_ukbiobank["variable"] == variable]['cases_std_1']))

        # Second time
        controls_variable_2.append(float(df_ukbiobank.loc[df_ukbiobank["variable"] == variable]['controls_mean_2']))
        controls_error_variable_2.append(
            float(df_ukbiobank.loc[df_ukbiobank["variable"] == variable]['controls_std_2']))

        cases_variable_2.append(float(df_ukbiobank.loc[df_ukbiobank["variable"] == variable]['cases_mean_2']))
        cases_error_variable_2.append(float(df_ukbiobank.loc[df_ukbiobank["variable"] == variable]['cases_std_2']))

    values = [controls_variable_1, cases_variable_1,
                           controls_variable_1, cases_variable_2]

    errors = [controls_error_variable_1, cases_error_variable_1,
                          controls_error_variable_1, cases_error_variable_2]

    return values, errors






if __name__ == '__main__':
    ukbiobank_data = "/home/sol/COVID/UKBiobank/GroupAverages_0101_new.txt"
    ukbiobank_plots = "/home/sol/COVID/UKBiobank/Analysis"

    headers = ["controls_mean_1","controls_std_1",
               "controls_mean_2", "controls_std_2",
               "cases_mean_1", "cases_std_1",
               "cases_mean_2", "cases_std_2",
               "variable"]

    df_ukbiobank = pd.read_csv(ukbiobank_data, sep=" ", skiprows=[0])
    df_ukbiobank.columns = headers

    # All important variables
    compare_variables = ["IDP_T1_SIENAX_brain-normalised_volume", "IDP_T1_SIENAX_brain-unnormalised_volume",
                         "IDP_T1_SIENAX_grey_normalised_volume", "IDP_T1_SIENAX_grey_unnormalised_volume",
                         "aseg_global_volume_BrainSeg", "aseg_global_volume_BrainSegNotVent", "aseg_global_volume_TotalGray",
                         "aseg_global_volume_EstimatedTotalIntraCranial"]

    # Total Brain Volume variables
    brain_volume_variables = ["IDP_T1_SIENAX_brain-normalised_volume", "IDP_T1_SIENAX_brain-unnormalised_volume",
                       "aseg_global_volume_BrainSegNotVent", "aseg_global_volume_BrainSegNotVent"]


    # Gray Matter Volume variables
    grey_matter_variables = ["IDP_T1_SIENAX_grey_normalised_volume","IDP_T1_SIENAX_grey_unnormalised_volume",
                              "aseg_global_volume_TotalGray"]


    freesurfer_ratios = ["aseg_global_volume-ratio_BrainSegVol-to-eTIV",
                         "aseg_global_volume-ratio_MaskVol-to-eTIV"]


    values_brain_volume, error_brain_volume = extract_values_errors_ukbiobank(df_ukbiobank, brain_volume_variables)
    values_GM_volume, error_GM_volume = extract_values_errors_ukbiobank(df_ukbiobank, grey_matter_variables)

    name_bar_charts = ["Controls 1", "Cases 1", "Controls 2", "Cases 2"]
    name_variables_brain_volume = ["SIENAX_brain_volume\n(normalised)", "SIENAX_brain_volume\n(unnormalised)",
                       "Freesurfer BrainSegNotVent", "Freesurfer BrainSeg"]

    name_variables_grey_matter_volume = ["SIENAX_grey_matter_volume\n(normalised)", "SIENAX_grey_matter_volume\n(unnormalised)",
                       "Freesurfer Total Gray Matter"]

    # All
    plt.figure()
    intensity_regions_bar_chart(brain_volume_variables, values_brain_volume, error_brain_volume,
                                general_legend = name_variables_brain_volume, specific_legend=name_bar_charts,
                                title="FreeSurfer vs Sienax Brain Volume")
    # Set figure size
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(16, 8)  # Set the size in inches (width, height)


    plt.savefig(os.path.join(ukbiobank_plots, "all_charts_brain_volume.png"))

    plt.figure()
    # Gray Matter
    intensity_regions_bar_chart(grey_matter_variables, values_GM_volume, error_GM_volume,
                                general_legend = name_variables_grey_matter_volume, specific_legend=name_bar_charts,
                                title="FreeSurfer vs Sienax Grey Matter Volume")

    # Set figure size
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(16, 8)  # Set the size in inches (width, height)

    plt.savefig(os.path.join(ukbiobank_plots, "all_charts_gm_volume.png"))



    # Only Second point in the time
    values_brain_volume2 = values_brain_volume[2:4]
    error_brain_volume2 = error_brain_volume[2:4]

    values_GM_volume2 = values_GM_volume[2:4]
    error_GM_volume2 =error_GM_volume[2:4]

    plt.figure()
    intensity_regions_bar_chart(brain_volume_variables, values_brain_volume2, error_brain_volume2,
                                general_legend = name_variables_brain_volume, specific_legend=name_bar_charts,
                                title="FreeSurfer vs Sienax Brain Volume")
    # Set figure size
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(16, 8)  # Set the size in inches (width, height)


    plt.savefig(os.path.join(ukbiobank_plots, "second_point_brain_volume.png"))




    plt.figure()
    intensity_regions_bar_chart(grey_matter_variables, values_GM_volume2, error_GM_volume2,
                                general_legend = name_variables_grey_matter_volume, specific_legend=name_bar_charts,
                                title="FreeSurfer vs Sienax Grey Matter Volume")

    # Set figure size
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(16, 8)  # Set the size in inches (width, height)

    plt.savefig(os.path.join(ukbiobank_plots, "second_point_gm_volume.png"))

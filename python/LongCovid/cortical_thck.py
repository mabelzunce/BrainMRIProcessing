import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

# Read CSV

parcellation_thick = pd.read_csv("/home/sol/COVID/CSV_subjects/parcellation_thick.csv")

# regions
parcellation_regions = list(parcellation_thick.columns[1:-2])

parcellation_plot = []

for i, region in enumerate(parcellation_regions):

    try:
        control_group_rh = parcellation_thick.loc[
            (parcellation_thick["Data"] == "Grupo Control") & (parcellation_thick["Hemisphere"] == "rh"), region]
        covid_group_rh = parcellation_thick.loc[
            (parcellation_thick["Data"] == "COVID Prolongado") & (parcellation_thick["Hemisphere"] == "rh"), region]

        control_group_lh = parcellation_thick.loc[
            (parcellation_thick["Data"] == "Grupo Control") & (parcellation_thick["Hemisphere"] == "lh"), region]
        covid_group_lh = parcellation_thick.loc[
            (parcellation_thick["Data"] == "COVID Prolongado") & (parcellation_thick["Hemisphere"] == "lh"), region]

        # Perform the t-test for independent samples (right hemisphere)
        t_statistic, p_value = stats.ttest_ind(control_group_rh, covid_group_rh)

        # Check if the difference is statistically significant (usually with a significance level of 0.05)
        if p_value < 0.05:
            parcellation_plot.append((region, "rh"))

            # Print the result
            print("T-value:", t_statistic)
            print("P-value:", p_value)

            print(f"{region} RH : The difference between the groups is statistically significant")

        # Perform the t-test for independent samples (left hemisphere)
        t_statistic, p_value = stats.ttest_ind(control_group_lh, covid_group_lh)

        # Check if the difference is statistically significant (usually with a significance level of 0.05)
        if p_value < 0.05:
            parcellation_plot.append((region, "lh"))

            # Print the result
            print("T-value:", t_statistic)
            print("P-value:", p_value)

            print(f"{region} LH : The difference between the groups is statistically significant")

    except:
        pass

# Boxplot Parcellation (only signficative)
for region, hemisphere in parcellation_plot:
    if hemisphere == "rh":
        control_group = parcellation_thick.loc[
            (parcellation_thick["Data"] == "Grupo Control") & (parcellation_thick["Hemisphere"] == "rh"), region]
        covid_group = parcellation_thick.loc[
            (parcellation_thick["Data"] == "COVID Prolongado") & (parcellation_thick["Hemisphere"] == "rh"), region]
    else:
        control_group = parcellation_thick.loc[
            (parcellation_thick["Data"] == "Grupo Control") & (parcellation_thick["Hemisphere"] == "lh"), region]
        covid_group = parcellation_thick.loc[
            (parcellation_thick["Data"] == "COVID Prolongado") & (parcellation_thick["Hemisphere"] == "lh"), region]

    data = [control_group, covid_group]

    fig, ax = plt.subplots(figsize=(3, 5))

    # Create a boxplot
    plt.boxplot(data, labels=["Control", "Long\nCOVID"])

    # title
    # ax.set_title(f'{name_regions[i]}')
    ax.set_title(f'{region}_{hemisphere}')

    ax.set_ylabel(r"Average Thickness (mm)", fontsize=10)

    # Adjust layout
    plt.tight_layout()

#    i += 1

    # # Save fig
    plt.savefig(f'/home/sol/COVID/Graficos/Parcellation_Thick/Significative/{region}_{hemisphere}.png')



# Boxplot Parcellation (all)
for region in parcellation_regions:
    # Right Hemisphere
    hemisphere = "rh"
    control_group = parcellation_thick.loc[
        (parcellation_thick["Data"] == "Grupo Control") & (parcellation_thick["Hemisphere"] == "rh"), region]
    covid_group = parcellation_thick.loc[
        (parcellation_thick["Data"] == "COVID Prolongado") & (parcellation_thick["Hemisphere"] == "rh"), region]

    data_rh = [control_group, covid_group]


    fig0, ax0 = plt.subplots(figsize=(3, 5))

    # Create a boxplot
    plt.boxplot(data_rh, labels=["Control", "Long\nCOVID"])

    # title
    # ax.set_title(f'{name_regions[i]}')
    ax0.set_title(f'{region}_{hemisphere}')

    ax0.set_ylabel(r"Average Thickness (mm)", fontsize=10)

    # Adjust layout
    plt.tight_layout()

    #    i += 1

    # # Save fig
    plt.savefig(f'/home/sol/COVID/Graficos/Parcellation_Thick/All/{region}_RH.png')

    # Left Hemisphere
    hemisphere = "lh"
    control_group = parcellation_thick.loc[
        (parcellation_thick["Data"] == "Grupo Control") & (parcellation_thick["Hemisphere"] == "lh"), region]
    covid_group = parcellation_thick.loc[
        (parcellation_thick["Data"] == "COVID Prolongado") & (parcellation_thick["Hemisphere"] == "lh"), region]

    data_lh = [control_group, covid_group]

    fig1, ax1 = plt.subplots(figsize=(3, 5))

    # Create a boxplot
    plt.boxplot(data_lh, labels=["Control", "Long\nCOVID"])

    # title
    # ax.set_title(f'{name_regions[i]}')
    ax1.set_title(f'{region}_{hemisphere}')

    ax1.set_ylabel(r"Average Thickness (mm)", fontsize=10)

    # Adjust layout
    plt.tight_layout()

    # # Save fig
    plt.savefig(f'/home/sol/COVID/Graficos/Parcellation_Thick/All/{region}_LH.png')


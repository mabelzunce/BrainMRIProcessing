
import statsmodels.api as sm
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt



# Read CSV
# general volumes
brain_volumes = pd.read_csv("/home/sol/COVID/CSV_subjects/brain_volumes.csv")

# parcellation
parcellation = pd.read_csv("/home/sol/COVID/CSV_subjects/parcellation.csv")

# segmentation
segmentation = pd.read_csv("/home/sol/COVID/CSV_subjects/segmentation.csv")

# regions
segmentation_regions = list(segmentation.columns[1:])
parcellation_regions = list(parcellation.columns[1:])
brain_volumes_values = list(brain_volumes.columns[1:])

# Q Q plot segmentation
for i, region in enumerate(segmentation_regions):
    try:
        control_group = segmentation.loc[segmentation["Data"] == "Grupo Control", region]
        covid_group = segmentation.loc[segmentation["Data"] == "COVID Prolongado", region]

        # Create Q-Q plot for variable 1
        qqplot1 = sm.qqplot(control_group, line='s', markerfacecolor='r', markeredgecolor='r')
        plt.title(f"Q-Q Plot {region} Control Group")

        # Save the Q-Q plot for variable 1 as an image file
        plt.savefig(f'/home/sol/COVID/Graficos/QQplot/Segmentation/qqplot_{region}_controlgroup.png')

        # Clear the current figure to create a new plot
        plt.clf()

        # Create Q-Q plot for variable 2
        qqplot2 = sm.qqplot(covid_group, line='s', markerfacecolor='g', markeredgecolor='g')
        plt.title(f"Q-Q Plot {region} COVID Group")

        # Save the Q-Q plot for variable 2 as an image file
        plt.savefig(f'/home/sol/COVID/Graficos/QQplot/qqplot_{region}_COVIDgroup.png')
    except TypeError:
        pass

# Q Q plot parcellation
for i, region in enumerate(parcellation_regions):
    try:
        control_group = parcellation.loc[parcellation["Data"] == "Grupo Control", region]
        covid_group = parcellation.loc[parcellation["Data"] == "COVID Prolongado", region]

        # Create Q-Q plot for variable 1
        qqplot1 = sm.qqplot(control_group, line='s', markerfacecolor='r', markeredgecolor='r')
        plt.title(f"Q-Q Plot {region} Control Group")

        # Save the Q-Q plot for variable 1 as an image file
        plt.savefig(f'/home/sol/COVID/Graficos/QQplot/Parcellation/qqplot_{region}_controlgroup.png')

        # Clear the current figure to create a new plot
        plt.clf()

        # Create Q-Q plot for variable 2
        qqplot2 = sm.qqplot(covid_group, line='s', markerfacecolor='g', markeredgecolor='g')
        plt.title(f"Q-Q Plot {region} COVID Group")

        # Save the Q-Q plot for variable 2 as an image file
        plt.savefig(f'/home/sol/COVID/Graficos/QQplot/Parcellation/qqplot_{region}_COVIDgroup.png')

    except TypeError:
        pass
    except AttributeError:
        pass

# Q Q plot brainvol

# Q Q plot parcellation
for i, region in enumerate(brain_volumes_values):
    try:
        control_group = parcellation.loc[parcellation["Data"] == "Grupo Control", region]
        covid_group = parcellation.loc[parcellation["Data"] == "COVID Prolongado", region]

        # Create Q-Q plot for variable 1
        qqplot1 = sm.qqplot(control_group, line='s', markerfacecolor='r', markeredgecolor='r')
        plt.title(f"Q-Q Plot {region} Control Group")

        # Save the Q-Q plot for variable 1 as an image file
        plt.savefig(f'/home/sol/COVID/Graficos/QQplot/BrainVol/qqplot_{region}_controlgroup.png')

        # Clear the current figure to create a new plot
        plt.clf()

        # Create Q-Q plot for variable 2
        qqplot2 = sm.qqplot(covid_group, line='s', markerfacecolor='g', markeredgecolor='g')
        plt.title(f"Q-Q Plot {region} COVID Group")

        # Save the Q-Q plot for variable 2 as an image file
        plt.savefig(f'/home/sol/COVID/Graficos/QQplot/BrainVol/qqplot_{region}_COVIDgroup.png')

    except TypeError:
        pass
    except AttributeError:
        pass




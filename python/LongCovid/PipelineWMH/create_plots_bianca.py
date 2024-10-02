import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.plotting import boxplot

colors = ['#fa5563', '#5988ff']
bianca_data = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/bianca_data_thr_0_7.0_cluster_5.csv"
boxplot_bianca = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/PlotsBianca/boxplot_bianca_0_7_thr_0_7_cluster_5.png"
histogram_bianca = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/PlotsBianca/histogram_bianca_0_7_thr_0_7_cluster_5.png"

# Read df
df_bianca = pd.read_csv(bianca_data)

# Create the figure
plt.figure(figsize=(8,8))

# Create the boxplots
plt.subplot(1, 2, 1)
sns.boxplot(
    x='Grupo',
    y='Number of clusters of WMH',
    data=df_bianca,  # Your dataframe
    palette=colors   # Apply custom colors
)

sns.swarmplot(x='Grupo', y='Number of clusters of WMH', data=df_bianca, color='black', alpha=0.5, size=8)


plt.title('Clusters of WMH', fontsize=12)
plt.xlabel('Group', fontsize=10)
plt.ylabel('Number of clusters', fontsize=10)
plt.grid()

plt.subplot(1, 2, 2)

sns.boxplot(
    x='Grupo',
    y='Total volume of clusters',
    data=df_bianca,  # Your dataframe
    palette=colors   # Apply custom colors
)

sns.swarmplot(x='Grupo', y='Total volume of clusters', data=df_bianca, color='black', alpha=0.5, size=8)


plt.title('Total volume of clusters', fontsize=12)
plt.xlabel('Group', fontsize=10)
plt.ylabel(r'Volume (mm$^3$)', fontsize=10)
plt.grid()

# Adjust layout and space between subplots
plt.subplots_adjust(wspace=0.4)  # Adjust horizontal space between plots
plt.savefig(boxplot_bianca)



df_control = df_bianca[df_bianca['Grupo'] == 'CONTROL']
df_covid = df_bianca[df_bianca['Grupo'] == 'COVID']

plt.figure(figsize=(10, 6))

plt.subplot(1, 2, 1)

# Define bin edges with a width of 5
bin_edges = np.arange(df_control['Number of clusters of WMH'].min(),
                      df_control['Number of clusters of WMH'].max() + 2, 2)

# Histogram for Number of lesions
sns.histplot(
    df_control['Number of clusters of WMH'],
    bins=bin_edges,  # Number of bins
    color=colors[0],  # Choose a color for the histogram
    kde=True,# Kernel Density Estimate
    stat='density'
)

plt.title('Distribution of Number of WMH clusters(CONTROL)', fontsize=12)
plt.xlabel('Number of WMH clusters', fontsize=10)
plt.ylabel('Frequency per Bin (n)', fontsize=10)
plt.grid()

# Set limits for x and y axes
plt.xlim(0, 100)  # Adjust these values based on your data
#plt.ylim(0, 6)  # Adjust these values based on your data

plt.subplot(1, 2, 2)



# Histogram for Number of lesions
sns.histplot(
    df_covid['Number of clusters of WMH'],
    bins=bin_edges,  # Use calculated bin edges for a width of 5
    #bins=20,  # Number of bins
    color=colors[1],  # Choose a color for the histogram
    kde=True,  # Kernel Density Estimate
stat='density'
)

plt.title('Distribution of Number of WMH clusters(COVID)', fontsize=12)
plt.xlabel('Number of clusters of WMH', fontsize=10)
plt.ylabel('Frequency per Bin (n)', fontsize=10)
plt.grid()

# Set limits for x and y axes
plt.xlim(0, 100)  # Adjust these values based on your data
#plt.ylim(0, 30)  # Adjust these values based on your data



# Show the histogram plot
plt.tight_layout()
plt.savefig(histogram_bianca)
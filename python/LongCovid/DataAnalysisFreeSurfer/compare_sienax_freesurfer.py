import os
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import ols

def exclude_subject_by_id(df, subjects_exclude):

    for subject in subjects_exclude:
        df = df[df['ID'] != subject]

    return df

def violin_plot(df_plot, region, title, y_label_plot, size=13):

    # initialize figure with 3 subplots in a row
    fig, ax = plt.subplots(figsize=(6, 5))

    sns.violinplot(data=df_plot, x="Grupo", y=region, palette='Blues')

    # title
    ax.set_title(title)
    ax.title.set_size(size)

    # y_label
    ax.set_ylabel(y_label_plot)


    return fig, ax


path_sienax_csv_2 = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Datos COVID - Disco/DataAnalysis/SienaxResults2.csv'
path_CSV_covid_group = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Datos COVID - Disco/DataAnalysis/RespuestasCuestionarioEvaluacionCognitivaResonancia.xlsx" # Subjects groups

# Freesurfer CSV
path_freesurfer_brainvol = "/home/sol/COVID/CSV_subjects/brain_volumes.csv"
path_freesurfer_aseg = "/home/sol/COVID/CSV_subjects/segmentation.csv"
path_freesurfer_aparc = "/home/sol/COVID/CSV_subjects/parcellation.csv"
path_freesurfer_cortical_thickness = "/home/sol/COVID/CSV_subjects/parcellation_thick.csv"


path_plot = "/home/sol/COVID/New_Plot"

if not os.path.exists(path_plot):
    os.mkdir(path_plot)


subjects_to_exclude = ['CP0011', 'CP0015', 'CP0144', 'CP0035', 'CP0106']
variables_names_sienax = ['GreyMatterVolume', 'WhiteMatterVolume', 'BrainVolume',
                   'GreyMatterUnNormVolume', 'WhiteMatterUnNNormVolume', 'BrainUnNNormVolume']

variables_names_freesurfer = ['Brain Segmentation Volume', 'Brain Segmentation Volume Without Ventricles',
                              'Supratentorial volume', 'Supratentorial volume NotVent', 'Subcortical gray matter volume', 'Left hemisphere cortical gray matter volume', 'Right hemisphere cortical gray matter volume', 'Total cortical gray matter volume', 'Total gray matter volume', 'Left hemisphere cerebral white matter volume', 'Right hemisphere cerebral white matter volume',
                              'Total cerebral white matter volume', 'Mask Volume', 'Supratentorial volume voxel count', 'Brain Segmentation Volume Without Ventricles from Surf', 'Volume of ventricles and choroid plexus', 'Estimated Total Intracranial Volume', 'Mean Thickness']

variables_important_freesurfer = ['Total gray matter volume', 'Total cerebral white matter volume', 'Brain Segmentation Volume' ]
variables_ancova_freesurfer = ['Total gray matter volume', 'Total cerebral white matter volume', 'Brain Segmentation Volume' ]
variables_ancova_sienax = 'GreyMatterUnNormVolume', 'WhiteMatterUnNNormVolume', 'BrainUnNNormVolume'
variables_important_sienax = ['GreyMatterVolume', 'WhiteMatterVolume', 'BrainVolume']

# Read and merge df (groups and sienax results)
df_groups = pd.read_excel(path_CSV_covid_group)
df_sienax = pd.read_csv(path_sienax_csv_2)

# Freesurfer
df_brainvol_fs = pd.read_csv(path_freesurfer_brainvol)
df_segmentation = pd.read_csv(path_freesurfer_aseg)
df_aparc = pd.read_csv(path_freesurfer_aparc)
df_cortical_thickness = pd.read_csv(path_freesurfer_cortical_thickness)

df_merged_sienax = pd.merge(df_groups, df_sienax, left_on='ID', right_on='SubjectNames')
df_merged_sienax = exclude_subject_by_id(df_merged_sienax, subjects_to_exclude)

df_merged_freesurfer_brainvol = pd.merge(df_groups, df_brainvol_fs, left_on='ID', right_on='subject')
df_merged_freesurfer_brainvol = exclude_subject_by_id(df_merged_freesurfer_brainvol, subjects_to_exclude)

df_merged_freesurfer_segmentation = pd.merge(df_groups, df_segmentation, left_on='ID', right_on='subject')
df_merged_freesurfer_segmentation = exclude_subject_by_id(df_merged_freesurfer_segmentation, subjects_to_exclude)

variables_fs_segmentation = list(df_segmentation.columns[1:-2])
variables_fs_parcellation = list(df_aparc.columns[1:-3])
# variables_fs_cortical_thickness = df_aparc.columns[1:-3]

# Sienax
for variable in variables_names_sienax:
    control_group = df_merged_sienax[df_merged_sienax['Grupo'] == 'CONTROL'][variable]
    control_group_ID = df_merged_sienax[df_merged_freesurfer_brainvol['Grupo'] == 'CONTROL']['ID']
    control_group_etiv = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'CONTROL'][
        'Estimated Total Intracranial Volume']
    control_group_age = df_merged_sienax[df_merged_sienax['Grupo'] == 'CONTROL']['Edad']

    covid_group = df_merged_sienax[df_merged_sienax['Grupo'] == 'COVID'][variable]
    covid_group_etiv = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'COVID'][
        'Estimated Total Intracranial Volume']
    covid_group_ID = df_merged_sienax[df_merged_freesurfer_brainvol['Grupo'] == 'COVID']['ID']
    covid_group_age = df_merged_sienax[df_merged_sienax['Grupo'] == 'COVID']['Edad']

    # Do Kruskal Wallis
    result = stats.kruskal(control_group, covid_group)

    # Violin plot
    # Create Dataframe
    df_control = pd.DataFrame({
        'Grupo': 'CONTROL',
        'ID': control_group_ID,
        'Age': control_group_age,
        'Value': control_group,
        'Variable': variable,
        'etiv': control_group_etiv,
    })

    df_covid = pd.DataFrame({
        'Grupo': 'COVID',
        'ID': covid_group_ID,
        'Age': covid_group_age,
        'Value': covid_group,
        'Variable': variable,
        'etiv': covid_group_etiv,

    })

    # Combine dataframe
    df_combined = pd.concat([df_control, df_covid], ignore_index=True)

    if variable in variables_important_sienax:
        fig, ax = violin_plot(df_combined, "Value", variable, "Variable", size=13)
        plt.savefig(os.path.join(path_plot, f'{variable}_normalized_sienax.png'))
    else:
        fig, ax = violin_plot(df_combined, "Value", variable, "Variable", size=13)
        plt.savefig(os.path.join(path_plot, f'{variable}_unnormalized_sienax.png'))

    if variable in variables_important_sienax:
        # ANCOVA
        model = ols(f'Value ~ Grupo + Age', data=df_combined).fit()

        # Imprimir resultados del modelo
        print(f'Sienax {variable}:', model.summary())

try:
    for variable in variables_names_freesurfer:
        control_group = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'CONTROL'][variable]
        control_group_etiv = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'CONTROL']['Estimated Total Intracranial Volume']
        control_group_age = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'CONTROL']['Edad']
        control_group_ID = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'CONTROL']['ID']
        control_group_norm_var = control_group /control_group_etiv

        covid_group = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'COVID'][variable]
        covid_group_age = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'COVID']['Edad']
        covid_group_etiv = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'COVID']['Estimated Total Intracranial Volume']
        covid_group_ID = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'COVID']['ID']
        covid_group_norm_var = covid_group / covid_group_etiv

        # Do Kruskal Wallis
        result = stats.kruskal(control_group_norm_var, covid_group_norm_var)


        # Violin plot
        # Create Dataframe
        df_control = pd.DataFrame({
            'Grupo': 'CONTROL',
            'etiv': control_group_etiv,
            'Age': control_group_age,
            'ID': control_group_ID,
            'Variable': variable,
            'Value': control_group,
            'Normalized_Variable': control_group_norm_var
        })

        df_covid = pd.DataFrame({
            'Grupo': 'COVID',
            'ID': covid_group_ID,
            'etiv': covid_group_etiv,
            'Age': covid_group_age,
            'Variable': variable,
            'Value': covid_group,
            'Normalized_Variable': covid_group_norm_var
        })

        # Combine dataframe
        df_combined = pd.concat([df_control, df_covid], ignore_index=True)


        if variable in variables_important_freesurfer:
            fig, ax = violin_plot(df_combined, "Normalized_Variable", variable, "Normalized_Variable", size=13)
            plt.savefig(os.path.join(path_plot, f'{variable}_normalized_freesurfer.png'))

            fig, ax = violin_plot(df_combined, "Value", variable, "Value", size=13)
            plt.savefig(os.path.join(path_plot, f'{variable}_unnormalized_freesurfer.png'))

        if variable in variables_ancova_freesurfer:
            # ANCOVA
            model = ols(f'Value ~ Grupo + Age + etiv', data=df_combined).fit()


            # Imprimir resultados del modelo
            print(f'Freesurfer {variable}:', model.summary())


except ValueError:
    pass
    #print(variable)


# Segmentation
p_values_segmentation = {}

for variable in variables_fs_segmentation:

    control_group_value = df_merged_freesurfer_segmentation[df_merged_freesurfer_segmentation['Grupo'] == 'CONTROL'][variable]
    control_group_etiv = df_merged_freesurfer_brainvol[df_merged_freesurfer_brainvol['Grupo'] == 'CONTROL'][
        'Estimated Total Intracranial Volume']
    control_group_age = df_merged_freesurfer_segmentation[df_merged_freesurfer_segmentation['Grupo'] == 'CONTROL']['Edad']
    control_group_ID = df_merged_freesurfer_segmentation[df_merged_freesurfer_segmentation['Grupo'] == 'CONTROL']['ID']
    control_group_normalized_value = control_group_value / control_group_etiv

    covid_group_value = df_merged_freesurfer_segmentation[df_merged_freesurfer_segmentation['Grupo'] == 'COVID'][variable]
    covid_group_etiv = df_merged_freesurfer_brainvol[df_merged_freesurfer_segmentation['Grupo'] == 'COVID'][
        'Estimated Total Intracranial Volume']
    covid_group_age = df_merged_freesurfer_segmentation[df_merged_freesurfer_segmentation['Grupo'] == 'COVID']['Edad']
    covid_group_ID = df_merged_freesurfer_segmentation[df_merged_freesurfer_segmentation['Grupo'] == 'COVID']['ID']
    covid_group_normalized_value = covid_group_value / covid_group_etiv

    # Create Dataframe
    df_control = pd.DataFrame({
        'Group': 'CONTROL',
        'Age': control_group_age,
        'ID': control_group_ID,
        'Variable': variable,
        'Value': control_group,
        'Normalized_value': control_group_normalized_value
    })

    df_covid = pd.DataFrame({
        'Group': 'COVID',
        'ID': covid_group_ID,
        'Age': covid_group_age,
        'Variable': variable,
        'Value': covid_group,
        'Normalized_value': covid_group_normalized_value
    })

    # Combine dataframe
    df_combined = pd.concat([df_control, df_covid], ignore_index=True)

    # ANCOVA
    model = ols(f'Normalized_value ~ Group + Age', data=df_combined).fit()

    # Get the summary of the model
    summary = model.summary()

    # Extract p-values for each variable
   # p_values = summary.tables[1]['P>|t|']

    # Imprimir resultados del modelo
    print(f'Freesurfer Segmentation {variable}:',summary.tables[1]['P>|t|'])
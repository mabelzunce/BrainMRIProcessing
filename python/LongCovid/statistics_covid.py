import pandas as pd 
import scipy.stats as stats
import os 

def do_test_t(covid_group, control_group, name_variable=None, verbose=False):
    '''Do Test T'''

    # Perform the independent t-test
    t_statistic, p_value = stats.ttest_ind(covid_group, control_group)

    # Interpret the results
    alpha = 0.05
    if verbose:
        if p_value < alpha:
            print(f"Test T: {name_variable}  There is a significant difference between the two groups.")

    
    return t_statistic, p_value


def do_mann_whitneyu(covid_group, control_group, name_variable=None, verbose=False):
    '''Do Mann Whitney U test'''

    # Perform the independent t-test
    t_statistic, p_value = stats.mannwhitneyu(covid_group, control_group)

    # Interpret the results
    alpha = 0.05
    if verbose:
        if p_value < alpha:
            print(f"Mann-Whitney U test {name_variable}  There is a significant difference between the two groups.")

    
    return t_statistic, p_value

def write_a_row_with_statistics(df_statistics, covid_values, control_values, name_variable, name_atlas, hemisphere='', verbose=False):
    '''Append to df a new row with statistic of variable "name_variable" '''
    control_values_mean = control_values.mean()
    control_values_std = control_values.std()

    covid_values_mean = covid_values.mean()
    covid_values_std = covid_values.std()

    t_statistic, p_value= do_test_t(covid_values, control_values, f'{name_variable} {hemisphere}', verbose)
    statistic, p_value_mw = do_mann_whitneyu(covid_values, control_values, f'{name_variable} {hemisphere}', verbose)
    

    # Create new row with rh values
    df_stat = pd.DataFrame({ 'variable': [f'{name_atlas}_{name_variable}_{hemisphere}'], 
                            'CONTROL mean': [control_values_mean],
                            'CONTROL std': [control_values_std],
                            'COVID mean': [covid_values_mean], 
                            'COVID std': [covid_values_std], 
                            'p value(test T)': [p_value],
                            'p value (test Mann Whitneyu)': [p_value_mw]
                                })
    
    df_statistics = pd.concat([df_statistics, df_stat])
    
    return df_statistics
    
    



def segmentation_statistics(df_segmentation, name_atlas='aseg'): 
    '''"A DataFrame with statistics for each region of the segmentation atlas'''

    df_statistics = pd.DataFrame()    
      
    # Columns with variables
    df_segmentation_variables = list(df_segmentation.columns)[1:-2]
    df_segmentation_volumes = []
    df_segmentation_global = []


    # Divide into global and volumes
    for variable in df_segmentation_variables:
       
        
        if variable.split('-')[0] == 'Right' or variable.split('-')[0] == 'Left':
            if variable.split('-')[-1] == 'hypointensities':
                continue
            name_variable = '-'.join(variable.split('-')[1:])
            if not name_variable in df_segmentation_volumes:
                df_segmentation_volumes.append(name_variable)
            
        else:
            df_segmentation_global.append(variable)
        
    
    # Statistics volumes both hemisphere (lh, rh)
    for variable in df_segmentation_volumes:
        name_atlas_volume  = f'{name_atlas}-volume'
        
        # Right Hemisphere
        control_values_rh = df_segmentation.loc[(df_segmentation['Data'] == 'Grupo Control'), f'Right-{variable}']
        control_values_rh = control_values_rh.reset_index(drop=True)

        covid_values_rh = df_segmentation.loc[(df_segmentation['Data'] == 'COVID Prolongado'), f'Right-{variable}']
        covid_values_rh = covid_values_rh.reset_index(drop=True)

        control_values_lh = df_segmentation.loc[(df_segmentation['Data'] == 'Grupo Control'), f'Left-{variable}']
        control_values_lh = control_values_lh.reset_index(drop=True)

        covid_values_lh = df_segmentation.loc[(df_segmentation['Data'] == 'COVID Prolongado'), f'Left-{variable}']
        covid_values_lh = covid_values_lh.reset_index(drop=True)

                
        # Mean Hemisphere 
        control_values = (control_values_lh + control_values_rh) / 2
        covid_values = (covid_values_lh + covid_values_rh) / 2

        df_statistics = write_a_row_with_statistics(df_statistics, covid_values_rh, control_values_rh, hemisphere= "rh", 
                                                    name_variable=variable, name_atlas=name_atlas_volume)

        df_statistics = write_a_row_with_statistics(df_statistics, covid_values_lh, control_values_lh, hemisphere= "lh", 
                                                    name_variable=variable, name_atlas=name_atlas_volume)
        
        df_statistics = write_a_row_with_statistics(df_statistics, covid_values, control_values, hemisphere= "mean", 
                                                    name_variable=variable, name_atlas=name_atlas_volume)
    
    # Statistics global variables
    for variable in df_segmentation_global:
        
        name_atlas_global  = f'{name_atlas}-global'
        control_values = df_segmentation.loc[(df_segmentation['Data'] == 'Grupo Control'), variable]
        control_values = control_values.reset_index(drop=True)

        covid_values = df_segmentation.loc[(df_segmentation['Data'] == 'COVID Prolongado'), variable]
        
        covid_values= covid_values.reset_index(drop=True)
        
        df_statistics = write_a_row_with_statistics(df_statistics, covid_values, control_values, 
                                                    name_variable=variable, name_atlas=name_atlas_global, verbose=True)

    
    return df_statistics
    
   





def parcellation_statistics(df_parcellation, name_atlas='aparc'): 
    '''"A DataFrame with statistics for each region of the parcellation atlas'''

    df_statistics = pd.DataFrame(df_parcellation)
    # Columns with variables
    df_parcellation_variables = list(df_parcellation.columns)[1:-3]
    
    df_statistics = pd.DataFrame()

    for variable in df_parcellation_variables:

        # Right Hemisphere
        control_values_rh = df_parcellation.loc[(df_parcellation['Data'] == 'Grupo Control')
                                                & (df_parcellation['Hemisphere'] == 'rh'), variable]
        control_values_rh = control_values_rh.reset_index(drop=True)

        covid_values_rh = df_parcellation.loc[(df_parcellation['Data'] == 'COVID Prolongado')
                                                & (df_parcellation['Hemisphere'] == 'rh'), variable]
        covid_values_rh = covid_values_rh.reset_index(drop=True)
        
        # Left Hemisphere
        control_values_lh = df_parcellation.loc[(df_parcellation['Data'] == 'Grupo Control')
                                                & (df_parcellation['Hemisphere'] == 'lh'), variable]
        control_values_lh = control_values_lh.reset_index(drop=True)

        covid_values_lh = df_parcellation.loc[(df_parcellation['Data'] == 'COVID Prolongado')
                                                & (df_parcellation['Hemisphere'] == 'lh'), variable]
        covid_values_lh = covid_values_lh.reset_index(drop=True)
        
        # Mean Hemisphere 
        control_values = (control_values_lh + control_values_rh) / 2
        covid_values = (covid_values_lh + covid_values_rh) / 2
        
        # Write rows
        df_statistics = write_a_row_with_statistics(df_statistics, covid_values_rh, control_values_rh, hemisphere= "rh", 
                                                    name_variable=variable, name_atlas=name_atlas)

        df_statistics = write_a_row_with_statistics(df_statistics, covid_values_lh, control_values_lh, hemisphere= "lh", 
                                                    name_variable=variable, name_atlas=name_atlas)
        
        df_statistics = write_a_row_with_statistics(df_statistics, covid_values, control_values, hemisphere= "mean", 
                                                    name_variable=variable, name_atlas=name_atlas)
    return df_statistics


path_csv = "data"
# Volumes


path_brain_volumes = os.path.join(path_csv, 'brain_volumes.csv')

# General
# Parcellation
path_parcellation = os.path.join(path_csv, "parcellation.csv")
path_parcellation_DKT = os.path.join(path_csv, "parcellation_DKT.csv")

# Segmentation
path_segmentation = os.path.join(path_csv, "segmentation.csv")

# Parcellation normalized to etiv 
path_parcellation_etiv = os.path.join(path_csv, "parcellation_etiv.csv")
path_parcellation_DKT_etiv = os.path.join(path_csv, "parcellation_DKT_etiv.csv")

# Segmentation normalized to etiv
path_segmentation_etiv = os.path.join(path_csv, "segmentation_etiv.csv")

# Cortical Thickness
path_parcellation_thick = os.path.join(path_csv, "parcellation_thick.csv")
path_parcellation_DKT_thick = os.path.join(path_csv, "parcellation_thick_DKT.csv")

# Read dataframes 
df_brainvol = pd.read_csv(path_brain_volumes)

# Parcellation
# Volumes
df_parcellation = pd.read_csv(path_parcellation)
df_parcellation_DKT = pd.read_csv(path_parcellation_DKT)

# Volumes to eTIV
df_parcellation_etiv = pd.read_csv(path_parcellation_etiv)
df_parcellation_DKT_etiv = pd.read_csv(path_parcellation_DKT_etiv)

# Cortical Thickness
df_parcellation_thick = pd.read_csv(path_parcellation_thick)
df_parcellation_DKT_thick = pd.read_csv(path_parcellation_DKT_thick)

# Segmentation
df_segmentation = pd.read_csv(path_segmentation)
df_segmentation_etiv = pd.read_csv(path_segmentation_etiv)


# Parcellation
# Columns with variables
df_statistics = pd.DataFrame()

# Segmentation Statistics
df_statistics_segmentation=segmentation_statistics(df_segmentation)
df_statistics = pd.concat([df_statistics,df_statistics_segmentation])

# Brain Volumes
df_statistics_brainvol=segmentation_statistics(df_brainvol)
df_statistics = pd.concat([df_statistics,df_statistics_brainvol])

# Parcellation Statistics
# Volumes
df_statistics_parcellation_volumes = parcellation_statistics(df_parcellation, name_atlas="aparc-a2009s-volume")
df_statistics = pd.concat([df_statistics, df_statistics_parcellation_volumes])


df_statistics_parcellation_DKT_volumes = parcellation_statistics(df_parcellation_DKT, name_atlas="aparc-Desikan-volume")
df_statistics = pd.concat([df_statistics, df_statistics_parcellation_DKT_volumes])

# Cortical Thickness
df_statistics_parcellation_thick = parcellation_statistics(df_parcellation_thick, name_atlas="aparc-a2009s-thickness")
df_statistics = pd.concat([df_statistics, df_statistics_parcellation_thick])

df_statistics_parcellation_DKT_thick = parcellation_statistics(df_parcellation_DKT_thick, name_atlas="aparc-Desikan-thickness")
df_statistics = pd.concat([df_statistics,df_statistics_parcellation_DKT_thick])

df_statistics.to_csv('statistics.csv')


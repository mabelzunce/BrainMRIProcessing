import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import numpy as np
import os,re
from scipy.stats import kruskal


participants_csv= '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/VoluntariosProyectoCovidProlongado.csv'
variables_participants_csv ='/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/RespuestasCuestionariosCovid_24_10.csv'

# White matter
wmh_total_csv = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/WMH/bianca_data_thr_0_7.0_cluster_5.csv"
wmh_deep_csv = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/WMH/bianca_data_thr_0_7.0_cluster_5_deep.csv"

folder_path_bids_1 = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/BIDS/derivatives/ExploreASL/Population/Stats'
folder_path_bids_2 = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/BIDS_2/derivatives/ExploreASL/Population/Stats'
output_data = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/Final_data_csv'

output_plots = os.path.join(output_data, 'Plots')

output_plots_atlas = ['DKT', 'MNI_cortical', 'Hammers', 'Harvard-Oxford_cortical', 'Harvard-Oxford_subcortical']
import pandas as pd
import numpy as np
import statsmodels.api as sm

def fit_linear_model(data, dependent_var, independent_vars):
    """
    Ajusta un modelo de regresión lineal a los datos con un ajuste opcional por edad.
    
    Parámetros:
    - data: DataFrame que contiene los datos.
    - dependent_var: str, el nombre de la variable dependiente (Y).
    - independent_vars: lista de str, nombres de las variables independientes (X).
    - adjust_for_age: bool, si es True, incluye 'Edad' como variable de control.
    
    Retorna:
    - model: modelo OLS ajustado.
    """
    # Selecciona las variables independientes
    X = data[independent_vars].copy()

    # Convertir a dummy solo si la variable es categórica
    for col in X.columns:
        if pd.api.types.is_categorical_dtype(X[col]) or X[col].dtype == 'object':
            X = pd.get_dummies(X, columns=[col], drop_first=True)


    # Agrega una constante para la intersección
    X = sm.add_constant(X)

    # Define la variable dependiente
    y = data[dependent_var]

    # Cuenta y muestra los valores NaN en cada columna
    nan_counts = pd.concat([X, y], axis=1).isna().sum()

    # Elimina filas con NaN o valores infinitos en X o y
    valid_data = pd.concat([X, y], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    X_clean = valid_data[X.columns]
    y_clean = valid_data[dependent_var]

    # Ajusta el modelo
    model = sm.OLS(y_clean, X_clean).fit()

    return model




def select_variables_by_pearson(data, target_var, threshold=0.3):
    """
    Calcula el coeficiente de Pearson y el p-value entre la variable dependiente
    y cada variable independiente. Retorna las variables cuyo p-value es > 0.3
    o el coeficiente de Pearson es < -0.3.

    Parameters:
    - data: DataFrame con los datos.
    - target_var: str, nombre de la variable dependiente.
    - threshold: float, valor umbral para p-value o coeficiente de Pearson.

    Returns:
    - selected_vars: lista de variables que cumplen la condición.
    """
    selected_vars = []
    pearson_coef = []


    # Filtrar solo las columnas numéricas
    numerical_data = data.select_dtypes(include=['float64', 'int64'])

    for col in numerical_data.columns:
        if col != target_var:  # Evitar correlacionar la variable consigo misma
            clean_data = numerical_data[[col, target_var]].dropna()

            if len(clean_data) >= 2:
                coef, p_value = pearsonr(clean_data[col], clean_data[target_var])

                # Verificar la condición de p-value o coeficiente de correlación
                if coef > threshold or coef < -threshold and p_value < 0.05:
                    selected_vars.append(col)
                    pearson_coef.append(coef)

    return pearson_coef, selected_vars


filenames = ['median_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_PVC0.csv',
                   'median_qCBF_StandardSpace_TotalWM_PVC0.csv',
                   'median_qCBF_StandardSpace_Hammers_PVC0.csv',
                   'median_qCBF_StandardSpace_WholeBrain_PVC0.csv',
                   'median_qCBF_StandardSpace_HOsub_CONN_PVC0.csv',
                   'median_qCBF_StandardSpace_MNI_Structural_PVC0.csv',
                   'median_qCBF_StandardSpace_DeepWM_2024_PVC0.csv',
                   'median_qCBF_StandardSpace_HOcort_CONN_2024_PVC0.csv',
                   'median_qCBF_StandardSpace_TotalGM_2024_PVC0.csv',
                   'median_qCBF_StandardSpace_Thalamus_2024_PVC0.csv'
                   ]

filenames_2 = ['mean_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_PVC2.csv',
                   'mean_qCBF_StandardSpace_TotalWM_PVC2.csv',
                   'mean_qCBF_StandardSpace_Hammers_PVC2.csv',
                   #'mean_qCBF_StandardSpace_WholeBrain_PVC2.csv',
                   'mean_qCBF_StandardSpace_HOsub_CONN_PVC2.csv',
                   'mean_qCBF_StandardSpace_MNI_Structural_PVC2.csv',
                   'mean_qCBF_StandardSpace_DeepWM_2024_PVC2.csv',
                   'mean_qCBF_StandardSpace_HOcort_CONN_2024_PVC2.csv',
                   'mean_qCBF_StandardSpace_TotalGM_2024_PVC2.csv',
                   'mean_qCBF_StandardSpace_Thalamus_2024_PVC2.csv'
                   ]
templates = ['mean_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_PVC2.csv',
            'mean_qCBF_StandardSpace_MNI_Structural_PVC2.csv',
            'mean_qCBF_StandardSpace_Hammers_PVC2.csv',
             'mean_qCBF_StandardSpace_HOcort_CONN_2024_PVC2.csv',
            'mean_qCBF_StandardSpace_HOsub_CONN_PVC2.csv',
            ]

ids_to_exclude = ['CP0062', 'CP0105', 'CP0216',
                  'CP0075', 'CP0076', 'CP0078','CP0079', 'CP0140','CP0196']

kruskal_wallis = False

#ids_to_exclude = ['CP0011', 'CP0015'] ?????
# IDs que quieres ignorar (artefactos)

significative_rois = []
#ids_to_exclude = ['CP0062', 'CP0105',
#                  'CP0044', 'CP0053',
#                  'CP0054', 'CP0061',
#                  'CP0107', 'CP0108',
#                  'CP0123', 'CP0142',
#                  'CP0147', 'CP0154']
# IDs que quieres ignorar (artefactos + vasculares)

df_variables_participants = pd.read_csv(variables_participants_csv)
df_white_matter_total = pd.read_csv(wmh_total_csv)
df_white_matter_total = df_white_matter_total.rename(columns={'Number of clusters of WMH': 'Number of clusters of WMH total',
                                                              'Total volume of clusters': 'Total volume of clusters total'})


df_white_matter_deep = pd.read_csv(wmh_deep_csv)
df_white_matter_deep = df_white_matter_deep.rename(columns={'Number of clusters of WMH': 'Number of clusters of WMH deep',
                                                            'Total volume of clusters': 'Total volume of clusters deep'})

columns_total= df_white_matter_total.columns[1:]
columns_deep = df_white_matter_deep.columns[1:]

df_merged = df_variables_participants.merge(df_white_matter_total, on='ID', how='left')
df_merged = df_merged.merge(df_white_matter_deep, on='ID', how='left')

# for roi in columns_deep:
#     # Mantener todas las columnas menos las que estamos iterando, pero conservar la roi actual
#     df_merged_filtered = df_merged.loc[:,
#                          ~df_merged.columns.isin(columns_total) | (df_merged.columns == roi)].copy()
#     try:
#         pearson_coef, variables = select_variables_by_pearson(df_merged_filtered, roi, threshold=0.3)
#         if variables:
#             print(roi, variables,'\n', pearson_coef)
#     except KeyError:
#         pass
#
# for roi in columns_total:
#     # Mantener todas las columnas menos las que estamos iterando, pero conservar la roi actual
#     df_merged_filtered = df_merged.loc[:,
#                          ~df_merged.columns.isin(columns_total) | (df_merged.columns == roi)].copy()
#     try:
#         pearson_coef, variables = select_variables_by_pearson(df_merged_filtered, roi, threshold=0.3)
#         if variables:
#             print(roi, variables,'\n', pearson_coef)
#     except KeyError:
#         pass



for i, template in enumerate(templates):
    print(templates[i])
    roi_df = pd.read_csv(os.path.join(output_data, template))
    roi_df_filtered = roi_df[~roi_df['ID'].isin(ids_to_exclude)]
    columns_to_iterate = roi_df_filtered.columns[1:-1]  # Excluye la primera y la última columna

    n_long_covid = len(roi_df_filtered[roi_df_filtered['Grupo'] == 'Covid Prolongado'])
    n_control = len(roi_df_filtered[roi_df_filtered['Grupo'] == 'Control'])
    n_total = len(roi_df_filtered)



    # Limpiar columnas con más de 10 NaN (estas van a ser las columnas en las que no tenemos información del CBF
    columns_with_more_than_10_nan = roi_df_filtered.columns[roi_df_filtered.isna().sum() > 50]
    roi_df_filtered = roi_df_filtered.drop(columns = columns_with_more_than_10_nan)
    columns_to_iterate = list(set(columns_to_iterate)- set(columns_with_more_than_10_nan))

    df_merged = roi_df_filtered.merge(df_merged, on='ID', how='left', suffixes=('', '_duplicado'))
    df_merged = df_merged.loc[:, ~df_merged.columns.str.endswith('_duplicado')]

    all_columns = list(df_merged.columns)

    for roi in columns_to_iterate:
        model = fit_linear_model(df_merged, roi, ['Grupo', 'Edad', 'Género', 'MOCA valor de referencia'])
        grupo_p_values = model.pvalues['Grupo_Covid Prolongado']

        if grupo_p_values < 0.05:
            print(roi)
            for variable, p_value in model.pvalues.items():
                print(f"{variable}: {p_value:.5f} ", end='')  # Imprimir cada variable y su valor p en la misma línea
            print()  # Salto de línea al final

        if kruskal_wallis:
            grouped_data = [group[roi].dropna().values for name, group in roi_df_filtered.groupby('Grupo')]
            stat, p_value = kruskal(*grouped_data)
            if p_value < 0.05:


                for name, group in roi_df_filtered.groupby('Grupo'):
                    ids_with_nan_in_group = group[group[
                        roi].isna()].index.tolist()
                    if ids_with_nan_in_group:
                        print(f"Grupo: {name}, IDs con NaN: {ids_with_nan_in_group}")

                print(f'Kruskal-Wallis para {roi}: estadístico={stat}, p-valor={p_value}')

                characters = ['(', ')']

                name_roi = [part.replace(char, '') for part in roi.split(',') for char in characters]

                # Crear un boxplot usando seaborn
                plt.figure(figsize=(4, 6))
                sns.boxplot(x='Grupo', y=roi,
                            data=roi_df_filtered)
                plt.title(f'Boxplot para {roi}')
                plt.xlabel('Grupo')
                plt.ylabel('qCBF')



            output_plot_template = os.path.join(output_plots, output_plots_atlas[i])

            if not os.path.exists(output_plot_template):
                os.makedirs(output_plot_template)

            output_filename = os.path.join(output_plot_template, f'Boxplot_{name_roi[0]}.png')


            plt.savefig(output_filename, format='png', dpi=300)  # dpi=300 para alta resolución






# # Mantener todas las columnas menos las que estamos iterando, pero conservar la roi actual
        # df_merged_filtered = df_merged.loc[:,
        #                      ~df_merged.columns.isin(columns_to_iterate) | (df_merged.columns == roi)].copy()
        #
        # try:
        #     pearson_coef, variables = select_variables_by_pearson(df_merged_filtered, roi, threshold=0.5)
        #     if variables:
        #         print(roi, variables, '\n', pearson_coef)
        # except KeyError:
        #     pass
        #

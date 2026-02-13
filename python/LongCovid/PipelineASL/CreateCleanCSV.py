import pandas as pd
import numpy as np
import os,re

participants_csv= '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/VoluntariosProyectoCovidProlongado.csv'

folder_path_bids_1 = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/BIDS/derivatives/ExploreASL/Population/Stats'
folder_path_bids_2 = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/BIDS_2/derivatives/ExploreASL/Population/Stats'
output_data = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/Final_data_csv'



regex_filenames = [r'median_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_.*?_PVC0.tsv',
                   r'median_qCBF_StandardSpace_TotalWM_.*?_PVC0.tsv',
                   r'median_qCBF_StandardSpace_Hammers_.*?_PVC0.tsv',
                   r'median_qCBF_StandardSpace_WholeBrain_.*?_PVC0.tsv',
                   r'median_qCBF_StandardSpace_HOsub_CONN_.*?_PVC0.tsv',
                   r'median_qCBF_StandardSpace_MNI_Structural_.*?_PVC0.tsv',
                   r'median_qCBF_StandardSpace_DeepWM_.*?_PVC0.tsv',
                   r'median_qCBF_StandardSpace_HOcort_CONN_.*?_PVC0.tsv',
                   r'median_qCBF_StandardSpace_TotalGM_.*?_PVC0.tsv',
                   r'median_qCBF_StandardSpace_Thalamus_.*?_PVC0.tsv'
                   ]

regex_filenames_2 = [r'mean_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_.*?_PVC2.tsv',
                   r'mean_qCBF_StandardSpace_TotalWM_.*?_PVC2.tsv',
                   r'mean_qCBF_StandardSpace_Hammers_.*?_PVC2.tsv',
                   r'mean_qCBF_StandardSpace_WholeBrain_.*?_PVC2.tsv',
                   r'mean_qCBF_StandardSpace_HOsub_CONN_.*?_PVC2.tsv',
                   r'mean_qCBF_StandardSpace_MNI_Structural_.*?_PVC2.tsv',
                   r'mean_qCBF_StandardSpace_DeepWM_.*?_PVC2.tsv',
                   r'mean_qCBF_StandardSpace_HOcort_CONN_.*?_PVC2.tsv',
                   r'mean_qCBF_StandardSpace_TotalGM_.*?_PVC0.tsv',
                   r'mean_qCBF_StandardSpace_Thalamus_.*?_PVC0.tsv'
                   ]


output_filenames = ['median_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_PVC0.csv',
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

output_filenames_2 = ['mean_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_PVC2.csv',
                   'mean_qCBF_StandardSpace_TotalWM_PVC2.csv',
                   'mean_qCBF_StandardSpace_Hammers_PVC2.csv',
                   'mean_qCBF_StandardSpace_WholeBrain_PVC2.csv',
                   'mean_qCBF_StandardSpace_HOsub_CONN_PVC2.csv',
                   'mean_qCBF_StandardSpace_MNI_Structural_PVC2.csv',
                   'mean_qCBF_StandardSpace_DeepWM_2024_PVC2.csv',
                   'mean_qCBF_StandardSpace_HOcort_CONN_2024_PVC2.csv',
                   'mean_qCBF_StandardSpace_TotalGM_2024_PVC2.csv',
                   'mean_qCBF_StandardSpace_Thalamus_2024_PVC2.csv'
                   ]

df_participants = pd.read_csv(participants_csv)

for i, regex_filename in enumerate(regex_filenames_2):
    path_filename_bids_1 = None
    path_filename_bids_2 = None

    for filename in os.listdir(folder_path_bids_1):
        if re.match(regex_filename, filename):
            path_filename_bids_1 = os.path.join(folder_path_bids_1, filename)

    for filename in os.listdir(folder_path_bids_2):
        if re.match(regex_filename, filename):
            path_filename_bids_2 = os.path.join(folder_path_bids_2, filename)

    if path_filename_bids_1 and path_filename_bids_2:
        header_1 = pd.read_csv(path_filename_bids_1, header=None, nrows=1, sep='\t')
        header_2 = pd.read_csv(path_filename_bids_1, header=None, skiprows=1, nrows=1, sep='\t')  # Second row

        df1 = pd.read_csv(path_filename_bids_1,skiprows=[1], sep='\t')
        df2 = pd.read_csv(path_filename_bids_2, skiprows=[1], sep='\t')

        merged_df = pd.concat([df1, df2], ignore_index=True, axis=0)
        merged_df.columns = ['_'.join(col).strip() for col in zip(header_1.values.flatten(), header_2.values.flatten())]
        multi_header = pd.MultiIndex.from_arrays([header_1.values.flatten(), header_2.values.flatten()])
        merged_df.columns = multi_header

        merged_df.iloc[:, 0] = merged_df.iloc[:, 0].apply(lambda x: x[4:-2])
        merged_df = merged_df.drop(merged_df.columns[1:5], axis=1)

        # Extraer el segundo nivel del MultiIndex como una nueva fila
        second_level_header = merged_df.columns.get_level_values(1).to_numpy()
        new_row = pd.DataFrame([second_level_header], columns=merged_df.columns)

        # Insertar la nueva fila al inicio del DataFrame
        merged_df = pd.concat([new_row, merged_df], ignore_index=True)

        # Cambiar el índice a una sola columna para poder fusionar
        merged_df.columns = merged_df.columns.get_level_values(0)

        # Realizar la fusión
        final_df = pd.merge(merged_df, df_participants, left_on='participant_id', right_on='ID', how='inner')

        # Eliminar la columna 'ID' y renombrar
        final_df = final_df.drop(columns='ID')
        final_df.rename(columns={'participant_id': 'ID'}, inplace=True)

        # Guardar el resultado
        final_df.to_csv(os.path.join(output_data, output_filenames_2[i]), index=False)

        #
        # #final_df = pd.merge(merged_df, df_participants, left_on=[('participant_id', 'StudyID')], right_on='ID')
        # final_df = pd.merge(merged_df, df_participants, left_on=('participant_id', 'StudyID'), right_on='ID', how='inner')
        # final_df = final_df.drop(columns='ID')
        # final_df.rename(columns={('participant_id', 'StudyID'): 'ID'}, inplace=True)
        #
        # final_df.to_csv(os.path.join(output_data, output_filenames_2[i]), index=False)




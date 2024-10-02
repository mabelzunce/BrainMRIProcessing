import os, csv
from unittest.mock import inplace

import pandas as pd
from pycparser.c_ast import Break

nifti_images = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Nifti/NewNifti"
processed_images = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/NewProcessedBIANCA/"
respuestas_cuestionarios_2_csv = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Respuestas_cuestionarios.csv"
respuestas_cuestionarios_csv= '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/VoluntariosProyectoCovidProlongado.csv'

all_bianca_data_csv = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/all_bianca_data.csv"
imagenes_sin_hallazgos_relevantes_csv = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/imagenes_sin_hallazgos_relevantes.csv"

def create_dataframe(number_thr):
    """Create a dataframe file with headers."""
    headers = ['Subject',
            f'Total Volume thr {number_thr}',
            f'Number Lesions thr {number_thr} cluster 3', f'Total Vol Lesions {number_thr} cluster 3',
            f'Number Lesions thr {number_thr} cluster 5', f'Total Vol Lesions {number_thr} cluster 5',
            f'Number Lesions thr {number_thr} cluster 7', f'Total Vol Lesions {number_thr} cluster 7',
            f'Number Lesions thr {number_thr} cluster 10', f'Total Vol Lesions {number_thr} cluster 10',
            ]
    df_lesions_volumes = pd.DataFrame(columns=headers)
    return df_lesions_volumes



if __name__ == '__main__':

    subjects_to_process = [subject for subject in sorted(os.listdir(processed_images)) if subject.startswith("CP")]
    subjects_to_exclude = ['CP0107','CP0106' 'CP0116', 'CP0163','CP0199']
    num_cluster_csv = 5
    num_thr_csv = 0.7
    deep = False

    bianca_data_thr_cluster = (f"/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/bianca_data_thr_0_"
                               f"{num_thr_csv * 10}_cluster_{num_cluster_csv}.csv")

    bianca_data_thr_cluster_deep = (f"/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/bianca_data_thr_0_"
                               f"{num_thr_csv * 10}_cluster_{num_cluster_csv}_deep.csv")

    # Read respuesta cuestionarios CSV
    df_respuestas_cuestionarios_1 = pd.read_csv(respuestas_cuestionarios_csv)
    df_respuestas_cuestionarios_2 = pd.read_csv(respuestas_cuestionarios_2_csv)

    number_thr = [0.7, 0.8, 0.9]
    num_clusters = [3, 5, 7, 10]

    df_lesions = pd.DataFrame()

    number_lesions = []
    volume_lesions = []

    for subject in subjects_to_process:
        if subject in subjects_to_exclude:
            continue

        df_subject_lesions = pd.DataFrame()

        df_subject_lesions['Subject'] = [subject]

        subject_path = os.path.join(processed_images, subject)
        T2_dir = os.path.join(subject_path,"T2", "lesions")


        if os.path.exists(T2_dir):

            for thr in number_thr:
                T2_volume_file = os.path.join(T2_dir, f"total_volume_0_{int(thr * 10)}.txt")

                if os.path.exists(T2_volume_file):
                    with open(T2_volume_file, 'r') as file_1:
                        lines_total_volume = file_1.readlines()
                        total_volume = int(lines_total_volume[0].strip())

                df_subject_lesions[f'Total Volume thr {thr}'] = [total_volume]

                for num_cluster in num_clusters:
                    T2_cluster = os.path.join(T2_dir, f"volume_cluster_0_{int(thr * 10)}_{num_cluster}.txt")

                    if os.path.exists(T2_cluster):
                        with open(T2_cluster, 'r') as file_2:
                            lines_clusters_and_volumes = file_2.readlines()

                        lines_clusters_and_volumes = [line.strip() for line in lines_clusters_and_volumes]
                        lines_clusters_and_volumes = [line.split() for line in lines_clusters_and_volumes]

                        row_lesions = lines_clusters_and_volumes[0]
                        row_volume_lesions = lines_clusters_and_volumes[1]

                        if row_lesions[-1] != "Cluster":
                            number_lesion = int(row_lesions[-1])
                        else:
                            number_lesion = 0

                        volume_lesion = float(row_volume_lesions[-1])

                        df_subject_lesions[f'Number Lesions thr {thr} cluster {num_cluster}'] = [number_lesion]
                        df_subject_lesions[f'Total Vol Lesions thr {thr} cluster {num_cluster}'] = [volume_lesion]
            if deep:
                T2_cluster_deep = os.path.join(T2_dir, f"volume_cluster_0_7_5_vent_mask.txt")

                if os.path.exists(T2_cluster_deep):
                    with open(T2_cluster_deep, 'r') as file_3:
                        lines_clusters_and_volumes = file_3.readlines()

                    lines_clusters_and_volumes = [line.strip() for line in lines_clusters_and_volumes]
                    lines_clusters_and_volumes = [line.split() for line in lines_clusters_and_volumes]

                    row_lesions = lines_clusters_and_volumes[2]
                    row_volume_lesions = lines_clusters_and_volumes[3]

                    if row_lesions[-1] != "Cluster":
                        number_lesion = int(row_lesions[-1])
                    else:
                        number_lesion = 0

                    volume_lesion = float(row_volume_lesions[-1])

                    df_subject_lesions[f'Number Lesions thr {num_thr_csv} cluster 5 deep'] = [number_lesion]
                    df_subject_lesions[f'Total Vol Lesions thr {num_thr_csv} cluster 5 deep'] = [volume_lesion]



            df_lesions = pd.concat([df_lesions, df_subject_lesions], ignore_index=True)


    # Extract two column and create new df
    df_final_lesions = df_lesions[['Subject',
                                   f'Number Lesions thr {num_thr_csv} cluster {num_cluster_csv}',
                                   f'Total Vol Lesions thr {num_thr_csv} cluster {num_cluster_csv}']]

    # Renaming columns
    if deep:
        df_final_lesions.rename(columns={
            f'Number Lesions thr {num_thr_csv} cluster {num_cluster_csv} deep': 'Number of clusters of WMH',
            f'Total Vol Lesions thr {num_thr_csv} cluster {num_cluster_csv} deep': 'Total volume of clusters',
        }, inplace=True)
    else:
        df_final_lesions.rename(columns={
            f'Number Lesions thr {num_thr_csv} cluster {num_cluster_csv}': 'Number of clusters of WMH',
            f'Total Vol Lesions thr {num_thr_csv} cluster {num_cluster_csv}': 'Total volume of clusters',
        }, inplace=True)

    # Merge with respuestas csv
    df_respuestas_cuestionarios_1 = df_respuestas_cuestionarios_1[['ID', 'Grupo']]
    df_respuestas_cuestionarios_2 = df_respuestas_cuestionarios_2[['ID','Impresión diagnóstica: ']]

    df_respuestas_cuestionarios_2.rename(columns={
        'ID': 'Subject'
    }, inplace = True)


    # Final Lesions
    df_final_lesions = df_respuestas_cuestionarios_1.merge(df_final_lesions, left_on='ID', right_on='Subject', how='inner')

    df_final_lesions = df_final_lesions.drop(columns=['Subject'])
    df_final_lesions_2 = df_respuestas_cuestionarios_2.merge(df_final_lesions, left_on='Subject', right_on='ID',
                                                             how='inner')
    df_final_lesions = df_final_lesions.replace('Covid Prolongado', 'COVID')
    df_final_lesions = df_final_lesions.replace('Control', 'CONTROL')

    df_final_lesions_imagenes_sin_lesiones = df_final_lesions_2.loc[
        df_final_lesions_2['Impresión diagnóstica: '] == 'Imágenes anatómicas de encéfalo sin hallazgos relevantes.']


    # Lesions
    df_lesions = df_respuestas_cuestionarios_1.merge(df_lesions, left_on='ID', right_on='Subject', how='inner')
    df_lesions = df_lesions.drop(columns=['Subject'])
    df_lesions = df_lesions.replace('Covid Prolongado', 'COVID')
    df_lesions = df_lesions.replace('Control', 'CONTROL')


    if deep:
        df_final_lesions.to_csv(bianca_data_thr_cluster_deep, index=False)
    else:
        df_final_lesions.to_csv(bianca_data_thr_cluster, index=False)

    df_final_lesions_imagenes_sin_lesiones.to_csv(imagenes_sin_hallazgos_relevantes_csv, index=False)
    df_lesions.to_csv(all_bianca_data_csv, index=False)








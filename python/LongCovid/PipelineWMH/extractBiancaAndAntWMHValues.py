import os
import pandas as pd


nifti_images = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Nifti/NewNifti"
processed_images = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/NewProcessedBIANCA/"
respuestas_cuestionarios_2_csv = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Respuestas_cuestionarios.csv"
respuestas_cuestionarios_csv= '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/VoluntariosProyectoCovidProlongado.csv'

if __name__ == '__main__':

    subjects_to_process = [subject for subject in sorted(os.listdir(processed_images)) if subject.startswith("CP")]
    subjects_to_exclude = ['CP0001','CP0035','CP0106', 'CP0163', 'CP0199']
    num_cluster_csv = 5
    num_thr_csv = 0.7


    respuestas_cuestionarios_2_csv = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Respuestas_cuestionarios.csv"
    respuestas_cuestionarios_csv = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/VoluntariosProyectoCovidProlongado.csv'

    # Read respuesta cuestionarios CSV
    df_respuestas_cuestionarios_1 = pd.read_csv(respuestas_cuestionarios_csv)
    df_respuestas_cuestionarios_2 = pd.read_csv(respuestas_cuestionarios_2_csv)

    # Outputs
    bianca_total_and_pv_values = os.path.join(processed_images, 'TotalPvAndDeepValuesBianca.csv')
    ant_total_and_pv_values = os.path.join(processed_images, 'TotalPvAndDeepValuesAnt.csv')

    number_thr_bianca = [0.7]
    number_thr_ant = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    df_respuestas_cuestionarios_1 = df_respuestas_cuestionarios_1[['ID', 'Grupo']]

    desired_order = ['Subject', 'Grupo', 'Total', 'Periventricular', 'Deep']

    # BIANCA results
    for thr in number_thr_bianca:
        df_subjects_lesions_bianca = pd.DataFrame()

        bianca_total_and_pv_values = os.path.join(processed_images, f'TotalPvAndDeepValuesBiancaThr0_{int(thr*10)}.csv')

        for subject in subjects_to_process:
            if subject in subjects_to_exclude:
                continue

            subject_path = os.path.join(processed_images, subject)
            T2_lessions_Bianca = os.path.join(subject_path,"T2", "LesionsBianca")


            df_subject_lesions_bianca = pd.DataFrame()
            df_subject_lesions_bianca['Subject'] = [str(subject)]

            # Extract Bianca Values
            if os.path.exists(T2_lessions_Bianca):


                results_Bianca = os.path.join(T2_lessions_Bianca,
                                              "results",
                                              f"deep_and_pv_0_{int(thr * 10)}",
                                              "pwmh_dwmh_output",
                                              "WMH_tot_pvent_deep_10mm.txt")


                if os.path.exists(results_Bianca):


                    with open(results_Bianca, 'r') as file_1:
                        lines_total_volume = file_1.readlines()
                        total_volume_pv_deep = [float(num) for num in lines_total_volume[0].split()]


                        total_volume = total_volume_pv_deep[0]
                        total_pv = total_volume_pv_deep[1]
                        total_deep = total_volume_pv_deep[2]

                        df_subject_lesions_bianca['Total'] = [total_volume]
                        df_subject_lesions_bianca['Periventricular'] = [total_pv]
                        df_subject_lesions_bianca['Deep'] = [total_deep]

            df_subjects_lesions_bianca = pd.concat([df_subjects_lesions_bianca, df_subject_lesions_bianca], ignore_index=True)

        # Merge with df cuestionarios
        df_subjects_lesions_bianca = df_respuestas_cuestionarios_1.merge(df_subjects_lesions_bianca, left_on='ID', right_on='Subject', how='inner')
        df_subjects_lesions_bianca = df_subjects_lesions_bianca.drop(columns=['ID'])
        df_subjects_lesions_bianca = df_subjects_lesions_bianca.replace('Covid Prolongado', 'COVID')
        df_subjects_lesions_bianca = df_subjects_lesions_bianca.replace('Control', 'CONTROL')

        df_subjects_lesions_bianca = df_subjects_lesions_bianca[desired_order]

        df_subjects_lesions_bianca.to_csv(bianca_total_and_pv_values, index=False)


    # ANT results
    for thr in number_thr_ant:
        df_subjects_lesions_ant = pd.DataFrame()

        ANT_total_and_pv_values = os.path.join(processed_images,
                                                  f'TotalPvAndDeepValuesAntThr0_{int(thr * 10)}.csv')



        for subject in subjects_to_process:
            if subject in subjects_to_exclude:
                continue
            subject_path = os.path.join(processed_images, subject)
            T2_lessions_ANTS = os.path.join(subject_path, "T2", "LesionsAntSysu")


            df_subject_lesions_ant = pd.DataFrame()
            df_subject_lesions_ant['Subject'] = [str(subject)]

            # Extract ANT Values
            if os.path.exists(T2_lessions_ANTS):

                results_ANT = os.path.join(T2_lessions_ANTS,
                                              "results",
                                              f"deep_and_pv_0_{int(thr * 10)}",
                                              "pwmh_dwmh_output",
                                              "WMH_tot_pvent_deep_10mm.txt")

                if os.path.exists(results_ANT):
                    with open(results_ANT, 'r') as file_1:
                        lines_total_volume = file_1.readlines()
                        total_volume_pv_deep = [float(num) for num in lines_total_volume[0].split()]

                        total_volume = total_volume_pv_deep[0]
                        total_pv = total_volume_pv_deep[1]
                        total_deep = total_volume_pv_deep[2]

                        df_subject_lesions_ant['Total'] = [total_volume]
                        df_subject_lesions_ant['Periventricular'] = [total_pv]
                        df_subject_lesions_ant['Deep'] = [total_deep]

            df_subjects_lesions_ant = pd.concat([df_subjects_lesions_ant, df_subject_lesions_ant],
                                                   ignore_index=True)

        df_subjects_lesions_ant = df_respuestas_cuestionarios_1.merge(df_subjects_lesions_ant, left_on='ID',
                                                                         right_on='Subject', how='inner')
        df_subjects_lesions_ant = df_subjects_lesions_ant.drop(columns=['ID'])
        df_subjects_lesions_ant = df_subjects_lesions_ant.replace('Covid Prolongado', 'COVID')
        df_subjects_lesions_ant = df_subjects_lesions_ant.replace('Control', 'CONTROL')

        df_subjects_lesions_ant = df_subjects_lesions_ant[desired_order]


        df_subjects_lesions_ant.to_csv(ANT_total_and_pv_values, index=False)

    # Compare
    # ANT 0.7 vs Bianca 0.7
    compare_ant_and_bianca_csv =  os.path.join(processed_images, f'CompareAntAndBianca0_7.csv')
    df_bianca = pd.read_csv(os.path.join(processed_images, f'TotalPvAndDeepValuesBiancaThr0_7.csv'))

    df_ant = pd.read_csv(os.path.join(processed_images,
                                      f'TotalPvAndDeepValuesAntThr0_7.csv'))
    df_compare_ant_and_bianca = pd.merge(df_bianca, df_ant,
                                         on='Subject',
                                         how='inner',
                                         suffixes=('_bianca', '_ant'))

    df_compare_ant_and_bianca.to_csv(compare_ant_and_bianca_csv, index=False)

    compare_ant_and_bianca_2_csv = os.path.join(processed_images, f'CompareAnt0_5AndBianca0_7.csv')

    # ANT 0.5 vs Bianca 0.7
    df_bianca_2 = pd.read_csv(os.path.join(processed_images, f'TotalPvAndDeepValuesBiancaThr0_7.csv'))

    df_ant_2 = pd.read_csv(os.path.join(processed_images,
                                      f'TotalPvAndDeepValuesAntThr0_5.csv'))
    df_compare_ant_and_bianca_2 = pd.merge(df_bianca_2, df_ant_2,
                                           on='Subject',
                                           how='inner',
                                           suffixes=('_bianca', '_ant'))

    df_compare_ant_and_bianca_2.to_csv(compare_ant_and_bianca_2_csv, index=False)
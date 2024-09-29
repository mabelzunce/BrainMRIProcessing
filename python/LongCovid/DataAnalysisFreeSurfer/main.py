import stat_to_csv, norm_to_etiv, statistics_covid, boxplots_covid
import csv_covid
import os
import pandas as pd

if __name__ == '__main__':
    path_FS_COVID_images = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/StructuralProcessed2/PreprocessedMRI/Freesurfer" #"/home/sol/COVID/FS_subjects"
    path_CSV_subjects = "/home/sol/COVID/CSV_subjects"
    path_plots = "/home/sol/COVID/Graficos"

    path_CSV_covid_group = "/home/sol/COVID/CSV_subjects/COVID_group.csv" # Subjects groups

    # Dirs vith FS stats in CSV files
    path_COVID_sanos = os.path.join(path_CSV_subjects, 'Control')
    path_COVID_prolongado = os.path.join(path_CSV_subjects, 'Covid Prolongado')

    # CSV files
    # # Volumes
    path_brain_volumes = os.path.join(path_CSV_subjects, 'brain_volumes.csv')
    path_parcellation = os.path.join(path_CSV_subjects, 'parcellation.csv')
    path_parcellation_DKT = os.path.join(path_CSV_subjects, 'parcellation_DKT.csv')
    path_segmentation = os.path.join(path_CSV_subjects, 'segmentation.csv')

    # # Normalized to eTIV
    path_parcellation_to_etiv = os.path.join(path_CSV_subjects, "parcellation_etiv.csv")
    path_parcellation_DKT_to_etiv = os.path.join(path_CSV_subjects, "parcellation_etiv.csv")
    path_segmentation_to_etiv = os.path.join(path_CSV_subjects, "parcellation_etiv.csv")

    # Thickness
    path_parcellation_thick = os.path.join(path_CSV_subjects, 'parcellation_thick.csv')
    path_parcellation_thick_DKT = os.path.join(path_CSV_subjects, 'parcellation_thick_DKT.csv')

    # Statistics
    path_statistics = os.path.join(path_CSV_subjects,'statistics.csv')


    if not os.path.exists(path_CSV_subjects):
        os.mkdir(path_CSV_subjects)

    # Subjects
    subjects_to_exclude = ['CP0011', 'CP0015', 'CP0144', 'CP0035', 'CP0106']
    subjects = [f for f in os.listdir(path_FS_COVID_images) if os.path.isdir(os.path.join(path_FS_COVID_images, f))]  # list of subjects
    subjects = sorted(subjects)


    # Read CSV with COVID subjects group
    CSV_covid_group = pd.read_csv(path_CSV_covid_group)

    # If covid files
    COVID = True

    # Convert .stat FS files to CSV
    for subject in subjects:
        if subject == "fsaverage":
            continue

        try:
            # Extract COVID GROUP
            group = CSV_covid_group.loc[CSV_covid_group["subject"] == subject]["group"].values[0]

            if group == "Prueba Piloto":
                group = "Control"

            output_path_COVID_subject = os.path.join(path_CSV_subjects, group) + "/"

            #  If it doesn't exist, create a directory with the name of the group
            if not os.path.exists(output_path_COVID_subject):
                os.mkdir(output_path_COVID_subject)

            stat_to_csv.create_csvs(subject, path_FS_COVID_images, output_path_COVID_subject)


        except IndexError: #If i dont know the CP group

            stat_to_csv.create_csvs(subject, path_FS_COVID_images, path_CSV_subjects + "/")

    # Create CSV with all COVID volumes and Cortical Thickness

    # Volumes
    df_brainvol = pd.DataFrame() # Brain Volumes General
    df_parcellation = pd.DataFrame() # Parcellation  (cortical)
    df_parcellation_DKT = pd.DataFrame() # Parcellation DKT (cortical)
    df_segmentation = pd.DataFrame() # Segmentation  (subcortical)

    # Cortical Thickness
    df_parcellation_thick = pd.DataFrame()
    df_parcellation_thick_DKT = pd.DataFrame()

    # Append CONTROL Group to dataframes
    df_brainvol, df_parcellation, df_segmentation, df_parcellation_thick, df_parcellation_DKT, df_parcellation_thick_DKT = csv_covid.create_df(
        path_COVID_sanos, df_brainvol, df_parcellation,
        df_segmentation, df_parcellation_thick, df_parcellation_DKT, df_parcellation_thick_DKT, "Grupo Control")

    # Append COVID Prolongado Group to dataframes
    df_brainvol, df_parcellation, df_segmentation, df_parcellation_thick, df_parcellation_DKT, df_parcellation_thick_DKT = csv_covid.create_df(
        path_COVID_prolongado, df_brainvol, df_parcellation,
        df_segmentation, df_parcellation_thick, df_parcellation_DKT, df_parcellation_thick_DKT, "COVID Prolongado")


    # Volumes: Convert mm3 to cm3 and Save CSV file
    df_brainvol = csv_covid.convert_mm3_to_cm3_brainvol(df_brainvol)
    df_brainvol.index.rename('subject', inplace=True)
    df_brainvol.to_csv(path_brain_volumes, decimal='.')
    df_brainvol.reset_index(drop=False, inplace=True)

    df_parcellation = csv_covid.convert_mm3_to_cm3_parcellation(df_parcellation)
    df_parcellation.index.rename('subject', inplace=True)
    df_parcellation.to_csv(path_parcellation, decimal='.')
    df_parcellation.reset_index(drop=False, inplace=True)

    df_parcellation_DKT = csv_covid.convert_mm3_to_cm3_parcellation(df_parcellation_DKT)
    df_parcellation_DKT.index.rename('subject', inplace=True)
    df_parcellation_DKT.to_csv(path_parcellation_DKT, decimal='.')
    df_parcellation_DKT.reset_index(drop=False, inplace=True)

    df_segmentation = csv_covid.convert_mm3_to_cm3_segmentation(df_segmentation)
    df_segmentation.index.rename('subject', inplace=True)
    df_segmentation.to_csv(path_segmentation, decimal='.')
    df_segmentation.reset_index(drop=False, inplace=True)

    # Cortical Thickness
    df_parcellation_thick.index.rename('subject', inplace=True)
    df_parcellation_thick.to_csv(path_parcellation_thick, decimal='.')
    df_parcellation_thick.reset_index(drop=False, inplace=True)

    df_parcellation_thick_DKT.index.rename('subject', inplace=True)
    df_parcellation_thick_DKT.to_csv(path_parcellation_thick_DKT, decimal='.')
    df_parcellation_thick_DKT.reset_index(drop=False, inplace=True)

    subject_brainvol = list(df_brainvol['subject'])
    extra_subject = [subject for subject in subjects if subject not in subject_brainvol]

    # Normalize volumes to eTIV
    df_parcellation_volumes_to_etiv = norm_to_etiv.parcellation_to_etiv(df_brainvol, df_parcellation)
    df_parcellation_DKT_volumes_to_etiv = norm_to_etiv.parcellation_to_etiv(df_brainvol, df_parcellation_DKT)
    df_segmentation_volumes_to_etiv = norm_to_etiv.segmentation_to_etiv(df_brainvol, df_segmentation)

    # Exclude subjects
    df_brainvol = csv_covid.exclude_subject_by_id(df_brainvol,subjects_to_exclude)
    df_parcellation = csv_covid.exclude_subject_by_id(df_parcellation, subjects_to_exclude)
    df_segmentation = csv_covid.exclude_subject_by_id(df_segmentation, subjects_to_exclude)
    df_parcellation_DKT = csv_covid.exclude_subject_by_id(df_parcellation_DKT, subjects_to_exclude)
    df_parcellation_thick = csv_covid.exclude_subject_by_id(df_parcellation_thick, subjects_to_exclude)
    df_parcellation_thick_DKT = csv_covid.exclude_subject_by_id(df_parcellation_thick_DKT, subjects_to_exclude)

    df_parcellation_volumes_to_etiv = csv_covid.exclude_subject_by_id(df_parcellation_volumes_to_etiv, subjects_to_exclude)
    df_parcellation_DKT_volumes_to_etiv = csv_covid.exclude_subject_by_id(df_parcellation_DKT_volumes_to_etiv,
                                                                      subjects_to_exclude)
    df_segmentation_volumes_to_etiv = csv_covid.exclude_subject_by_id(df_segmentation_volumes_to_etiv,
                                                                          subjects_to_exclude)
    # Do Statistics CSV: all mean and std with p value for t student test and mann whitneyu

    # Parcellation
    # Columns with variables
    df_statistics = pd.DataFrame()

    # Segmentation Statistics
    df_statistics_segmentation = statistics_covid.segmentation_statistics(df_segmentation)
    df_statistics = pd.concat([df_statistics, df_statistics_segmentation])

    # Brain Volumes
    df_statistics_brainvol = statistics_covid.segmentation_statistics(df_brainvol)
    df_statistics = pd.concat([df_statistics, df_statistics_brainvol])

    # Parcellation Statistics
    # Volumes
    df_statistics_parcellation_volumes = statistics_covid.parcellation_statistics(df_parcellation, name_atlas="aparc-a2009s-volume")
    df_statistics = pd.concat([df_statistics, df_statistics_parcellation_volumes])

    df_statistics_parcellation_DKT_volumes = statistics_covid.parcellation_statistics(df_parcellation_DKT,
                                                                     name_atlas="aparc-Desikan-volume")
    df_statistics = pd.concat([df_statistics, df_statistics_parcellation_DKT_volumes])

    # Cortical Thickness
    df_statistics_parcellation_thick = statistics_covid.parcellation_statistics(df_parcellation_thick,
                                                               name_atlas="aparc-a2009s-thickness")
    df_statistics = pd.concat([df_statistics, df_statistics_parcellation_thick])

    df_statistics_parcellation_DKT_thick = statistics_covid.parcellation_statistics(df_parcellation_thick_DKT,
                                                                   name_atlas="aparc-Desikan-thickness")
    df_statistics = pd.concat([df_statistics, df_statistics_parcellation_DKT_thick])

    # Volumes Normalized to eTIV
    df_statistics_parcellation_eTIV = statistics_covid.parcellation_statistics(df_parcellation_volumes_to_etiv,
                                                                                name_atlas="aparc-a2009s-normalized")
    df_statistics = pd.concat([df_statistics, df_statistics_parcellation_eTIV])

    df_statistics_parcellation_DKT_eTIV = statistics_covid.parcellation_statistics(df_parcellation_DKT_volumes_to_etiv,
                                                                               name_atlas="aparc-Desikan-normalized")
    df_statistics = pd.concat([df_statistics, df_statistics_parcellation_eTIV])

    df_statistics_segmentation_eTIV = statistics_covid.segmentation_statistics(df_segmentation_volumes_to_etiv, etiv=True)
    df_statistics = pd.concat([df_statistics, df_statistics_segmentation_eTIV])

    df_statistics.to_csv(path_statistics)

    print()

    # # Boxplots
    # # Segmentation Boxplot
    # path_segmentation_plots = os.path.join(path_plots, 'Segmentation')
    # if not os.path.exists(path_segmentation_plots):
    #     os.mkdir(path_segmentation_plots)
    #
    # regions_segmentation = list(df_segmentation.columns)[1:-2]
    # boxplots_covid.do_boxplot_segmentation(regions_segmentation, df_segmentation, path_segmentation_plots, normalized=False)
    #
    # # Parcellation Boxplot
    # path_parcellation_plots = os.path.join(path_plots, 'Parcellation')
    # if not os.path.exists(path_parcellation_plots):
    #     os.mkdir(path_parcellation_plots)
    #
    # regions_parcellation = list(df_parcellation.columns)[1:-3]
    # boxplots_covid.do_boxplot_parcellation(regions_parcellation, df_parcellation,  path_parcellation_plots,
    #                                                               normalized=False, hemisphere='lh')
    # boxplots_covid.do_boxplot_parcellation(regions_parcellation, df_parcellation, path_parcellation_plots,
    #                                                               normalized=False, hemisphere='rh')
    #
    #
    # # Total GM, Brain Volume and WM
    # regions = ["Brain Segmentation Volume Without Ventricles from Surf", "Estimated Total Intracranial Volume",
    #            "Total gray matter volume"]

    # # Normalized Gray Matter
    #
    # df_statistics = statistics_covid.write_a_row_with_statistics(df_statistics, covid_values, control_values, name_variable, name_atlas, hemisphere='',
    #                             verbose=False)
    # # white_matter_brainvol = ["Total cerebral white matter volume"]
    # # white_matter_segmentation = ["Left-Cerebellum-White-Matter", "Right-Cerebellum-White-Matter", "CC_Posterior"]
    #
    # # Search signficative differences
    # path_signficative_statistics = os.path.join(path_CSV_subjects, 'signficative_statistics.csv')
    # df_signficative_statistics = df_statistics.loc[(df_statistics['p value(test T)'] < 0.05) | (
    #             df_statistics['p value (test Mann Whitneyu)'] < 0.05), :]
    #
    #
    # df_signficative_statistics.to_csv(path_signficative_statistics)


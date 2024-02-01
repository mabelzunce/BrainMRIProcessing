# en la terminal de Ubuntu python3 stat_to_csv.py nombre_de_los_sujetos_que_quiero_parsear

import csv
import os
import sys

import pandas as pd


def aseg_to_csv(file, name_new_file, name_brainvol):
    '''Convert aseg.stats to csv'''

    with open(file, 'r') as f:
        lines = f.readlines()
        etiv = lines[33].strip().split(", ")[2::]
        headers = lines[78].split()[2::]
        data = lines[79::]

    with open(name_new_file, 'w') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(headers)

        for row in data:
            row = row.strip().split()
            csv_writer.writerow(row)

    # Write Estimated Total IntraCraneal Value in brainvol.csv
    with open(name_brainvol, 'a') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(etiv)




    return


def aparc_to_csv(file, name_new_file):
    '''Convert aparc.stats to csv'''

    with open(file, 'r') as f:
        lines = f.readlines()
        headers = lines[60].split()[2::]
        mean_cortical_thick = lines[20].strip().split(", ")[2::]
        data = lines[61::]


    with open(name_new_file, 'w') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(headers)

        for row in data:
            row = row.strip().split()
            csv_writer.writerow(row)


    return

def brainvol_to_csv(file, name_new_file):
    '''Convert brainvol.stats to csv'''

    with open(file, 'r') as f:
        data = f.readlines()

    headers = ['Region', 'Volume', 'Units']

    with open(name_new_file, 'w') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(headers)

        for row in data:
            row = row.strip().split(', ')
            if row[1] == "SupraTentorialVolNotVent":
                row[2] = "Supratentorial volume NotVent"

            if row[2] == "Brain Segmentation Volume":
                total_brain_vm = float(row[3])
            if row[2] == "Total cortical gray matter volume":
                total_gm = float(row[3])


            csv_writer.writerow(row[2::])


        # # normalizacion
        #
        # normalization = total_gm / total_brain_vm * 100
        # row = ["Grey Matter Normalization (notFS)", normalization, "s/u"]
        # csv_writer.writerow(row)

    return

def extract_cortical_thickness(file_lh, file_rh, name_brainvol):
    '''File lh: aparc left hemisphere
        File rh: aparc left hemisphere
        name_brainvol: brainvol csv
        '''

    with open(file_lh, 'r') as f:
        lines = f.readlines()
        mean_cortical_thick_lh = float(lines[20].strip().split(", ")[2::][1])

    with open(file_rh, 'r') as f:
        lines = f.readlines()
        mean_cortical_thick_rh = float(lines[20].strip().split(", ")[2::][1])

    mean_cortical_thick = (mean_cortical_thick_rh + mean_cortical_thick_lh) / 2

    # Write Mean Cortical Thickness in brainvol.csv
    with open(name_brainvol, 'a') as file:
        csv_writer = csv.writer(file)
        mean_cortical_thick_row = ["Mean Thickness", mean_cortical_thick, 'mm']
        csv_writer.writerow(mean_cortical_thick_row)

    return

def create_csvs(subject, path_subject, output_path
                ):
    subjects_dir_stats = path_subject + "/" + subject + "/stats"

    path_aseg = subjects_dir_stats + "/aseg.stats"

    # aparc
    path_lh_aparc = subjects_dir_stats + "/lh.aparc.stats"
    path_rh_aparc = subjects_dir_stats + "/rh.aparc.stats"

    # aparca2009s
    path_lh_aparca2009s = subjects_dir_stats + "/lh.aparc.a2009s.stats"
    path_rh_aparca2009s = subjects_dir_stats + "/rh.aparc.a2009s.stats"

    # aparcaDKT
    path_lh_aparcDKT = subjects_dir_stats + "/lh.aparc.DKTatlas.stats"
    path_rh_aparcDKT = subjects_dir_stats + "/rh.aparc.DKTatlas.stats"


    path_brainvol = subjects_dir_stats + "/brainvol.stats"



    subject_path = os.path.join(output_path, subject)

    if not os.path.exists(subject_path):  # create path of subject
        os.makedirs(subject_path)

    common_path = output_path + subject + '/' + subject

    # brainvol
    brainvol_to_csv(path_brainvol, common_path + '_brainvol.csv')

    aseg_to_csv(path_aseg, common_path + '_aseg.csv', common_path + '_brainvol.csv')
    # aparc
    aparc_to_csv(path_lh_aparc, common_path + '_lh_aparc.csv')
    aparc_to_csv(path_rh_aparc, common_path + '_rh_aparc.csv')

    # aparca2009s
    aparc_to_csv(path_lh_aparca2009s, common_path + '_lh_aparca2009s.csv')
    aparc_to_csv(path_rh_aparca2009s, common_path + '_rh_aparca2009s.csv')

    # aparcDKT
    aparc_to_csv(path_lh_aparcDKT, common_path + '_lh_aparcDKT.csv')
    aparc_to_csv(path_rh_aparcDKT, common_path + '_rh_aparcDKT.csv')

    # Mean cortical thickness
    extract_cortical_thickness(path_lh_aparc, path_rh_aparc, common_path + '_brainvol.csv')

    return


if __name__ == '__main__':
    path_COVID = "/home/sol/COVID/FS_subjects"
    path_ADNI_CN = "/home/sol/ADNI/ADNI_COVID/CN"
    path_ADNI_AD = "/home/sol/ADNI/ADNI_COVID/AD"

    output_path_COVID = "/home/sol/COVID/CSV_subjects"
    output_path_AD = "/home/sol/COVID/CSV_subjects/AD/"
    output_path_CN = "/home/sol/COVID/CSV_subjects/CN/"

    if not os.path.exists(output_path_COVID):
        os.mkdir(output_path_COVID)

    # Subjects
    subjects = [f for f in os.listdir(path_COVID) if os.path.isdir(os.path.join(path_COVID, f))]  # list of subjects
    subjects = sorted(subjects)


    # Read CSV with COVID subjects group
    path_CSV_covid_group = "/home/sol/COVID/CSV_subjects/COVID_group.csv"
    CSV_covid_group = pd.read_csv(path_CSV_covid_group)

    # If covid files
    COVID = True

    for subject in subjects:
        if subject == "fsaverage":
            continue

        try:    #COVID
            if COVID:

                group = CSV_covid_group.loc[CSV_covid_group["subject"] == subject]["group"].values[0]

                if group == "Prueba Piloto":
                    group = "Control"

                output_path_COVID_subject = os.path.join(output_path_COVID, group) + "/"

                if not os.path.exists(output_path_COVID_subject):
                    os.mkdir(output_path_COVID_subject)

                create_csvs(subject, path_COVID, output_path_COVID_subject)

            else:
                create_csvs(subject, path_ADNI_AD, output_path_AD)
        except IndexError:
            print(subject)

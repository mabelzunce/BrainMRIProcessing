# en la terminal de Ubuntu python3 stat_to_csv.py nombre_de_los_sujetos_que_quiero_parsear

import csv
import os
import sys


def aseg_to_csv(file, name_new_file):
    '''Convert aseg.stats to csv'''

    with open(file, 'r') as f:
        lines = f.readlines()
        headers = lines[78].split()[2::]
        data = lines[79::]

    with open(name_new_file, 'w') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(headers)

        for row in data:
            row = row.strip().split()
            csv_writer.writerow(row)

    return


def aparc_to_csv(file, name_new_file):
    '''Convert aparc.stats to csv'''

    with open(file, 'r') as f:
        lines = f.readlines()
        headers = lines[60].split()[2::]
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
            csv_writer.writerow(row[2::])

    return


def create_csvs(subject):
    subjects_dir_stats = "/home/sol/COVID/subjects/" + subject + "/stats"

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

    subject_path = os.path.join("/home/sol/COVID/CSV_subjects/", subject)

    if not os.path.exists(subject_path):  # create path of subject
        os.makedirs(subject_path)

    common_path = "/home/sol/COVID/CSV_subjects/" + subject + '/' + subject
    aseg_to_csv(path_aseg, common_path + '_aseg.csv')
    # aparc
    aparc_to_csv(path_lh_aparc, common_path + '_lh_aparc.csv')
    aparc_to_csv(path_rh_aparc, common_path + '_rh_aparc.csv')

    # aparca2009s
    aparc_to_csv(path_lh_aparca2009s, common_path + '_lh_aparca2009s.csv')
    aparc_to_csv(path_rh_aparca2009s, common_path + '_rh_aparca2009s.csv')

    # aparcDKT
    aparc_to_csv(path_lh_aparcDKT, common_path + '_lh_aparcDKT.csv')
    aparc_to_csv(path_rh_aparcDKT, common_path + '_rh_aparcDKT.csv')

    brainvol_to_csv(path_brainvol, common_path + '_brainvol.csv')

    return


if __name__ == '__main__':
    path = "/home/sol/COVID/subjects/"

    if len(sys.argv) > 1:
        subjects = sys.argv[1::]
    else:
        print(os.listdir(path))
        subjects = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]  # list of subjects

    for subject in subjects:
        create_csvs(subject)

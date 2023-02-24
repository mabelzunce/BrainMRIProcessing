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
    subjects_dir_stats = "/usr/local/freesurfer/7.3.2/subjects/" + subject + "/stats"

    path_aseg = subjects_dir_stats + "/aseg.stats"
    path_lh_aparc = subjects_dir_stats + "/lh.aparc.stats"
    path_rh_aparc = subjects_dir_stats + "/lh.aparc.stats"
    path_brainvol = subjects_dir_stats + "/brainvol.stats"

    if not os.path.exists(subject):  # create path of subject
        os.makedirs(subject)

    common_path = subject + '/' + subject

    aseg_to_csv(path_aseg, common_path + '_aseg.csv')
    aparc_to_csv(path_lh_aparc, common_path + '_lh_aparc.csv')
    aparc_to_csv(path_rh_aparc, common_path + '_rh_aparc.csv')
    brainvol_to_csv(path_brainvol, common_path + '_brainvol.csv')

    return


if __name__ == '__main__':

    if len(sys.argv) > 0:
        subjects = sys.argv[1::]

        for subject in subjects:
            create_csvs(subject)

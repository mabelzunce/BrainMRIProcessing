import SimpleITK as sitk
import numpy as np
import os
import pandas as pd
from collections import Counter
import shutil

subjects_asl_csv =  '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/AllImages4D/subjects_asl_total_covid_4d_2.csv'
variables_participants_csv ='/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/RespuestasCuestionariosCovid_31_10.csv'
output_csv ='/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/subject_covariate.csv'
output_csv_control ='/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/subject_covariate_control.csv'
output_csv_covid ='/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/subject_covariate_covid.csv'


output_images ='/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/CBFimages'

output_control = os.path.join(output_images, 'CONTROL')
output_covid = os.path.join(output_images, 'COVID')

if not os.path.exists(output_control):
    os.mkdir(output_control)

if not os.path.exists(output_covid):
    os.mkdir(output_covid)


# Load data
df_variables_participants = pd.read_csv(variables_participants_csv, skiprows=2)
df_variables_participants.rename(columns={df_variables_participants.columns[0]: "ID"}, inplace=True)
df_subjects_asl = pd.read_csv(subjects_asl_csv)

files_images = os.listdir(output_images)

data = []

data_control = []
data_covid = []

# ID Grupo Edad Género
# Loop through subjects
for index, row in df_subjects_asl.iterrows():
    subject_id = row['Subject']

    # Find the corresponding group
    group_row = df_variables_participants[df_variables_participants['ID'] == subject_id]
    group = group_row['Grupo'].values[0]
    age = group_row['Edad'].values[0]
    genre = group_row['Género'].values[0]

    data.append([subject_id, group, age, genre])
    if genre == 'M':
        genre = 1
    else:
        genre = 0
    if group == "CONTROL":
        data_control.append([subject_id, group, age, genre])
    else:
        data_covid.append([subject_id, group, age, genre])

    # for file in files_images:
    #     if subject_id in file:
    #         CBF_image = os.path.join(output_images, file)
    #
    #         if os.path.exists(CBF_image) and group == 'CONTROL':
    #             shutil.copy(CBF_image, output_control)
    #         elif os.path.exists(CBF_image) and group == 'COVID':
    #             shutil.copy(CBF_image, output_covid)
    #         else:
    #             print(subject_id)

df_output = pd.DataFrame(data, columns=['Subject', 'Group', 'Age', 'Genre'])
df_output_control = pd.DataFrame(data_control, columns=['Subject', 'Group', 'Age', 'Genre'])
df_output_covid = pd.DataFrame(data_covid, columns=['Subject', 'Group', 'Age', 'Genre'])

df_output.to_csv(output_csv, index=False)
df_output_control.to_csv(output_csv_control, index=False)
df_output_covid.to_csv(output_csv_covid, index=False)



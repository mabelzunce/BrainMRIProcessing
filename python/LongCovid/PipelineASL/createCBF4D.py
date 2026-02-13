import SimpleITK as sitk
import numpy as np
import os
import pandas as pd
from collections import Counter
import shutil

images_bids1 = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/BIDS/derivatives/ExploreASL/Population'
images_bids2 = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/BIDS_2/derivatives/ExploreASL/Population'

output_path_image =  '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/asl_total_covid_4d_2.nii.gz'
output_path_csv =  '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/subjects_asl_total_covid_4d_2.csv'
output_path_cbf_images =  '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/CBFimages'

variables_participants_csv ='/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/VoluntariosProyectoCovidProlongado.csv'

pd_variables_participants = pd.read_csv(variables_participants_csv)

ids_to_exclude = ['CP0062', 'CP0105', 'CP0216']

if not os.path.exists(output_path_cbf_images):
    os.mkdir(output_path_cbf_images)

files_1 = os.listdir(images_bids1)
files_2 = os.listdir(images_bids2)

all_files = files_1 + files_2

images = []
subjects = []
groups = []

for file in all_files:
    if file.startswith('qCBF_sub'):
        images.append(file)

images.sort()


# First image
first_image = sitk.ReadImage(os.path.join(images_bids1, images[0]))
first_image_size = first_image.GetSize()

# Crear un array 4D vacío
img_array_4d = np.zeros((*first_image_size, len(images)), dtype=np.float32)  # Ajusta el tipo de datos según sea necesario

for i,image in enumerate(images):
    subject = image[9:15]

    if subject in ids_to_exclude:
        continue

    group = (pd_variables_participants.loc[pd_variables_participants['ID'] == subject, 'Grupo']).squeeze()
    groups.append(group)

    path_image_bids1 = os.path.join(images_bids1, image)
    path_image_bids2 = os.path.join(images_bids2, image)

    if os.path.exists(path_image_bids1):
        path_image = path_image_bids1
    elif os.path.exists(path_image_bids2):
        path_image = path_image_bids2
    else:
        print(subject)
        continue

    img = sitk.ReadImage(path_image)
    img.CopyInformation(first_image)


    if img.GetSize() != first_image.GetSize():
        raise ValueError(f"Las dimensiones de {subject} no coinciden con la primera imagen.")

    img_array_4d[:, :, :, i] = sitk.GetArrayFromImage(img)
    subjects.append(subject)

    shutil.copy(path_image, os.path.join(output_path_cbf_images, image))



df = pd.DataFrame(subjects, columns=["Subject"])
df.to_csv(output_path_csv, index=False)

subject_counts = Counter(groups)

image_4d = sitk.GetImageFromArray(img_array_4d)
image_4d.CopyInformation(first_image)

sitk.WriteImage(image_4d, output_path_image)


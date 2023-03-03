import SimpleITK as sitk
import numpy as np
import pandas as pd

if __name__ == '__main__':

    # A path to a parcellated brain .nii image:
    aparc = '/home/sol/COVID/image_processing/aparc+aseg_CP0001.nii'
    aparc_DKT = '/home/sol/COVID/image_processing/aparc.DKTatlas+aseg_CP0001.nii'
    aparc_2009 = '/home/sol/COVID/image_processing/aparc.a2009s+aseg_CP0001.nii'

    #parcellation = [aparc, aparc_DKT, aparc_2009]

    #for image in parcellation:

    # Read the .nii image containing the volume with SimpleITK:
    sitk_aparc= sitk.ReadImage(aparc_2009) # SimpleITK object

    # and access the numpy array:
    np_aparc = sitk.GetArrayFromImage(sitk_aparc)

    # dimensions
    shape_aparc = np_aparc.shape

    # labels
    labels = np.unique(np_aparc)
    n_labels = len(labels)

    # read csv labels
    labels_path = "/home/sol/COVID/atlas/labels.csv"
    df_labels = pd.read_csv(labels_path)

    columns = ['n_label', 'structure', 'total_vol']
    new_df = pd.DataFrame(columns=columns)

    for n_label in labels:
        maskLabel = np_aparc == n_label
        volume_voxels = np.sum(maskLabel)

        label_row = df_labels.loc[df_labels['n_label'] == n_label] # find row with label in csv
        structure_name = str(label_row['structure'].values)[1:-1]# match with the structure name

        new_row = pd.Series({'n_label': n_label, 'structure': structure_name, 'total_vol': volume_voxels})
        new_df = pd.concat([new_df, new_row.to_frame().T], ignore_index=True)

    new_df.to_csv('/home/sol/COVID/image_processing/compare.csv')


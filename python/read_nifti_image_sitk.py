import SimpleITK as sitk
import numpy as np
import pandas as pd

if __name__ == '__main__':

    # A path to a parcellated brain .nii image:
    aparc = '/home/sol/COVID/image_processing/aparc+aseg_CP0001.nii'
    aparc_DKT = '/home/sol/COVID/image_processing/aparc.DKTatlas+aseg_CP0001.nii'
    aparc_2009 = '/home/sol/COVID/image_processing/aparc.a2009s+aseg_CP0001.nii'

    # Read the .nii image containing the volume with SimpleITK:
    sitk_aparc = sitk.ReadImage(aparc_2009)  # SimpleITK object

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

    # read aseg csv
    aseg_path = "/home/sol/COVID/CSV_subjects/CP0001/CP0001_aseg.csv"
    df_aseg = pd.read_csv(aseg_path)

    # read lh aparc csv
    aseg_lh_aparc = "/home/sol/COVID/CSV_subjects/CP0001/CP0001_lh_aparca2009s.csv"
    df_lh_aparc = pd.read_csv(aseg_lh_aparc)

    # read rh aparc csv
    aseg_rh_aparc = "/home/sol/COVID/CSV_subjects/CP0001/CP0001_rh_aparca2009s.csv"
    df_rh_aparc = pd.read_csv(aseg_rh_aparc)

    # create new df
    columns = ['n_label', 'structure', 'total_vol']
    new_df = pd.DataFrame(columns=columns)

    for n_label in labels:
        maskLabel = np_aparc == n_label
        volume_voxels = np.sum(maskLabel)

        label_row = df_labels.loc[df_labels['n_label'] == n_label]  # find row with label in csv
        structure_name = str(label_row['structure'].values)[2:-2]  # match with the structure name

        new_row = pd.Series({'n_label': n_label, 'structure': structure_name, 'total_vol': volume_voxels})

        if n_label <= 255:  # segmentation
            aseg_structure = df_aseg.loc[df_aseg['SegId'] == n_label]  # find row of the structure
            try:
                vol_stat = (aseg_structure['Volume_mm3'].values)[0]
                segmentation_measure = pd.Series({'vol_stat': vol_stat})
                new_row = pd.concat([new_row, segmentation_measure])
            except IndexError:
                pass

        elif n_label < 12100 : # parcellation lh
            aparc_lh_structure = df_lh_aparc.loc[df_lh_aparc['StructName'] == structure_name[7:]]  # find row of the structure
            vol_stat = aparc_lh_structure['GrayVol'].values[0]
            parcellation_measure = pd.Series({'vol_stat': vol_stat})
            new_row = pd.concat([new_row, parcellation_measure])

        else : # parcellation rh
            aparc_rh_structure = df_rh_aparc.loc[df_lh_aparc['StructName'] == structure_name[7:]]  # find row of the structure
            vol_stat = aparc_rh_structure['GrayVol'].values[0]
            parcellation_measure = pd.Series({'vol_stat': vol_stat})
            new_row = pd.concat([new_row, parcellation_measure])

        new_df = pd.concat([new_df, new_row.to_frame().T], ignore_index=True)

    new_df["Error"] = ((new_df["total_vol"] - new_df["vol_stat"]) / new_df['vol_stat']) * 100
    new_df.to_csv('/home/sol/COVID/image_processing/compare.csv')

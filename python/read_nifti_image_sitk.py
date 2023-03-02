import SimpleITK as sitk
import numpy as np

if __name__ == '__main__':

    # A path to a parcellated brain .nii image:
    aparc = '/home/sol/COVID/image_processing/aparc+aseg_CP0001.nii'

    # Read the .nii image containing the volume with SimpleITK:
    sitk_aparc= sitk.ReadImage(aparc)

    # and access the numpy array:
    np_aparc = sitk.GetArrayFromImage(sitk_aparc)

    # dimensions
    shape_aparc = np_aparc.shape

    # labels
    labels = np.unique(np_aparc)
    n_labels = len(labels)

    for label in labels:
        maskLabel = np_aparc == label
        volume_voxels = np.sum(maskLabel)
        print(label, volume_voxels)
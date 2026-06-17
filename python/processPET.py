import os
import SimpleITK as sitk
import ImageRegistration as reg
import matplotlib.pyplot as plt


if __name__ == '__main__':
    path_subject_MRI = "/home/sol/Subject/Procesado/14_t1_mprage_sag_p3_iso.nii.gz"

    path_subject_PET = "/home/sol/Subject/PET/"
    output_subject = "/home/sol/Subject/Procesado"

    # create output dir
    output_subject_PET = "/home/sol/Subject/Procesado/PET"

    sliceToShow = 160

    if not os.path.exists(output_subject_PET):
        os.mkdir(output_subject_PET)

    dir_PET_image = os.listdir(path_subject_PET) # every dir has a pet image
    paths_PET_image = [] # list of path of all pet image
    names_PET_image = [] # list of names of all pet image

    for dir in dir_PET_image:
        path_dir_PET_images = path_subject_PET + "/" + dir
        PET_image = os.listdir(path_dir_PET_images)[0]

        # add PET image to the list
        path_PET_image = path_dir_PET_images + "/" + PET_image
        paths_PET_image.append(path_PET_image)
        names_PET_image.append(PET_image[0:-7]) # delete .nii.gz extension

        print("Image in dir ..: ", PET_image)

    # print(paths_PET_image)
    # print(names_PET_image)

    # read PET images with SimpleITK
    PET_images_SITK = []
    for PET_image in paths_PET_image:
        new_image = sitk.ReadImage(PET_image)
        new_image = sitk.Cast(new_image, sitk.sitkFloat32) # cast
        PET_images_SITK.append(new_image)

    # read MRI image with SimpleITK
    MRI_image = sitk.ReadImage(path_subject_MRI)
    MRI_image = sitk.Cast(MRI_image, sitk.sitkFloat32) # cast


    # Registration
    ref_PET_Image = PET_images_SITK[0] #reference
    for i in range(1, len(PET_images_SITK)):
        resultReg = reg.RigidImageRegistration(PET_images_SITK[i], ref_PET_Image, printLog=True)
        PET_images_SITK[i] = resultReg['image']

    # Compute the sum
    sumImage = PET_images_SITK[0]
    for image in PET_images_SITK[1:]:
        sumImage = sitk.Add(sumImage, image)

    # Register to t1 using the sum
    resultReg = reg.RigidImageRegistration(sumImage, MRI_image, printLog=True)

    sumImage = resultReg['image']
    txPet2Mri = resultReg['tx']

    # Write sum image:
    sitk.WriteImage(sumImage, output_subject_PET + "Subject_sum_reg_t1.nii.gz")

    imageArraySum = sitk.GetArrayFromImage(sumImage)
    imageArrayT1 = sitk.GetArrayFromImage(MRI_image)
    sliceSum = imageArraySum[sliceToShow, :, :]
    sliceT1 = imageArrayT1[sliceToShow, :, :]

    plt.imshow(sliceT1, cmap='gray', alpha=1)
    plt.imshow(sliceSum, cmap='hot', alpha=0.4)
    plt.savefig(output_subject_PET + "pet_t1")

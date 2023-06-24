import os
import dicom2nifti
import dicom_parser


def write_line_freesurfer(path_file, name_subject, path_nifti):
    with open(path_file, 'a') as f:
        f.write(
            f'recon-all -s {name_subject} -i {path_nifti} -all -parallel '
            f'-openmp 16\n')
    pass


def write_line_fastsurfer(path_file, name_subject, path_nifti):
    with open(path_file, 'a') as f:
        f.write(
            f'sudo docker run -v /home/sol/COVID:/data \n'
            f'-v /home/sol/COVID:/output \n'
            f'-v /home/sol:/fs_license \n'
            f'--rm --user $(id -u):$(id -g) deepmi/fastsurfer:latest \n'
            f'  --fs_license /fs_license/license.txt \n'
            f' --t1 /data/t1_mprage_1x1x1.nii.gz \n'
            f'--device cpu \n'
            f'--sid prueba --sd /output\n'
            f'--parallel\n'
            f'--parallel\n'
        )
    pass


if __name__ == '__main__':
    pathCOVIDNifti = "/home/sol/COVID/MRIs/Nifti"
    #pathCOVIDProcessed = "/home/sol/COVID/MRIs/Processed"

    outputCOVID = "/home/sol/COVID/subjects"

    mriDataToProcess = ['t1_mprage_1x1x1', 't1_mprage_1x1x1_estric_']

    images_COVID_processed = sorted(os.listdir(outputCOVID))
    images_COVID_Nifti = sorted(os.listdir(pathCOVIDNifti))
    images_COVID_to_process = []

    # cuales son las imagenes a procesar
    for subject in images_COVID_Nifti:
        if not subject in images_COVID_processed:
            images_COVID_to_process.append(subject)

    # create a script to run FastSurfer
    FreeSurferPath = "/home/sol/COVID/script/freesurferScript2.sh"
    #FastSurferPath = "/home/sol/ADNI/ADNIproc2/ADNIFastSurfer/fastsurferADNI.sh"

    # clear the data in the FreeSurfer file and write first line
    with open(FreeSurferPath, 'w') as file:
        file.write('#!/bin/bash\n')

    # add a line for every subject
    for MRI_subject in images_COVID_to_process:
        outputPath = "/home/sol/COVID/script/freesurferScript2.sh"
        path_nifti = pathCOVIDNifti + "/" + MRI_subject
        path_image_mri = None
        # archivos en la carpeta nifti
        files = os.listdir(path_nifti)

        # busco la imagen NIFTI
        for possible_mri in mriDataToProcess:
            mri = possible_mri + ".nii.gz"

            if mri in files:
                path_image_mri = path_nifti + "/" + mri

        if path_image_mri == None:
            print("Error!! ", MRI_subject)
            continue

        write_line_freesurfer(outputPath, MRI_subject, path_image_mri)


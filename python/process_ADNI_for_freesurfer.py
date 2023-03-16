import os
import dicom2nifti

def subjects_in_dir(path):
    '''Return a list of subjects into the ADNI directory'''

    subjects = []
    files = sorted(os.listdir(pathADNIdicom))  # list of files in path in alphabetical order

    for file in files:
        path_subject = path + '/' + file

        if os.path.isdir(path_subject):
            subjects.append(file)

    return subjects


def MRI_in_subject(path_subject, names_MRI):
    """
        Check if structural MRI is available

        path_subject: path of the subject
        names_MRI: list of possible names for structural MRI images. 
        return: name of MRI or False
        """
    files = os.listdir(path_subject)
    for file in files:
        if file in names_MRI:
            return file
    return False


if __name__ == '__main__':
    pathADNIproject = "/home/sol/ADNI"
    pathADNIdicom = "/home/sol/ADNI/ADNIdicom"

    mriDataToProcess = ['Accelerated_Sagittal_MPRAGE']
    ADNI_subjects = subjects_in_dir(pathADNIdicom)

    # Create output dir
    outputPath = pathADNIproject + "/ADNIproc/"
    if not os.path.exists(outputPath):
        os.mkdir(outputPath)
    else:
        pass

    # Search for every subject with structural MRI availabLe
    # and append  in a dict (subject: subject's name, date: MRI's date, path_subject: MRI's path)

    MRI_subjects = []

    for subject in ADNI_subjects:
        path_subject = pathADNIdicom + "/" + subject
        name_MRI = MRI_in_subject(path_subject, mriDataToProcess)

        if name_MRI:
            path_subject += "/" + name_MRI

            dates_MRI = os.listdir(path_subject)
            for date in dates_MRI:
                path_date_subject = path_subject + "/" + date
                MRI_subjects.append({"subject": subject, "date": date, "path_subject": path_date_subject})

    # Create a dir for every structural MRI
    for MRI_subject in MRI_subjects:

        # create a dir for the subject (if it doesn't exist)
        outputSubject = outputPath + MRI_subject["subject"]
        if not os.path.exists(outputSubject):
            os.mkdir(outputSubject)
        else:
            pass

        # create a dir for every date
        outputSubject += "/" + MRI_subject["date"]
        if not os.path.exists(outputSubject):
            os.mkdir(outputSubject)
        else:
            pass

        DICOM_file = MRI_subject["path_subject"] + "/" + os.listdir(MRI_subject["path_subject"])[0]
        print(f'Convert {DICOM_file} to Nifti in the path {outputSubject}')
        dicom2nifti.convert_directory(DICOM_file, outputSubject, compression=True, reorient=True)



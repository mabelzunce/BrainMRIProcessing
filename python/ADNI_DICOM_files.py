import os
def files_in_path(path_ADNI, path_txt):
    ''' Create a txt file (ADNIfiles) with the name of the DICOM files in the specified path '''

    # clear the data in the txt file
    with open(path_txt, 'w') as file:
        pass

    files = sorted(os.listdir(path_ADNI))  # list of files in dir in alphabetical order

    for file in files:
        path_subject = path_ADNI + '/' + file

        if os.path.isdir(path_subject):
            image = os.listdir(path_subject)  # image

            with open(path_txt, 'a') as f:
                f.write(file + ": " + ','.join(image) + "\n")

    return


if __name__ == "__main__":
    pathADNIdicom = "/home/sol/ADNI/ADNIdicom"
    pathTxt = "/home/sol/ADNI/ADNIfiles.txt"

    files_in_path(pathADNIdicom, pathTxt)

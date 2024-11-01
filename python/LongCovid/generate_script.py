import os

if __name__ == '__main__':
    nifti_images = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Nifti"
    processed_files = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ProcessedBIANCA/"

    subjects_to_process = sorted(os.listdir(nifti_images))

    # clear the data in the script and write first line
    with open("Pipeline/pipeline_script.sh", 'w') as file:
        file.write('#!/bin/bash\n')

    for subject in subjects_to_process:
        subject_path = os.path.join(nifti_images, subject)

        if os.path.isdir(subject_path):

            # add a line in script
            new_line = f'python3 main.py {subject_path} -o {processed_files}\n'

            with open("Pipeline/pipeline_script.sh", 'a') as file:
                file.write(new_line)

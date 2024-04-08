import argparse
from tools import loggingTool as LT
import os

processed_files = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/StructuralProcessed/CP0001"
subject_dir = "/home/sol/COVID/Seguro"

def main():
    parser = argparse.ArgumentParser(description='Structural Pipeline')
    parser.add_argument("subjectDir", help='Subject Directory')
    parser.add_argument('--outputDir', '-o', help='Subject Directory')

    argsa = parser.parse_args()

    subjectDir = argsa.subjectDir
    outputDir = argsa.outputDir

    subject = subjectDir.split("/")[-1]

    outputSubjectDir = os.path.join(argsa.outputDir, subject)
    if not os.path.exis ts(outputSubjectDir):
        os.mkdir(outputSubjectDir)

    # Start Logging
    logger = LT.initLogging(__file__, subject, outputDir)

    # Find T1 image
    t1_suffixes = ["t1_mprage_1x1x1.nii.gz", "t1_mprage_1x1x1_estric_.nii.gz"]  # Possible T1 Suffix
    t1 = None

    for t1_suffix in t1_suffixes:
        possible_image = os.path.join(subjectDir, t1_suffix)

        if os.path.exists(possible_image):
            t1 = possible_image
            break

    if t1 == None:
        logger.error('There is no T1. Subject ' + subject + ' cannot be processed.')
    else:
        LT.runCommand(logger, f'chmod +x ./structural_pipeline/struct_init') # Give permission to Bash Script
        LT.runCommand(logger, f'./structural_pipeline/struct_init {t1} {outputSubjectDir}')

    LT.finishLogging(logger)

    # Find T2 image
    t2_suffixes = ["t2_space_dark_fluid_sag_p3_iso.nii.gz", "t2_space_dark_fluid_sag_p3_iso_estric_.nii.gz"]  # Possible T2 Suffix
    t2 = None

    for t2_suffix in t2_suffixes:
        possible_image = os.path.join(subjectDir, t2_suffix)

        if os.path.exists(possible_image):
            t2 = possible_image
            break

    if t2 == None:
        logger.error('There is no T2. Subject ' + subject + ' cannot be processed.')
    else:
        LT.runCommand(logger, f'chmod +x ./white_matter_pipeline/white_matter')  # Give permission to Bash Script
        LT.runCommand(logger, f'./white_matter_pipeline/white_matter {t2} {outputSubjectDir}')

    LT.finishLogging(logger)


if __name__ == "__main__":
    main()

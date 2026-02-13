import os
import subprocess
import shutil
import tensorflow as tf
import ants
import antspynet


def create_dir(path_dir):
    if not os.path.exists(path_dir):
        os.mkdir(path_dir)
    pass

subjects_nifti_dir = "/home/sol/COVID/Nifti"
subjects_processed_dir = "/home/sol/COVID/image_processing/StructuralProcessed"

t1_suffix0 = "t1_mprage_1x1x1.nii.gz"
t1_suffix1 = "t1_mprage_1x1x1_estric_.nii.gz"
t2_suffix0 = "t2_space_dark_fluid_sag_p3_iso.nii.gz"
t2_suffix1 = "t2_space_dark_fluid_sag_p3_iso_estric_.nii.gz"


# ANTs template
template_ants_0 = "/home/sol/BrainTools/Templates/Oasis/T_template0.nii.gz"
template_ants_1 = "/home/sol/BrainTools/Templates/Oasis/T_template0_BrainCerebellumProbabilityMask.nii.gz"


subjects = sorted(os.listdir(subjects_nifti_dir))
subjects = [subject for subject in subjects if os.path.isdir(os.path.join(subjects_nifti_dir, subject))]

for subject in subjects:
    dir_subject = os.path.join(subjects_nifti_dir, subject)
    mri_files_subject = os.listdir(dir_subject)

    t1, t2 = None, None

    # T1 MRI image
    if t1_suffix0 in mri_files_subject:
        t1 = os.path.join(dir_subject, t1_suffix0)

    elif t1_suffix1 in mri_files_subject:
        t1 = os.path.join(dir_subject, t1_suffix1)

    elif subject == "CP0015":
        t1 = os.path.join(dir_subject, "t1_mprage_1x1x1_s0005.nii.gz")

    elif subject == "CP0074":
        t1 = os.path.join(dir_subject, "t1_mprage_1x1x1_s0007.nii.gz")

    # T2 MRI image
    if t2_suffix0 in mri_files_subject:
        t2 = os.path.join(dir_subject, t2_suffix0)
    elif t2_suffix1 in mri_files_subject:
        t2 = os.path.join(dir_subject, t2_suffix1)


    # create subject dir in image processing dir
    image_processing_dir = os.path.join(subjects_processed_dir, subject)
    if not os.path.exists(image_processing_dir):
        os.mkdir(image_processing_dir)

    # Freesurfer SAMSEG
    samsegdir = os.path.join(image_processing_dir, "Freesurfer")
    if not os.path.exists(samsegdir):
        os.mkdir(samsegdir)
        
    # MRI_coreg command (register t2 to t1)
    reg_mri_coreg = os.path.join(samsegdir, 'flairToT1.lta')
    vol2vol = os.path.join(samsegdir, 'flair_reg.nii')

    mri_coreg_command = f'mri_coreg --mov {t2} --ref {t1} --reg {reg_mri_coreg} --threads 16'
    mri_vol2vol_command = f'mri_vol2vol --mov {t2} --reg {reg_mri_coreg} --o {vol2vol} --targ {t1} '

    if not os.path.exists(reg_mri_coreg):
       subprocess.run([mri_coreg_command], shell=True)

    if not os.path.exists(vol2vol):
       subprocess.run([mri_vol2vol_command], shell=True)

    # run samseg
    samseg_output = os.path.join(samsegdir, 'SAMSEG')
    samseg_command = (f'run_samseg --input {t1} {vol2vol} --pallidum-separate --lesion --lesion-mask-pattern 0 1 '
                      f'--output {samseg_output} --threads 16')

    if not os.path.exists(samseg_output):
        os.mkdir(samseg_output)

    if not os.listdir(samseg_output):
        subprocess.run([samseg_command], shell=True)

    # ANTs Brain extraction
    ants_dir = os.path.join(image_processing_dir, 'ANT')
    if not os.path.exists(ants_dir):
        os.mkdir(ants_dir)

    # T1
    ants_t1_dir = os.path.join(ants_dir, f"T1_{subject}_")
    ants_t1_brain_extracted = f'{ants_t1_dir}BrainExtractionBrain.nii.gz'
    command_ants_brain_extraction_t1 = (f'antsBrainExtraction.sh -d 3 -a {t1} -e {template_ants_0} -m {template_ants_1} '
                                        f'-o {ants_t1_dir} -c 3x1x2x3')

    if not os.path.exists(ants_t1_brain_extracted):
        subprocess.run([command_ants_brain_extraction_t1], shell=True)

    # T2
    ants_t2_dir = os.path.join(ants_dir, f"T2_{subject}_")

    ants_n4_bias_field_correction_image = f'{ants_t2_dir}N4BiasFieldCorrection.nii.gz'
    ants_t2_brain_extracted_mask = f'{ants_t2_dir}BrainExtractionMask.nii.gz'
    ants_t2_brain_extracted_brain = f'{ants_t2_dir}BrainExtractionBrain.nii.gz'

    fsl_maths_command = f'fslmaths {ants_n4_bias_field_correction_image} -mul {ants_t2_brain_extracted_mask} {ants_t2_brain_extracted_brain}'

    # Segmentation
    ants_segmentation_dir = os.path.join(ants_dir, 'Segmentation')
    ants_kmeans_dir = os.path.join(ants_segmentation_dir, 'KMeans')
    ants_wmh = os.path.join(ants_segmentation_dir, 'WMH')

    # KMeans
    ants_t1_kmeans_segmentation =os.path.join(ants_kmeans_dir, f'KMeans_{subject}_T1.nii.gz')
    ants_t2_kmeans_segmentation =os.path.join(ants_kmeans_dir, f'KMeans_{subject}_T2.nii.gz')

    ants_wmh_segmentations_sysu = os.path.join(ants_wmh, f'WMH_segmentation_sysu_{subject}.nii.gz')
    ants_wmh_segmentations_hypermapper = os.path.join(ants_wmh, f'WMH_segmentation_hypermapper_{subject}.nii.gz')
    ants_wmh_segmentations_antsxnet = os.path.join(ants_wmh, f'WMH_segmentation_antsxnet_{subject}.nii.gz')

    # Process T2 image
    if t2 != None and t1 != None:

        ants_t2_image = ants.image_read(vol2vol, dimension=3)

        # N4 bias field correction
        if not os.path.exists(ants_n4_bias_field_correction_image):
            ants_t2_bias_field_correction = ants.n4_bias_field_correction(ants_t2_image)
            ants.image_write(ants_t2_bias_field_correction, ants_n4_bias_field_correction_image)
        else:
            ants_t2_bias_field_correction = ants.image_read(ants_n4_bias_field_correction_image)


        # Brain Mask
        if not os.path.exists(ants_t2_brain_extracted_mask):
            ants_t2_brain_extracted_image = antspynet.utilities.brain_extraction(ants_t2_bias_field_correction, modality="t2", verbose=True)
            ants.image_write(ants_t2_brain_extracted_image, ants_t2_brain_extracted_mask)
        else:
            ants_t2_brain_extracted_image = ants.image_read(ants_t2_brain_extracted_mask)

        # Apply Brain Mask
        if not os.path.exists(ants_t2_brain_extracted_brain):
            subprocess.run([fsl_maths_command], shell=True)

        # Segmentation
        create_dir(ants_segmentation_dir)
        ants_t1 = ants.image_read(t1)
        ants_t1_image_processed = ants.image_read(ants_t1_brain_extracted)
        ants_t2_image_processed = ants_t2_brain_extracted_image

        # TO DO: K Means Segmentation
        create_dir(ants_kmeans_dir)

        # WMH
        create_dir(ants_wmh)

        # Sysu
        ants_wmh_segmentation_image_sysu = antspynet.sysu_media_wmh_segmentation(ants_t2_image, t1=ants_t1,verbose=True)
        ants.image_write(ants_wmh_segmentation_image_sysu, ants_wmh_segmentations_sysu)

        # Hypermapp3r
        ants_wmh_segmentation_image_hypermapper =  antspynet.utilities.white_matter_hyperintensity_segmentation.hypermapp3r_segmentation(ants_t1_image_processed, ants_t2_image_processed, do_preprocessing=False, verbose=True)
        ants.image_write(ants_wmh_segmentation_image_hypermapper, ants_wmh_segmentations_hypermapper)

        # ANTsXNet
        # ants_t1_image_resample = ants.resample_image(ants_t1_image_processed, (64, 64, 64), use_voxels=True)
        # ants_t2_image_resample  = ants.resample_image(ants_t2_image_processed, (64, 64, 64), use_voxels=True)
        #
        # ants_wmh_segmentation_image_antsxnet = antspynet.wmh_segmentation(ants_t2_image_resample, ants_t1_image_resample,do_preprocessing=False,
        #                                                                   use_combined_model=True, verbose=True)
        #
        # ants.image_write(ants_wmh_segmentation_image_antsxnet, ants_wmh_segmentations_antsxnet)



    # FSL FAST
    fsl_dir = os.path.join(image_processing_dir, 'FSL')

    create_dir(fsl_dir)

    # T1
    fsl_t1_dir = os.path.join(fsl_dir, 'T1')
    fsl_t1_brain_extracted = os.path.join(fsl_t1_dir, f'T1_{subject}_BrainExtractionBrain.nii.gz')
    fsl_fast_command_t1 = f'fast -t 1 -n 3 - o fast_T1_{subject} -g --verbose {fsl_t1_brain_extracted}'

    create_dir(fsl_t1_dir)

    # Copy Brain extracted image and run fsl fast
    try:
        if not os.path.exists(fsl_t1_brain_extracted):
            shutil.copyfile(ants_t1_brain_extracted, fsl_t1_brain_extracted)

        if len(os.listdir(fsl_t1_dir)) < 2:
            subprocess.run([fsl_fast_command_t1], shell=True)


    except FileNotFoundError:
        print(f'T1: Error in {subject}')

    # T2
    fsl_t2_dir = os.path.join(fsl_dir, 'T2')
    fsl_t2_brain_extracted = os.path.join(fsl_t2_dir, f'T2_{subject}_BrainExtractionBrain.nii.gz')
    fsl_fast_command_t2 = f'fast -t 2 -n 3 - o fast_T2_{subject} -g --verbose {fsl_t2_brain_extracted}'

    if not os.path.exists(fsl_t2_dir):
        os.mkdir(fsl_t2_dir)

    # Copy Brain extracted image and run fsl fast
    try:
        if not os.path.exists(fsl_t2_brain_extracted):
            shutil.copyfile(ants_t2_brain_extracted_brain, fsl_t2_brain_extracted)


        if len(os.listdir(fsl_t2_dir)) < 2:
            subprocess.run([fsl_fast_command_t2], shell=True)
            print('Fast:ok')

    except FileNotFoundError:
        print(f'T2: Error in {subject}')
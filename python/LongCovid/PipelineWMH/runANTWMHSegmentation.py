import os, time, subprocess
import antspynet, ants
import numpy as np
import SimpleITK as sitk
import tensorflow as tf
import matplotlib.pyplot as plt
from fsl.utils.fslsub import output

# Set environment variables for TensorFlow and ITK to use all available threads
os.environ["TF_NUM_INTRAOP_THREADS"] = str(os.cpu_count())
os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(os.cpu_count())
tf.keras.backend.clear_session()

processed_images = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/NewProcessedBIANCA"

subjects_to_process = [subject for subject in sorted(os.listdir(processed_images)) if subject.startswith("CP")]

for subject in subjects_to_process:
    subject_dir = os.path.join(processed_images, subject)

    t1_image = os.path.join(subject_dir, 'T1', 'T1_brain.nii.gz')
    t2_image = os.path.join(subject_dir, 'T2', 'T2_FLAIR_brain.nii.gz')

    old_output_dir_path = os.path.join(subject_dir, 'T2', 'ANT SYSU')
    output_dir_path = os.path.join(subject_dir, 'T2', 'LesionsAntSysu')

    mask_sysu_path = os.path.join(output_dir_path, 'mask')
    results_sysu_path = os.path.join(output_dir_path, 'results')



    output_probabilistic_map = os.path.join(mask_sysu_path, "output_sysu_wmh_map.nii.gz")

    if os.path.exists(t1_image) and os.path.exists(t2_image):

        if os.path.exists(old_output_dir_path):
            os.rename(old_output_dir_path, output_dir_path)

        if not os.path.exists(output_dir_path):
            os.mkdir(output_dir_path)

        if not os.path.exists(mask_sysu_path):
            os.mkdir(mask_sysu_path)

        if not os.path.exists(results_sysu_path):
            os.mkdir(results_sysu_path)

        t1 = ants.image_read(t1_image)
        t2 = ants.image_read(t2_image)

    #    gm_mask = ants.image_read(gm_mask_image)  # Read the GM mask

        if not os.path.exists(output_probabilistic_map):
            print(f"Subject {subject}: start processing")
            start = time.time()
            probability_mask_wmh= antspynet.sysu_media_wmh_segmentation(t2, t1=t1, use_ensemble=True, verbose=False)
            end = time.time()

            print("Elapsed time: ", end - start)

            ants.image_write(probability_mask_wmh, output_probabilistic_map)

        if os.path.exists(
                os.path.join(subject_dir, 'T2', 'LesionsBianca' 'T1_unbiased_ventmask.nii.gz')):
            vent_mask = os.path.join(subject_dir, 'T2', 'LesionsBianca', 'T1_unbiased_ventmask.nii.gz')
        elif os.path.join(subject_dir, 'T2', 'LesionsBianca' 'mask', 'T1_unbiased_ventmask.nii.gz'):
            vent_mask = os.path.join(subject_dir, 'T2', 'LesionsBianca', 'mask', 'T1_unbiased_ventmask.nii.gz')
        else:
            vent_mask = None


        for i in range(5,10):

            probability_mask_wmh = ants.image_read(output_probabilistic_map)

            thr = i / 10
            output_wmh_mask_image = os.path.join(mask_sysu_path, f"output_sysu_wmh_0_{i}.nii.gz")
            results_deep_and_pv = os.path.join(results_sysu_path, f"deep_and_pv_0_{i}")

            if not os.path.exists(results_deep_and_pv):
                os.mkdir(results_deep_and_pv)

            wmh_mask = ants.threshold_image(probability_mask_wmh, thr, 1.1, 1, 0)
            ants.image_write(wmh_mask, output_wmh_mask_image)


            if vent_mask:
                # TO DO: Use that: https://open.win.ox.ac.uk/pages/fsl/fslpy/
                bianca_command = [
                    "bianca_perivent_deep",
                    output_wmh_mask_image,
                    vent_mask,
                    str(5),
                    str(2),
                    results_deep_and_pv
                ]

                bianca_command = f"bianca_perivent_deep {output_wmh_mask_image} {vent_mask} 5 2 {results_deep_and_pv}"
                subprocess.run(bianca_command, shell=True)



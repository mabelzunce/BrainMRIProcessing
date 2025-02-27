import os, time, subprocess
import antspynet, ants
import numpy as np
import SimpleITK as sitk
import tensorflow as tf
from fsl.wrappers.bianca import bianca_perivent_deep
import matplotlib.pyplot as plt
from fsl.utils.fslsub import output

def create_dir_ant(subject_dir, name_dir):
    output_dir = str(os.path.join(subject_dir, 'T2', name_dir))
    mask_dir = os.path.join(output_dir, 'mask')
    results_dir = os.path.join(output_dir, 'results')

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if not os.path.exists(mask_dir):
        os.mkdir(mask_dir)

    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    return mask_dir, results_dir

def mask_and_results(output_probabilistic_map, mask_dir, results_dir, name_cnn, vent_mask_path):

    for i in range(1, 10):
        thr = i / 10

        if os.path.exists(output_probabilistic_map):
            probability_mask_wmh = ants.image_read(output_probabilistic_map)
            output_wmh_mask_image = os.path.join(mask_dir, f"output_{name_cnn}`_wmh_0_{i}.nii.gz")

            if not os.path.exists(output_wmh_mask_image):
                wmh_mask = ants.threshold_image(probability_mask_wmh, thr, 1.1, 1, 0)
                ants.image_write(wmh_mask, output_wmh_mask_image)

            results_deep_and_pv = os.path.join(results_dir, f"deep_and_pv_0_{i}")

            if os.path.exists(vent_mask_path):
                if not os.path.exists(os.path.join(output_wmh_mask_image, 'pwmh_dwmh_output')):
                    bianca_perivent_deep(output_wmh_mask_image, vent_mask, 5, results_deep_and_pv, do_stats=2)

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

    if os.path.exists(t1_image) and os.path.exists(t2_image):

        mask_dir_sysu, results_dir_sysu = create_dir_ant(subject_dir, 'LesionsAntSysu')
        mask_dir_hippermapper, results_dir_hippermapper = create_dir_ant(subject_dir, 'LesionsAntHippermaper')
        mask_dir_shiva, results_dir_shiva = create_dir_ant(subject_dir, 'LesionsAntShiva')
        mask_dir_antsxnet, results_dir_antsxnet = create_dir_ant(subject_dir, 'LesionsAntsXnet')


        output_probabilistic_map_sysu = os.path.join(mask_dir_sysu, "output_sysu_wmh_map.nii.gz")
        output_probabilistic_map_hipermapper = os.path.join(mask_dir_hippermapper, "output_sysu_hypermapper_map.nii.gz")
        output_probabilistic_map_shiva = os.path.join(mask_dir_shiva, "output_sysu_shiva_map.nii.gz")
        output_probabilistic_map_antsxnet = os.path.join(mask_dir_antsxnet, "output_sysu_antsxnet_map.nii.gz")

        output_probabilistic_map = output_probabilistic_map_sysu


        t1 = ants.image_read(t1_image)
        t2 = ants.image_read(t2_image)


        # Sysu
        if not os.path.exists(output_probabilistic_map_sysu):
            print(f"Subject {subject} Sysu: start processing")
            start = time.time()
            probability_mask_wmh= antspynet.sysu_media_wmh_segmentation(t2, t1=t1, use_ensemble=True, verbose=False)
            end = time.time()

            print("Elapsed time: ", end - start)

            ants.image_write(probability_mask_wmh, output_probabilistic_map_sysu)

        if not os.path.exists(output_probabilistic_map_hipermapper):
            print(f"Subject {subject}  HyperMapper : start processing")
            start = time.time()

            probability_mask_wmh= antspynet.hypermapp3r_segmentation(t1, t2, verbose=False)
            end = time.time()

            print("Elapsed time: ", end - start)

            ants.image_write(probability_mask_wmh, output_probabilistic_map_hipermapper)


        if not os.path.exists(output_probabilistic_map_antsxnet):
            t1_resample = ants.resample_image(t1, (240, 240, 64), use_voxels=True)
            t2_resample = ants.resample_image(t2, (240, 240, 64), use_voxels=True)

            print(f"Subject {subject}  AntsXNet : start processing")
            start = time.time()

            probability_mask_wmh = antspynet.wmh_segmentation(t2_resample, t1_resample, verbose=False)
            end = time.time()

            print("Elapsed time: ", end - start)

            ants.image_write(probability_mask_wmh, output_probabilistic_map_antsxnet)

        vent_mask = os.path.join(subject_dir, 'T2', 'LesionsBianca', 'mask', 'T1_unbiased_ventmask.nii.gz')


        for i in range(1,10):
            thr = i /10

            output_probabilistic_map_sysu = os.path.join(mask_dir_sysu, "output_sysu_wmh_map.nii.gz")
            output_probabilistic_map_hipermapper = os.path.join(mask_dir_hippermapper,
                                                                "output_hypermapper_map.nii.gz")
            output_probabilistic_map_shiva = os.path.join(mask_dir_shiva,
                                                          "output_shiva_map.nii.gz")
            output_probabilistic_map_antsxnet = os.path.join(mask_dir_antsxnet,
                                                          "output_antsxnet_map.nii.gz")

            mask_and_results(output_probabilistic_map_sysu, mask_dir_sysu, results_dir_sysu, 'sysu', vent_mask)
            mask_and_results(output_probabilistic_map_hipermapper, mask_dir_hippermapper, 'hippermapper', vent_mask)
            mask_and_results(output_probabilistic_map_antsxnet, mask_dir_hippermapper, results_dir_antsxnet, 'antsxnet', vent_mask)

            if os.path.exists(output_probabilistic_map_sysu):
                probability_mask_wmh_sysu = ants.image_read(output_probabilistic_map_sysu)
                output_sysu_wmh_mask_image = os.path.join(mask_dir_sysu, f"output_sysu_wmh_0_{i}.nii.gz")

                if not os.path.exists(output_sysu_wmh_mask_image):
                    wmh_mask_sysu = ants.threshold_image(probability_mask_wmh_sysu, thr, 1.1, 1, 0)
                    ants.image_write(wmh_mask_sysu, output_sysu_wmh_mask_image)

                results_deep_and_pv_sysu = os.path.join(results_dir_sysu, f"deep_and_pv_0_{i}")

                if os.path.exists(vent_mask):
                    if not os.path.exists(os.path.join(results_deep_and_pv_sysu, 'pwmh_dwmh_output')):
                        bianca_perivent_deep(output_sysu_wmh_mask_image, vent_mask, 5, results_deep_and_pv_sysu, do_stats=2)

            if os.path.exists(output_probabilistic_map_hipermapper):
                probability_mask_wmh_hipermapper = ants.image_read(output_probabilistic_map_hipermapper)
                output_hipermapper_wmh_mask_image = os.path.join(mask_dir_hippermapper,f"output_hippermapper_wmh_0_{i}.nii.gz")

                if not os.path.exists(output_hipermapper_wmh_mask_image):
                    wmh_mask_hippermapper = ants.threshold_image(probability_mask_wmh_hipermapper, thr, 1.1, 1, 0)
                    ants.image_write(wmh_mask_hippermapper, output_hipermapper_wmh_mask_image)

                results_deep_and_pv_hippermapper = os.path.join(results_dir_hippermapper, f"deep_and_pv_0_{i}")

                if os.path.exists(vent_mask):
                    if not os.path.exists(os.path.join(results_deep_and_pv_sysu, 'pwmh_dwmh_output')):
                        bianca_perivent_deep(output_hipermapper_wmh_mask_image, vent_mask, 5, results_deep_and_pv_hippermapper,
                                             do_stats=2)

    #
        #     thr = i / 10
        #     output_wmh_mask_image = os.path.join(mask_sysu_path, f"output_sysu_wmh_0_{i}.nii.gz")
        #     results_deep_and_pv = os.path.join(results_sysu_path, f"deep_and_pv_0_{i}")
        #
        #     if not os.path.exists(os.path.join(results_deep_and_pv, 'pwmh_dwmh_output')):
        #         if not os.path.exists(results_deep_and_pv):
        #             os.mkdir(results_deep_and_pv)
        #
        #
        #
        #
        #         if vent_mask:
        #             bianca_perivent_deep(output_wmh_mask_image, vent_mask, 5, results_deep_and_pv, do_stats=2)

        # if not os.path.exists(output_probabilistic_map_shiva):
        #     print(f"Subject {subject}  Shiva : start processing")
        #     start = time.time()
        #
        #     #probability_mask_wmh= antspynet.shiva_wmh_segmentation(t2, t1, which_model = "all", verbose=True)
        #     end = time.time()
        #
        #     print("Elapsed time: ", end - start)
        #
        #     ants.image_write(probability_mask_wmh, output_probabilistic_map_shiva)


        # if not os.path.exists(
        #         vent_mask):
        #     vent_mask = None

        # for i in range(1,10):
        #
        #     probability_mask_wmh = ants.image_read(output_probabilistic_map_sysu)
        #
        #     thr = i / 10
        #     output_wmh_mask_image = os.path.join(mask_sysu_path, f"output_sysu_wmh_0_{i}.nii.gz")
        #     results_deep_and_pv = os.path.join(results_sysu_path, f"deep_and_pv_0_{i}")
        #
        #     if not os.path.exists(os.path.join(results_deep_and_pv, 'pwmh_dwmh_output')):
        #         if not os.path.exists(results_deep_and_pv):
        #             os.mkdir(results_deep_and_pv)
        #
        #         wmh_mask = ants.threshold_image(probability_mask_wmh, thr, 1.1, 1, 0)
        #         ants.image_write(wmh_mask, output_wmh_mask_image)
        #
        #
        #         if vent_mask:
        #             bianca_perivent_deep(output_wmh_mask_image, vent_mask, 5, results_deep_and_pv, do_stats=2)



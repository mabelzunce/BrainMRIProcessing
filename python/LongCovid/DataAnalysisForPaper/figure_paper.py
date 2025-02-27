import SimpleITK as sitk
import numpy as np
import pandas as pd
from nilearn import plotting, datasets, surface
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


path_mni_152 ="/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Figure_Paper/MNI152_T1_1mm_brain.nii.gz"
path_aparc_image = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/UTEC_images/068_S_0127/New/Aseg_Norm_MNI_152_ANT.nii.gz"
path_mask_image_volume = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Figure_Paper/Aparc_volume.nii.gz"
path_mask_image_thickness = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Figure_Paper/Aparc_thickness.nii.gz"
path_mask_image_volume_binary = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Figure_Paper/Aparc_volume_binary.nii.gz"
path_mask_image_thickness_binary = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Figure_Paper/Aparc_thickness_binary.nii.gz"

# Path Plots
path_plot_volume = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Figure_Paper/plot_volume.png"
path_plot_thickness = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Figure_Paper/plot_thickness.png"



# Volumes
lh_lingual_gyrus_label = 1013
lh_inferior_parietal_label = 1008
lh_cerebellum_label = 8
rh_cerebellum_label = 47

lh_lingual_gyrus_p_value_volume = 0.05
lh_inferior_parietal_p_value_volume = 0.03
lh_cerebellum_p_value_volume  = 0.03
rh_cerebellum_p_value_volume  = 0.03

volume_values = [lh_lingual_gyrus_label,
                 lh_inferior_parietal_label,
                 lh_cerebellum_label,
                 rh_cerebellum_label]

volume_p_values =[lh_lingual_gyrus_p_value_volume,
                  lh_inferior_parietal_p_value_volume,
                  lh_cerebellum_p_value_volume,
                  rh_cerebellum_p_value_volume]

# Cortical Thickness
lh_postcentral_gyrus_label = 1022
lh_precuneus_label = 1025
lh_supramarginal_gyrus_label = 1031
rh_postcentral_gyrus_label = 2022
rh_superior_temporal_gyrus_label = 2030
rh_precuneus_label = 2025

lh_lingual_gyrus_p_value_thickness = 0.003
lh_postcentral_gyrus_p_value_thickness = 0.02
lh_precuneus_p_value_thickness = 0.01
lh_supramarginal_gyrus_p_value_thickness = 0.03
rh_postcentral_gyrus_p_value_thickness = 0.03
rh_superior_temporal_gyrus_p_value_thickness = 0.049
rh_precuneus_p_value_thickness = 0.01

cortical_thickness_values = [lh_lingual_gyrus_label,
                             lh_postcentral_gyrus_label,
                             lh_precuneus_label,
                             lh_supramarginal_gyrus_label,
                             rh_postcentral_gyrus_label,
                             rh_superior_temporal_gyrus_label,
                             rh_precuneus_label]

cortical_thickness_p_value = [lh_lingual_gyrus_p_value_thickness,
                              lh_postcentral_gyrus_p_value_thickness,
                              lh_precuneus_p_value_thickness,
                              lh_supramarginal_gyrus_p_value_thickness,
                              rh_postcentral_gyrus_p_value_thickness,
                              rh_superior_temporal_gyrus_p_value_thickness,
                              rh_precuneus_p_value_thickness]

segmentation_image = sitk.ReadImage(path_aparc_image)
array_segmentation_image = sitk.GetArrayFromImage(segmentation_image)

# Masks
array_volume = np.zeros_like(array_segmentation_image)
array_volume_binary = np.zeros_like(array_segmentation_image)
array_thickness = np.zeros_like(array_segmentation_image)
array_thickness_binary = np.zeros_like(array_segmentation_image)

# Volume
for i,volume in enumerate(volume_values):
    mask_label = array_segmentation_image == volume
    array_volume[mask_label] = volume_p_values[i]
    array_volume_binary[mask_label] = 1

# Thickness
for i,thick in enumerate(cortical_thickness_values):
    mask_label = array_segmentation_image == thick
    array_thickness[mask_label] = cortical_thickness_p_value[i]
    array_thickness_binary[mask_label] = 1

volume_image =  sitk.GetImageFromArray(array_volume)
thickness_image =  sitk.GetImageFromArray(array_thickness)
volume_binary_image =  sitk.GetImageFromArray(array_volume_binary)
thickness_binary_image =  sitk.GetImageFromArray(array_thickness_binary)



volume_image.CopyInformation(segmentation_image)
thickness_image.CopyInformation(segmentation_image)
volume_binary_image.CopyInformation(segmentation_image)
thickness_binary_image.CopyInformation(segmentation_image)

sitk.WriteImage(volume_image, path_mask_image_volume)
sitk.WriteImage(thickness_image, path_mask_image_thickness)
sitk.WriteImage(volume_binary_image, path_mask_image_volume_binary)
sitk.WriteImage(thickness_binary_image, path_mask_image_thickness_binary)


# Plots
# Crear un mapa de colores que vaya de azul a transparente
colors = [(0, 0, 1), (0, 0, 1, 0)]  # (R, G, B, Alpha)
cmap = LinearSegmentedColormap.from_list("blue_transparent", colors)

# Definir los colores para el colormap NIH
colors = [
    (0, 0, 0, 1),       # Negro
    (0.1, 0.1, 0.4, 1), # Azul oscuro
    (0.2, 0.2, 0.6, 1), # Azul
    (0.4, 0.4, 0.8, 1), # Azul claro
    (0.6, 0.6, 1, 1),   # Casi blanco
    (1, 1, 1, 1)        # Blanco
]
cmap_nih= LinearSegmentedColormap.from_list("nih_colormap", colors)
cmap_jet = plt.cm.jet_r


# Crear el colormap
# Crear el colormap (yellow)
colors = [(1, 1, 0), (1, 1, 0, 0)]  # (R, G, B, Alpha)
cmap_y = LinearSegmentedColormap.from_list("yellow_transparent", colors)


plotting.plot_stat_map(path_mask_image_volume, bg_img=path_mni_152,
                       annotate=False,
                       cmap=cmap_jet,black_bg=False, cut_coords=(-18,4,19), display_mode="ortho",
                       threshold = 0, vmax=0.05, dim= 0, alpha=1.0)
plt.savefig(path_plot_volume)




plotting.plot_stat_map(path_mask_image_thickness, bg_img=path_mni_152,
                       annotate=False,
                       cmap=cmap_jet,black_bg=False, cut_coords=(-18,4,19), display_mode="ortho",
                       threshold = 0, vmax=0.05, dim= 0, alpha=1.0)
plt.savefig(path_plot_thickness)


path_plot_volume_nih = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Figure_Paper/plot_volume_nih.png"
path_plot_thickness_nih = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/Figure_Paper/plot_thickness_nih.png"

plotting.plot_stat_map(path_mask_image_volume, bg_img=path_mni_152,
                       annotate=False,
                       cmap=cmap_nih,black_bg=False, cut_coords=(-18,4,19), display_mode="ortho",
                       threshold = 0, vmax=0.05, dim= 0, alpha=1.0)
plt.savefig(path_plot_volume_nih)


plotting.plot_stat_map(path_mask_image_thickness, bg_img=path_mni_152,
                       annotate=False,
                       cmap=cmap_nih,black_bg=False, cut_coords=(-18,4,19), display_mode="ortho",
                       threshold = 0, vmax=0.05, dim= 0, alpha=1.0)
plt.savefig(path_plot_thickness_nih)

import time
import numpy as np
import SimpleITK as sitk
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np


def view_slice_all(image_original=None, mask_ant_original=None, mask_bianca_original=None,
               slice_idx=107, alpha = 0.7):
    fig, axs = plt.subplots(3, 3, figsize=(12, 12))

    for i in range(3):
        print(i)
        if i == 0:
            image = np.flipud(image_original[slice_idx, :, :])
            mask_bianca = np.flipud(mask_bianca_original[slice_idx, :, :])
            mask_ant = np.flipud(mask_ant_original[slice_idx, :, :])
        elif i == 1:
            image = np.flipud(image_original[:, slice_idx, :])
            mask_bianca = np.flipud(mask_bianca_original[:, slice_idx, :])
            mask_ant = np.flipud(mask_ant_original[:, slice_idx, :])
        elif i == 2:
            image = np.flipud(image_original[:,:,slice_idx])
            mask_bianca = np.flipud(mask_bianca_original[:, :, slice_idx])
            mask_ant = np.flipud(mask_ant_original[:, :, slice_idx,])

        # Mostrar la imagen original con la máscara 1
        axs[i][0].imshow(image, cmap='gray', alpha=1.0)  # Mostrar imagen original
        color_mask_bianca = np.zeros((*image.shape, 4))  # Crear imagen RGBA para la máscara 1
        color_mask_bianca[..., 3] = 0  # Inicialmente, el canal alfa (transparencia) es 0 para todo
        color_mask_bianca[mask_bianca == 1] = [1, 0, 0, alpha]  # Rojo para máscara 1
        color_mask_bianca[mask_bianca == 0] = [1, 1, 1, 0]  # Transparente para el resto
        axs[i][0].imshow(color_mask_bianca)  # Mostrar máscara 1
        axs[i][0].set_title(f'Slice {slice_idx} - Image with BIANCA')
        axs[i][0].axis('off')

        # Mostrar la imagen original con la máscara 2
        axs[i][1].imshow(image, cmap='gray', alpha=1.0)  # Mostrar imagen original
        color_mask_ant = np.zeros((*image.shape, 4))  # Crear imagen RGBA para la máscara 2
        color_mask_ant[..., 3] = 0  # Inicialmente, el canal alfa (transparencia) es 0 para todo
        color_mask_ant[..., 3] = 0
        color_mask_ant[mask_ant == 1] = [0, 0, 1, alpha]  # Azul para máscara 2
        color_mask_ant[mask_ant == 0] = [1, 1, 1, 0]  # Transparente para el resto
        axs[i][1].imshow(color_mask_ant)  # Mostrar máscara 2
        axs[i][1].set_title(f'Slice {slice_idx} - Image with ANTs CNN (Sysu)')
        axs[i][1].axis('off')

        # Corte de la segunda imagen
        axs[i][2].imshow(image, cmap='gray')
        axs[i][2].set_title(f'Image 2 - Slice {slice_idx}')
        axs[i][2].axis('off')


def view_slice_one_slice(image_original=None,
                               slice_idx=107,
                               mask_ant_original=None,
                               mask_bianca_original=None,
                               type="a",
                               alpha=None,):
    fig, axs = plt.subplots(1, 3, figsize=(12, 6))  # Crear dos subgráficos

    # Voltear la imagen verticalmente
    if type == "a":
      image = np.flipud(image_original[slice_idx, :, :])
      mask_ant= np.flipud(mask_ant_original[slice_idx, :, :])
      mask_bianca = np.flipud(mask_bianca_original[slice_idx, :, :])
    elif type == "c":
      image = np.flipud(image_original[:, slice_idx, :])
      mask_ant = np.flipud(mask_ant_original[:, slice_idx, :])
      mask_bianca = np.flipud(mask_bianca_original[:, slice_idx, :])
    elif type == "s":
      image  = np.flipud(image_original[:,:,slice_idx])
      mask_ant = np.flipud(mask_ant_original[:, :, slice_idx])
      mask_bianca = np.flipud(mask_bianca_original[:, :, slice_idx])



    # Mostrar la imagen original con la máscara 2 (invertida)
    axs[0].imshow(image, cmap='gray', alpha=1.0)  # Mostrar imagen original
    color_mask_bianca = np.zeros((*image.shape, 4))  # Crear imagen RGBA para la máscara 2
    color_mask_bianca[..., 3] = 0  # Inicialmente, el canal alfa (transparencia) es 0 para todo
    color_mask_bianca[mask_bianca == 1] = [0, 0, 1, alpha]  # Azul para máscara 2
    color_mask_bianca[mask_bianca == 0] = [1, 1, 1, 0]  # Transparente para el resto
    axs[0].imshow(color_mask_bianca)  # Mostrar máscara 2
    axs[0].set_title(f'Slice {slice_idx} - Image with BIANCA')
    axs[0].axis('off')

    # Mostrar la imagen original con la máscara 1
    axs[1].imshow(image, cmap='gray', alpha=1.0)  # Mostrar imagen original
    color_mask_ant = np.zeros((*image.shape, 4))  # Crear imagen RGBA para la máscara 1
    color_mask_ant[..., 3] = 0  # Inicialmente, el canal alfa (transparencia) es 0 para todo
    color_mask_ant[mask_ant == 1] = [1, 0, 0, alpha]  # Rojo para máscara 1
    color_mask_ant[mask_ant == 0] = [1, 1, 1, 0]  # Transparente para el resto
    axs[1].imshow(color_mask_ant)  # Mostrar máscara 1
    axs[1].set_title(f'Slice {slice_idx} - Image with ANTs CNN (Sysu)')
    axs[1].axis('off')

    # Corte de la segunda imagen
    axs[2].imshow(image, cmap='gray')
    axs[2].set_title(f'Image 2 - Slice {slice_idx}')
    axs[2].axis('off')




if __name__ == "__main__":
    path_image_t2 = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/NewProcessedBIANCA/CP0087/T2/T2_FLAIR_brain.nii.gz'
    path_mask_ant = ('/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/NewProcessedBIANCA/'
                     'CP0087/T2/LesionsAntSysu/mask/output_sysu_wmh_0_5.nii.gz')
    path_mask_bianca =('/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/NewProcessedBIANCA/'
                       'CP0087/T2/LesionsBianca/mask/final_mask_0_7.nii.gz')

    image_t2 = sitk.ReadImage(path_image_t2)
    image_mask_ant = sitk.ReadImage(path_mask_ant)
    image_mask_bianca = sitk.ReadImage(path_mask_bianca)

    image_array_t2 = sitk.GetArrayFromImage(image_t2)
    image_array_mask_ant = sitk.GetArrayFromImage(image_mask_ant)
    image_array_mask_bianca = sitk.GetArrayFromImage(image_mask_bianca)


    view_slice_all(image_original=image_array_t2,
               mask_ant_original=image_array_mask_ant,
               mask_bianca_original=image_array_mask_bianca,
               slice_idx=120,
               alpha = 1)


    view_slice_one_slice(image_original=image_array_t2,
                   mask_ant_original=image_array_mask_ant,
                   mask_bianca_original=image_array_mask_bianca,
                   slice_idx=90,
                   type="s",
                   alpha = 1)

    plt.tight_layout()
    plt.show()

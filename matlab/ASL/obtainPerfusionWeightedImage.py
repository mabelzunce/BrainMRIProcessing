import SimpleITK as sitk

# Leer la imagen 4D
image_4d = sitk.ReadImage('/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Nifti/CP0001/asl_3d_tra_iso_3_0_highres.nii.gz')
size = image_4d.GetSize()

# Usar el primer volumen como referencia
reference_image = image_4d[:, :, :, 0]

# Inicializar imágenes para las sumas y restas utilizando la referencia de la primera imagen
sum_pares_1 = sitk.Image(reference_image.GetSize(), sitk.sitkFloat32)
sum_pares_1.CopyInformation(reference_image)

resta_impares_1 = sitk.Image(reference_image.GetSize(), sitk.sitkFloat32)
resta_impares_1.CopyInformation(reference_image)

sum_impares_2 = sitk.Image(reference_image.GetSize(), sitk.sitkFloat32)
sum_impares_2.CopyInformation(reference_image)

resta_pares_2 = sitk.Image(reference_image.GetSize(), sitk.sitkFloat32)
resta_pares_2.CopyInformation(reference_image)

# Contar los volúmenes pares e impares
num_pares = (size[3] + 1) // 2  # Número de volúmenes pares
num_impares = size[3] // 2      # Número de volúmenes impares

# Función para reamostrar imágenes
def resample_to_reference(input_image, reference_image):
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(reference_image)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetTransform(sitk.Transform())
    return resampler.Execute(input_image)

# Iterar sobre la cuarta dimensión
for i in range(size[3]):
    # Extraer la imagen en la cuarta dimensión (volumen i)
    volume = image_4d[:, :, :, i]

    # resample el volumen al espacio físico de la imagen de referencia
    resampled_volume = resample_to_reference(volume, reference_image)

    if i % 2 == 0:
        # Si es par, sumar (caso 1) restar (caso2)
        sum_pares_1 = sitk.Add(sum_pares_1, sitk.Cast(resampled_volume, sitk.sitkFloat32))
        resta_pares_2 = sitk.Subtract(resta_pares_2, sitk.Cast(resampled_volume, sitk.sitkFloat32))
    else:
        # Si es impar, restar (caso 1) rsumarr(caso2)
        resta_impares_1 = sitk.Subtract(resta_impares_1, sitk.Cast(resampled_volume, sitk.sitkFloat32))
        sum_impares_2 = sitk.Add(sum_impares_2, sitk.Cast(resampled_volume, sitk.sitkFloat32))

# Combinar el resultado de la suma y la resta
resultado_final_1 = sitk.Add(sum_pares_1, resta_impares_1)
resultado_final_2 = sitk.Add(sum_impares_2, resta_pares_2)

# Calcular el promedio
promedio_final_1 = sitk.Divide(resultado_final_1, num_pares + num_impares)
promedio_final_2 = sitk.Divide(resultado_final_2, num_pares + num_impares)

# Guardar el resultado final y el promedio
sitk.WriteImage(resultado_final_1, '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/perfusion_weighted_pares_impares.nii')
sitk.WriteImage(resultado_final_2, '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/perfusion_weighted_impares_pares.nii')

sitk.WriteImage(promedio_final_1, '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/promedio_pares_impares.nii')
sitk.WriteImage(promedio_final_2, '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/promedio_impares_pares.nii')

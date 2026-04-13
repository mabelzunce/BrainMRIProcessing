#Librerias estandar de neuroimagen
import os # genera rutas segun el sistema operativo
import glob # buscador de archivos
import numpy as np # para calculos
import nibabel as nib # convierte los archivos .nii en nros para np
import matplotlib.pyplot as plt #grafica
import pandas as pd # DataFrames
from scipy.ndimage import binary_erosion # scipy trae dif. herramientas, en este caso p/ erosionar
from scipy.stats import iqr
from scipy.io import savemat
import glob

# RUTAS
dataPartitionPath = '/data/' # raiz del disco
dataPath = '/data/belen/Proyecto_Alzheimer/derivatives/' # carpeta con todos los sujetos
preprocessedImageFilenames = 'desc-preproc_bold.nii.gz' # func preprocesada
preprocessedFolder = ''

preprocessedDataPath = os.path.join(dataPath, preprocessedFolder) # os.path.join construye rutas de archivos, por ejem0plo, subcarpetas

#Rutas de salida
outputNormalizedImageSubdir = 'Imag_NormalizadasIntensidad' # nombre de la carpeta para las imagenes nuevas
outputPathImages = os.path.join(dataPath, outputNormalizedImageSubdir) # ruta completa donde se guardan las imag NIfTI
outputPathSignals = os.path.join(dataPath, 'Results', f'ROISignals_{outputNormalizedImageSubdir}') # rutas donde se guardan los CSV/Tablas

#Crear directorios si no existen
os.makedirs(outputPathImages, exist_ok=True) # os.makedirs crea las carpetas y si ya estan creadas, continua sin error
os.makedirs(outputPathSignals, exist_ok=True)

#PARAMETROS 
# de normalizacion
saveNormalizedImages = False
scaleMeanValue = 10000 # mismo valor para todas las imagenes 

# listas para almacenar métricas
metrics = { # diccionario para guardar los resultados
    'grandMeanValue': [], 'gmMeanValue': [], 'stdMeanValue': [],
    'wmMeanValue': [], 'csfMeanValue': [],
    'grandMeanValuePost': [], 'gmMeanValuePost': [], 'stdMeanValuePost': [],
    'wmMeanValuePost': [], 'csfMeanValuePost': [],
    'tSNR': [], 'tSNR_gm': [], 'SNR': [], 'SNR_gm': [], 'CNR': []
}

# función auxiliar para cargar máscaras como booleanas
def load_mask(path): # def define una accion personalizada
    img = nib.load(path) # carga el archivo del disco
    data = img.get_fdata() # extrae los numero (de la matriz 3D)
    return img.get_fdata().astype(bool) # lo convierte a verdadero/falso a partir de astype(bool)

# lista para guardar los calculos de calidad
subject_names_processed = [] # creo lista vacia

# lista de sujetos a analizar
fmriSubjects = [d for d in os.listdir(dataPath) 
                if d.startswith('sub-') and os.path.isdir(os.path.join(dataPath, d))] # solo si es carpeta (isdir)
                 # solo busca en carpetas que empiezan con sub-
fmriSubjects.sort() # ordenar alfabéticamente los id de los sujetos

# ANALISIS POR SUJETO
for subj in fmriSubjects: # para cada sujeto de la lista...
    
    # rutas dinámicas para el sujeto actual
    print(f"--- Procesando: {subj} ---")
    
    #subj_func_dir guarda la ruta hacia la carpeta funcional => imagen BOLD + mascara del cerebro
    subj_func_dir = os.path.join(dataPath, subj, 'ses-01', 'func') 
    
    #subj_anat_dir guarda la ruta hacia la carpeta anatomica => mascaras de tejido                                                            
    subj_anat_dir = os.path.join(dataPath, subj, 'ses-01', 'anat')

    # Busqueda recursiva de archivos NIfTI y Mascaras
    try:
        bold_path = glob.glob(os.path.join(subj_func_dir, '*_res-2_desc-preproc_bold.nii.gz'))[0]

        # mascaras
        path_brain = glob.glob(os.path.join(subj_func_dir, '*_res-2_desc-brain_mask.nii.gz'))[0]
        path_gm    = glob.glob(os.path.join(subj_anat_dir, '*_res-2_label-GM_probseg.nii.gz'))[0]
        path_wm    = glob.glob(os.path.join(subj_anat_dir, '*_res-2_label-WM_probseg.nii.gz'))[0]
        path_csf   = glob.glob(os.path.join(subj_anat_dir, '*_res-2_label-CSF_probseg.nii.gz'))[0]

        # cargar datos en memoria
        fmri_img = nib.load(bold_path) # cargado de imagen
        fmri_data = fmri_img.get_fdata() # conversion a matriz 4D

        brain_mask = load_mask(path_brain)
        gm_mask    = load_mask(path_gm)
        wm_mask    = load_mask(path_wm)
        csf_mask   = load_mask(path_csf)
        
        #Mascara no-brain
        non_brain_mask = ~brain_mask # ~ invierte entre verdadero y falso
        struct_elem = np.ones((10, 10, 10)) # un cubo de 10x10x10
        non_brain_mask = binary_erosion(non_brain_mask, structure=struct_elem) # se eliminan lo que no corresponde a cerebro
                                                                                # para evitar tomar pedazos de craneo
    
        # 1) Pre Normalización
        # Calcular la media de todos los voxeles dentro de la máscara a través del tiempo
        metrics['grandMeanValue'].append(np.mean(fmri_data[brain_mask, :])) # fmri_data[brain_mask, :] toma solo los TRUE
                                                                        # np.mean() calcula el promedio
                                                                        # .append() agrega el resultado a la lista de metricas
        metrics['gmMeanValue'].append(np.mean(fmri_data[gm_mask, :])) #
        metrics['stdMeanValue'].append(np.std(fmri_data[gm_mask, :])) # desviacion global
        metrics['wmMeanValue'].append(np.mean(fmri_data[wm_mask, :]))
        metrics['csfMeanValue'].append(np.mean(fmri_data[csf_mask, :]))

        # 2) Normalizacion
        current_grand_mean = metrics['grandMeanValue'][-1]

        if current_grand_mean == 0: current_grand_mean = 1e-10 # Evita división por cero

        # todas las imagenes con el mismo valor de brillo base
        # se divide la matris fmri_data por un numero y se multiplica por 1000 (scaleMeanValue)
        fmri_data_norm = (fmri_data / current_grand_mean) * scaleMeanValue

        if saveNormalizedImages:
            output_filename_full = os.path.join(outputPathImages, subj, preprocessedImageFilenames)
            if not os.path.exists(output_filename_full):
                out_dir_subj = os.path.join(outputPathImages, subj)
                os.makedirs(out_dir_subj, exist_ok=True)
                #Crear nueva imagen Nifti
                new_img = nib.Nifti1Image(fmri_data_norm, fmri_img.affine, fmri_img.header)
                nib.save(new_img, output_filename_full)

        # 3) Metricas post-normalizacion
        metrics['grandMeanValuePost'].append(np.mean(fmri_data_norm[brain_mask, :]))
        gm_post = np.mean(fmri_data_norm[gm_mask, :])
        metrics['gmMeanValuePost'].append(gm_post)
        metrics['stdMeanValuePost'].append(np.std(fmri_data_norm[gm_mask, :]))
        wm_post = np.mean(fmri_data_norm[wm_mask, :])
        metrics['wmMeanValuePost'].append(wm_post)
        metrics['csfMeanValuePost'].append(np.mean(fmri_data_norm[csf_mask, :]))

        # 4) tSNR
        # axis=3 incluye las 4 dimensiones (0,1,2,3)

        # Mapas temporales
        mean_brain_time_courses = np.mean(fmri_data_norm, axis=3) # media (promedio) temporal == fuerza de la señal
                                                                # promedio debil == hueso/aire
                                                                # promedio alto == tejido cerebral
        std_brain_time_courses = np.std(fmri_data_norm, axis=3) # desviación estandar temporal == variabilidad

        # Evitar dividir por 0
        std_brain_time_courses[std_brain_time_courses == 0] = 1e-10 # datos nulos = 1e-10

        # tSNR (Temporal SNR) = Mean time / Std time
        # Variacion de la señal en el tiempo contemplando las 4 dimensiones (axis=3)
        tsnr_map = mean_brain_time_courses / std_brain_time_courses # Relacion señal-ruido temporal
                                                                    # buscamos tSNR alto (= señal potente y poco ruido)

        metrics['tSNR'].append(np.mean(tsnr_map[brain_mask]))
        metrics['tSNR_gm'].append(np.mean(tsnr_map[gm_mask]))

        # 5) SNR
        #Calculo de desviacion estandar del fondo (ruido espacial)
        # se usan todos los voxeles de non_brain_mask a lo largo del tiempo
        background_std = np.std(fmri_data_norm[non_brain_mask, :])  # calculo desviacion estandar

        if background_std == 0: background_std = 1e-10 # evito dividir por 0

        metrics['SNR'].append(metrics['grandMeanValuePost'][-1] / background_std) # calculo señal/ruido y guardado en la lista
        metrics['SNR_gm'].append(gm_post / background_std) # idem materia gris
        metrics['CNR'].append((gm_post - wm_post) / background_std) # CNR (Contrast to Noise Ratio)

        # se guarda el nombre solo si funciono
        subject_names_processed.append(subj)
        print(f"{subj} completado exitosamente.")

    except IndexError as e:
        print(f"Advertencia: Faltan archivos para {subj}. Error: {e}")
        continue

# Convertir listas a arrays numpy
for key in metrics:
  metrics[key] = np.array(metrics[key])

# PLOTS
# Gráfico: Valores medios POST Normalizacion
plt.figure(figsize=(10, 6))
plt.plot(metrics['grandMeanValuePost'], label='Grand Mean')
plt.plot(metrics['gmMeanValuePost'], label='Grey Matter')
plt.plot(metrics['wmMeanValuePost'], label='White Matter')
plt.plot(metrics['csfMeanValuePost'], label='CSF')
plt.title('Mean Values (Post-Normalization)')
plt.legend()
plt.show()

# Grafico: Desviacion global y outliners
plt.figure(figsize=(10, 6))
plt.plot(metrics['stdMeanValuePost'], label='Desviacion Global')

# Gráfico: SNR Metrics con Outliers
def detect_outliers_indices(data): # funcion que detectan valores alejados de la media
    median = np.nanmedian(data)
    mad = np.nanmedian(np.abs(data - median)) # calcula la distancia a la media y almacena mad
    scaled_mad = mad * 1.4826 # para consistencia con dist normal
    threshold = 3 * scaled_mad # limite para conciderar que tan lejos esta
    return np.where(np.abs(data - median) > threshold)[0]  # np.where detecta en que posicion de la lista

# almacenamiento de outliners
outliersSNR = detect_outliers_indices(metrics['SNR'])
outliersSNR_gm = detect_outliers_indices(metrics['SNR_gm'])
outlierstSNR_gm  = detect_outliers_indices(metrics['tSNR_gm'])

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# SNR
axes[0].plot(metrics['SNR'])
axes[0].plot(outliersSNR, metrics['SNR'][outliersSNR], 'x', color='red')
axes[0].set_title('SNR')

# SNR GM
axes[1].plot(metrics['SNR_gm'])
axes[1].plot(outliersSNR_gm, metrics['SNR_gm'][outliersSNR_gm], 'x', color='red')
axes[1].set_title('SNR Grey Matter')

# tSNR GM
axes[2].plot(metrics['tSNR_gm'])
axes[2].plot(outlierstSNR_gm, metrics['tSNR_gm'][outlierstSNR_gm], 'x', color='red')
axes[2].set_title('tSNR Grey Matter')

plt.tight_layout()
plt.show()

#Reporte
# Filtrar outliers: nos quedamos con los valores altos (+calidad) y descartamos los valores bajor (+ruido)
median_snr = np.nanmedian(metrics['SNR']) # calculo de la mediana ignorando los NaN
median_snr_gm = np.nanmedian(metrics['SNR_gm'])
median_tsnr_gm = np.nanmedian(metrics['tSNR_gm'])

# se crea una nueva lista: low_snr_idx
# se recorre la lista y se guardan los indices de
# solo aquellos que el SNR de ese sujeto es menor que la mediana
low_snr_idx = [i for i in outliersSNR if metrics['SNR'][i] < median_snr]
low_snr_gm_idx = [i for i in outliersSNR_gm if metrics['SNR_gm'][i] < median_snr_gm]
low_tsnr_gm_idx = [i for i in outlierstSNR_gm if metrics['tSNR_gm'][i] < median_tsnr_gm]

print("\n--- REPORTE DE OUTLIERS ---")
print(f"Sujetos con bajo SNR GM: {[subject_names_processed[i] for i in low_snr_gm_idx]}")
print(f"Sujetos con bajo tSNR GM: {[subject_names_processed[i] for i in low_tsnr_gm_idx]}")

nan_snr_gm = np.where(np.isnan(metrics['SNR_gm']))[0]
nan_tsnr_gm = np.where(np.isnan(metrics['tSNR_gm']))[0]

if len(nan_snr_gm) > 0:
    print(f"Sujetos con NaN SNR: {[subject_names_processed[i] for i in nan_snr_gm]}")
if len(nan_tsnr_gm) > 0:
    print(f"Sujetos con NaN tSNR: {[subject_names_processed[i] for i in nan_tsnr_gm]}")

# GUARDAR
# guardar tabla
df_metrics = pd.DataFrame({
    'subjectsNames': subject_names_processed,
    'grandMeanValue': metrics['grandMeanValue'],
    'grandMeanValueNorm': metrics['grandMeanValuePost'],
    'SNR': metrics['SNR'],
    'SNR_gm': metrics['SNR_gm'],
    'tSNR_gm': metrics['tSNR_gm'],
    'CNR': metrics['CNR']
})

# imprime el Dataframe por terminal
print("\n--- RESULTADOS CALCULADOS ---")
print(df_metrics.to_string(index=False))

# Guardar como CSV
csv_path = os.path.join(dataPath, 'tableImageQuality.csv') # nombre carpeta donde se guarda csv
df_metrics.to_csv(csv_path, index=False)

print(f"\nTabla guardada en: {csv_path}")

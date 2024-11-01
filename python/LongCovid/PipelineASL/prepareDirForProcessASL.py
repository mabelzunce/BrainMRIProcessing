import os, shutil
import pandas as pd

parent_dir = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti"
respuestas_cuestionarios_csv= '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/VoluntariosProyectoCovidProlongado.csv'

dir_processed = [
    'TissueVolume',
    'T2Check',
    'T1Check',
    'SliceGradientCheck',
    'SD_SNR',
    'RawSourceIMCheck',
    'MotionASL',
    'M0Reg_ASL',
    'FLAIRCheck',
    'ASLCheck',
]


directories = {
    'processed_total': "BIDS_2", # Cambiar
    'processed_long_covid': "ProcessedLC",
    'processed_control': "ProcessedControl",
    'raw_data': "rawdata",
    'derivatives': "derivatives/ExploreASL",
}


def filter_participants(participants_file):

    directory = os.path.dirname(participants_file)
    df = pd.read_csv(participants_file, sep='\t')
    existing_folders = {name for name in os.listdir(directory) if os.path.isdir(os.path.join(directory, name))}
    df_filtered = df[df['participant_id'].isin(existing_folders)]
    df_filtered.to_csv(participants_file, sep='\t', index=False)



def create_directory_structure(parent_dir, processed_name):

    processed_dir = os.path.join(parent_dir, processed_name)
    rawdata_dir = os.path.join(processed_dir, directories['raw_data'])
    derivatives_dir = os.path.join(processed_dir, directories['derivatives'])

    os.makedirs(processed_dir, exist_ok=True)
    os.makedirs(rawdata_dir, exist_ok=True)
    os.makedirs(derivatives_dir, exist_ok=True)

    return {
        'rawdata': rawdata_dir,
        'derivatives': derivatives_dir,
        'participants_derivatives': os.path.join(derivatives_dir, 'participants.tsv'),
        'dataset_description_derivatives': os.path.join(derivatives_dir, 'dataset_description.json'),
        #'participants_rawdata': os.path.join(rawdata_dir, 'participants.tsv'),
        'dataset_description_rawdata': os.path.join(rawdata_dir, 'dataset_description.json'),
    }

def copy_if_not_exists(source, destination):
    """Copia el archivo de participantes si no existe."""
    if not os.path.exists(destination):
        shutil.copy(source, destination)

def copy_files_with_sequence(source_folder, destination_folder, sequence):
    """
    Copia archivos de source_folder a destination_folder si contienen la sequence en su nombre.

    :param source_folder: Ruta de la carpeta fuente
    :param destination_folder: Ruta de la carpeta de destino
    :param sequence: Secuencia de strings a buscar en el nombre del archivo
    """
    os.makedirs(destination_folder, exist_ok=True)



    # Listar todos los archivos en la carpeta fuente
    for filename in os.listdir(source_folder):
        if os.path.isdir(os.path.join(source_folder, filename)):
            continue

        # Verificar si la secuencia está en el nombre del archivo
        if sequence in str(filename):
            source_path = os.path.join(source_folder, filename)
            destination_path = os.path.join(destination_folder, filename)

            # Copiar el archivo
            shutil.copy(source_path, destination_path)


if __name__ == '__main__':
    # Crear directorios procesados
    total_processed = create_directory_structure(parent_dir, directories['processed_total'])
    long_covid_processed = create_directory_structure(parent_dir, directories['processed_long_covid'])
    control_processed = create_directory_structure(parent_dir, directories['processed_control'])

    # Copiar archivos de participantes
    # copy_if_not_exists(total_processed['participants_derivatives'], control_processed['participants_derivatives'])
    # copy_if_not_exists(total_processed['participants_rawdata'], control_processed['participants_rawdata'])
    # copy_if_not_exists(total_processed['participants_derivatives'], long_covid_processed['participants_derivatives'])
    # copy_if_not_exists(total_processed['participants_rawdata'], long_covid_processed['participants_rawdata'])

    copy_if_not_exists(total_processed['dataset_description_derivatives'], control_processed['dataset_description_derivatives'])
    copy_if_not_exists(total_processed['dataset_description_rawdata'], control_processed['dataset_description_rawdata'])
    copy_if_not_exists(total_processed['dataset_description_derivatives'],
                       long_covid_processed['dataset_description_derivatives'])
    copy_if_not_exists(total_processed['dataset_description_rawdata'],
                       long_covid_processed['dataset_description_rawdata'])

    # Create aux dir
    for dir in dir_processed:
        control_dir = os.path.join(control_processed['derivatives'], 'Population',dir)
        covid_dir = os.path.join(long_covid_processed['derivatives'], 'Population',dir)
        os.makedirs(control_dir, exist_ok=True)
        os.makedirs(covid_dir, exist_ok=True)

    df_respuestas_cuestionarios_csv = pd.read_csv(respuestas_cuestionarios_csv)

    subjects_to_process = [subject for subject in sorted(os.listdir(total_processed['rawdata']))
                           if subject.startswith("sub")]
    subjects_to_exclude = ['sub-CP0075', 'sub-CP0076',
                           'sub-CP0062', 'sub-CP0078','sub-CP0079',
                            'sub-CP0105', 'sub-CP0216']

    for subject in subjects_to_process:
        if subject in subjects_to_exclude:
            continue

        total_subject_raw_data = os.path.join(total_processed['rawdata'], subject)
        total_subject_derivatives = os.path.join(total_processed['derivatives'], f'{subject}_1')


        subject_id = subject[4:]
        subject_group = df_respuestas_cuestionarios_csv[df_respuestas_cuestionarios_csv['ID'] == subject_id]['Grupo']

        if not subject_group.empty:
            subject_group_name = subject_group.iloc[0]
        else:
            print(f"No se encontró grupo para el sujeto con ID: {subject_id}")
            continue

        if subject_group_name == "Control":
            subject_control_raw_data = os.path.join(control_processed['rawdata'], subject)
            subject_control_derivatives = os.path.join(control_processed['derivatives'], f'{subject}_1')

            if not os.path.exists(subject_control_raw_data):
                shutil.copytree(total_subject_raw_data, subject_control_raw_data)

            if not os.path.exists(subject_control_derivatives):
                shutil.copytree(total_subject_derivatives, subject_control_derivatives)

            copy_files_with_sequence(os.path.join(total_processed['derivatives'], 'Population'),
                                     os.path.join(control_processed['derivatives'], 'Population'),
                                     subject[4:])

            for dir in dir_processed:
                if os.path.exists(os.path.join(total_processed['derivatives'], 'Population', dir)):
                    copy_files_with_sequence(os.path.join(total_processed['derivatives'], 'Population', dir),
                                             os.path.join(control_processed['derivatives'], 'Population', dir), subject[4:])

        elif subject_group_name == "Covid Prolongado":
            subject_lc_raw_data = os.path.join(long_covid_processed['rawdata'], subject)
            subject_lc_derivatives = os.path.join(long_covid_processed['derivatives'], f'{subject}_1')

            if not os.path.exists(subject_lc_raw_data):
                shutil.copytree(total_subject_raw_data, subject_lc_raw_data)

            if not os.path.exists(subject_lc_derivatives):
                shutil.copytree(total_subject_derivatives, subject_lc_derivatives)

            copy_files_with_sequence(os.path.join(total_processed['derivatives'], 'Population'),
                                     os.path.join(long_covid_processed['derivatives'], 'Population'),
                                     subject[4:])
            for dir in dir_processed:
                if os.path.exists(os.path.join(total_processed['derivatives'], 'Population', dir)):
                    copy_files_with_sequence(os.path.join(total_processed['derivatives'], 'Population', dir),
                                             os.path.join(long_covid_processed['derivatives'], 'Population', dir), subject[4:])


        else:
            print(subject_id, subject_group)

    # filter_participants(os.path.join(control_processed['derivatives'],'participants.tsv'))
    # filter_participants(os.path.join(control_processed['rawdata'],'participants.tsv'))
    # filter_participants(os.path.join(long_covid_processed['derivatives'],'participants.tsv'))
    # filter_participants(os.path.join(long_covid_processed['rawdata'],'participants.tsv'))

ruta_principal = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/BIDS_MRI/'; % Ruta de la base de datos o de prueba
carpeta_dicoms = '/media/delfi/Lexar/Estudio2/MRI/';
cd(ruta_principal);

carpetas = dir(carpeta_dicoms);
carpetas = carpetas([carpetas.isdir] & ~strcmp({carpetas.name}, '.') & ~strcmp({carpetas.name}, '..'));

for i = 1:length(carpetas)
%for i = 1:1
    nombre_sujeto = carpetas(i).name;
    fprintf('Procesando: %s\n', nombre_sujeto);
    comando = sprintf('bash -c "source /home/delfi/anaconda3/bin/activate dcm2bids && dcm2bids -d %s%s -p %s -c %scode/dcm2bids_config.json -o /mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/BIDS_MRI/"', carpeta_dicoms, nombre_sujeto, nombre_sujeto, ruta_principal);
    system(comando);
    
    % Copiar y renombrar aslcontext.tsv a la carpeta perf del sujeto
    origen = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/bids_ejemplo_MRI/code/_aslcontext.tsv';
    
    nombre_archivo_destino = sprintf('sub-%s_aslcontext.tsv', nombre_sujeto);
    destino = fullfile(['sub-' nombre_sujeto], 'perf', nombre_archivo_destino);
    
    if isfolder(fullfile(['sub-' nombre_sujeto], 'perf'))
        copyfile(origen, destino);
        fprintf(' - %s copiado para %s\n', nombre_archivo_destino, nombre_sujeto);
    else
        fprintf(' - Carpeta perf no encontrada para %s\n', nombre_sujeto);
    end
end

fprintf('Proceso completado\n');

ruta_principal = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/BIDS_PET/';
carpeta_dicoms = '/media/delfi/Lexar/Estudio2/PET/';

cd(ruta_principal);
carpetas = dir(carpeta_dicoms);
carpetas = carpetas([carpetas.isdir] & ~strcmp({carpetas.name}, '.') & ~strcmp({carpetas.name}, '..'));

for i = 1:length(carpetas)
%for i = 1:1
    nombre_sujeto = carpetas(i).name;
    fprintf('Procesando: %s\n', nombre_sujeto);
    comando = sprintf('bash -c "source /home/delfi/anaconda3/bin/activate dcm2bids && dcm2bids -d %s%s -p %s -c %scode/dcm2bids_config_PET.json -o /mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/BIDS_PET/"', carpeta_dicoms, nombre_sujeto, nombre_sujeto, ruta_principal);
    system(comando);
end

fprintf('Proceso completado\n');

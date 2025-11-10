% Este codigo busca los archivos PerfusionWeighted, ADC y TRACE de niftis
% sin convertir y los guarda en la carpeta derivatives
ruta_origen = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/BIDS_MRI/tmp_dcm2bids/'; %ruta con los niftis
ruta_destino_base = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/BIDS_MRI/derivatives/';


carpetas = dir(fullfile(ruta_origen, 'sub-CP*'));
carpetas = carpetas([carpetas.isdir]);

for i = 1:length(carpetas)
    nombre_sujeto = carpetas(i).name;
    carpeta_sujeto = fullfile(ruta_origen, nombre_sujeto);
    
    fprintf('Procesando: %s\n', nombre_sujeto);
    
    archivos_tensor = dir(fullfile(carpeta_sujeto, '*TENSOR*.nii.gz'));
    archivos_movidos = {};
    
    archivos_asl = dir(fullfile(carpeta_sujeto, '*asl*.nii.gz'));
    
    for j = 1:length(archivos_tensor)
        archivo_nii = archivos_tensor(j).name;
        archivo_json = strrep(archivo_nii, '.nii.gz', '.json');
        
        ruta_nii = fullfile(carpeta_sujeto, archivo_nii);
        ruta_json = fullfile(carpeta_sujeto, archivo_json);
        
        if ~isfile(ruta_json)
            continue;
        end
    
        tipo_archivo = '';
        
        try
            json_data = jsondecode(fileread(ruta_json));
            json_string = fileread(ruta_json);
            
            if isfield(json_data, 'ImageType') && iscell(json_data.ImageType)
                target_tracew = {"DERIVED", "PRIMARY", "DIFFUSION", "TRACEW", "ND", "NORM"};
                target_adc = {"DERIVED", "PRIMARY", "DIFFUSION", "ADC", "ND", "NORM"};
                
                if isequal(json_data.ImageType, target_tracew) || contains(json_string, 'TRACEW')
                    tipo_archivo = 'trace';
                elseif isequal(json_data.ImageType, target_adc) || contains(json_string, 'ADC')
                    tipo_archivo = 'adc';
                end
            else
                if contains(json_string, 'TRACEW')
                    tipo_archivo = 'trace';
                elseif contains(json_string, 'ADC')
                    tipo_archivo = 'adc';
                end
            end
            
        catch
            continue;
        end
        
        if ~isempty(tipo_archivo)
            % nombres nuevos
            nuevo_nombre_nii = sprintf('%s_acq-%s_dwi.nii.gz', nombre_sujeto, tipo_archivo);
            nuevo_nombre_json = sprintf('%s_acq-%s_dwi.json', nombre_sujeto, tipo_archivo);
            
            carpeta_destino = fullfile(ruta_destino_base, nombre_sujeto, 'dwi');
            if ~isfolder(carpeta_destino)
                mkdir(carpeta_destino);
            end
            
            destino_nii = fullfile(carpeta_destino, nuevo_nombre_nii);
            destino_json = fullfile(carpeta_destino, nuevo_nombre_json);

            try
                movefile(ruta_nii, destino_nii);
                movefile(ruta_json, destino_json);
                
                base_name = strrep(archivo_nii, '.nii.gz', '');
                archivo_bvec = [base_name '.bvec'];
                archivo_bval = [base_name '.bval'];
                
                ruta_bvec = fullfile(carpeta_sujeto, archivo_bvec);
                ruta_bval = fullfile(carpeta_sujeto, archivo_bval);
                
                nuevo_base = sprintf('%s_acq-%s_dwi', nombre_sujeto, tipo_archivo);
                destino_bvec = fullfile(carpeta_destino, [nuevo_base '.bvec']);
                destino_bval = fullfile(carpeta_destino, [nuevo_base '.bval']);
                
                if isfile(ruta_bvec)
                    movefile(ruta_bvec, destino_bvec);
                    archivos_movidos{end+1} = [nuevo_base '.bvec'];
                end
                
                if isfile(ruta_bval)
                    movefile(ruta_bval, destino_bval);
                    archivos_movidos{end+1} = [nuevo_base '.bval'];
                end
                
                archivos_movidos{end+1} = nuevo_nombre_nii;
                archivos_movidos{end+1} = nuevo_nombre_json;
                
            catch ME
            end
        end
    end
    
    
    for j = 1:length(archivos_asl)
        archivo_nii = archivos_asl(j).name;
        archivo_json = strrep(archivo_nii, '.nii.gz', '.json');
        
        ruta_nii = fullfile(carpeta_sujeto, archivo_nii);
        ruta_json = fullfile(carpeta_sujeto, archivo_json);
        
        if ~isfile(ruta_json)
            continue;
        end
        
        try
            json_data = jsondecode(fileread(ruta_json));
            
            if isfield(json_data, 'SeriesDescription') && strcmp(json_data.SeriesDescription, 'Perfusion_Weighted')
                nuevo_nombre_nii = sprintf('%s_PerfusionWeighted.nii.gz', nombre_sujeto);
                nuevo_nombre_json = sprintf('%s_PerfusionWeighted.json', nombre_sujeto);
                
                carpeta_destino = fullfile(ruta_destino_base, nombre_sujeto, 'perf');
                if ~isfolder(carpeta_destino)
                    mkdir(carpeta_destino);
                end
                
                destino_nii = fullfile(carpeta_destino, nuevo_nombre_nii);
                destino_json = fullfile(carpeta_destino, nuevo_nombre_json);
                
                try
                    movefile(ruta_nii, destino_nii);
                    movefile(ruta_json, destino_json);
                    
                    archivos_movidos{end+1} = nuevo_nombre_nii;
                    archivos_movidos{end+1} = nuevo_nombre_json;
                    
                catch ME
                end
            end
            
        catch
            continue;
        end
    end
    
    if ~isempty(archivos_movidos)
        fprintf('  Archivos movidos: %s\n', strjoin(archivos_movidos, ', '));
    end
    
end

fprintf('\nProceso completado.\n');

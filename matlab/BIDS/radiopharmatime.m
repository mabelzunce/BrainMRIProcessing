% Este código hace una tabla con RadiopharmaceuticalStartTime,
% AcquisitionTime y la diferencia en minutos y segundos entre ambos a
% partir de los archivos DICOM de PET
ruta_base = '/media/delfi/Lexar/Estudio2/PET/'; %carpetas DICOM PET
ruta_guard = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/BIDS_PET/';  
cd(ruta_guard);
tabla_resultados = table();

sujetos = dir(fullfile(ruta_base, 'CP*'));
sujetos = sujetos([sujetos.isdir]);

for i = 1:length(sujetos)
    nombre_sujeto = sujetos(i).name;
    ruta_sujeto = fullfile(ruta_base, nombre_sujeto, 'DICOM');
    nombre_sujeto_BIDS = strcat('sub-', nombre_sujeto);
    
    if ~exist(ruta_sujeto, 'dir')
        continue;
    end
    
    carpetas_finales = encontrar_carpetas_finales(ruta_sujeto);
    
    for j = 1:length(carpetas_finales)
        ruta_carpeta_final = carpetas_finales{j};
        subcarpetas = dir(ruta_carpeta_final);
        subcarpetas = subcarpetas([subcarpetas.isdir] & ~strcmp({subcarpetas.name}, '.') & ~strcmp({subcarpetas.name}, '..'));
        
        for k = 1:length(subcarpetas)
            ruta_subcarpeta = fullfile(ruta_carpeta_final, subcarpetas(k).name);
            archivos = dir(ruta_subcarpeta);
            archivos = archivos(~[archivos.isdir]);
            
            if isempty(archivos)
                continue;
            end
            
            archivo_dicom = fullfile(ruta_subcarpeta, archivos(1).name);
            
            try
                info_dicom = dicominfo(archivo_dicom);
                
                if isfield(info_dicom, 'SeriesDescription')
                    series_desc = info_dicom.SeriesDescription;
                    
                    if contains(series_desc, 'BR_CTAC')
                        radio_start_time = '';
                        acquisition_time = '';
                        
                        if isfield(info_dicom, 'RadiopharmaceuticalInformationSequence') && ...
                           isfield(info_dicom.RadiopharmaceuticalInformationSequence, 'Item_1') && ...
                           isfield(info_dicom.RadiopharmaceuticalInformationSequence.Item_1, 'RadiopharmaceuticalStartTime')
                            radio_start_time = info_dicom.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;
                        end

                       if isfield(info_dicom, 'AcquisitionTime')
                            acquisition_time = info_dicom.AcquisitionTime;
                        end
                        
                        nueva_fila = table({nombre_sujeto_BIDS}, {series_desc}, {radio_start_time}, {acquisition_time}, ...
                            'VariableNames', {'Sujeto', 'SeriesDescription', 'RadiopharmaceuticalStartTime', 'AcquisitionTime'});
                        tabla_resultados = [tabla_resultados; nueva_fila];
                    end
                end
            catch
            end
        end
    end
end

diferencias_minutos = zeros(height(tabla_resultados), 1);

for i = 1:height(tabla_resultados)
    radio_time = tabla_resultados.RadiopharmaceuticalStartTime(i);
    acq_time = tabla_resultados.AcquisitionTime(i);
 
    radio_minutos = tiempo_a_minutos(radio_time);
    acq_minutos = tiempo_a_minutos(acq_time);
    
    if ~isnan(radio_minutos) && ~isnan(acq_minutos)
        diferencia = acq_minutos - radio_minutos;
        
        diferencias_minutos(i) = diferencia;
    else
        diferencias_minutos(i) = NaN; 
    end
end

tabla_resultados.DiferenciaMinutos = diferencias_minutos;
tabla_resultados.DiferenciaSegundos = diferencias_minutos*60*(-1);

disp('Tabla con diferencias de tiempo:');
disp(tabla_resultados);

writetable(tabla_resultados, '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/BIDS_MRI/code/radiotime_ESTUDIO2.csv'); %Ruta guardado tabla

function carpetas_finales = encontrar_carpetas_finales(ruta)
    carpetas_finales = {};
    contenido = dir(ruta);
    carpetas = contenido([contenido.isdir] & ~strcmp({contenido.name}, '.') & ~strcmp({contenido.name}, '..'));
    
    if length(carpetas) > 1  
        carpetas_finales{end+1} = ruta;
    else
        for i = 1:length(carpetas)
            sub_carpetas = encontrar_carpetas_finales(fullfile(ruta, carpetas(i).name));
            carpetas_finales = [carpetas_finales, sub_carpetas];
        end
    end
end

function minutos = tiempo_a_minutos(tiempo_str)
    if isempty(tiempo_str) || isempty(tiempo_str{1})
        minutos = NaN;
        return;
    end
    
    tiempo = tiempo_str{1}; 
    tiempo = sprintf('%06s', tiempo);
    horas = str2double(tiempo(1:2));
    mins = str2double(tiempo(3:4));
    segs = str2double(tiempo(5:6));
    minutos = horas * 60 + mins + segs / 60;
end

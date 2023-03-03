%% TODO EL PROCESO JUNTO
% compilar solo las secciones que hagan falta

clc
close all 
clear all

%% Declaro variables utiles

tandas_AD = {'AD_tanda1','AD_tanda2','AD_tanda3','AD_tanda4','AD_tanda5','AD_tanda6','AD_tanda7'} ;
tandas_CN = {'CN_tanda1','CN_tanda2','CN_tanda3','CN_tanda4','CN_tanda5','CN_tanda6'} ;

headers = '/dcmHeaders.mat';
path = '/home/triniibar/Desktop/' ;
% todos estos dir ya estan creados
carepta_SIEMENS = strcat(path, 'SIEMENS') ;
AD = strcat ( carepta_SIEMENS , '/AD');
CN = strcat ( carepta_SIEMENS , '/CN');
disp('variables inicializadas')


%% [0] AGREGAR PATHS 

clc
cd /home/triniibar/
addpath spm12
cd /home/triniibar/
addpath dicm2
cd /home/triniibar/
addpath /home/triniibar/DPABI_V6.2_220915
cd DPABI_V6.2_220915/
addpath Analysis
addpath DPABINet
addpath DPABISurf
addpath DPARSF
addpath Logo
addpath QualityControl
addpath Standardization
addpath StatisticalAnalysis
addpath Subfunctions
addpath Templates
addpath Utilities
addpath Viewer
cd DPARSF
addpath Docs
addpath Jobmats
addpath SubGUIs
addpath Subfunctions
addpath bet
cd /home/triniibar/
addpath Documents

disp('paths agregados');


%% [1] CONVERTIR IMAGENES DICOM A NII
% reconvertir todos los FunRaw a FunImg
% porque despues de procesar las imagenes, DPARFS sobre-escribe FunImg
% con imagenes con los primeros 5 time points cortados 
% tengo que recuperar las imagenes nii con todos los time points 
% asi que las tengo que convertir devuelta 

path = '/home/triniibar/Desktop/' ;

% AD
% loop sin comentar
% exactamente igual que los CN
for i = 1: length(tandas_AD) 
    tanda = string( tandas_AD(i) ) ;
    disp ('-------------------');
    disp(tanda);
    FilesFun = dir( strcat(path,tanda,'/FunImg') );
    for sub=1:length(FilesFun)
        ext = string(FilesFun(sub).name) ; 
        if (ext~= '.') && (ext~= '..')
            path_DICOM = strcat(path, tanda, '/FunRaw/', string(ext));
            path_nii = strcat(path, tanda, '/FunImg/', string(ext) );    
            varargout = dicm2nii(path_DICOM, path_nii, 0) ;    
        end   
    end
end

% CN
% loop comentado y explicado paso por paso
for i = 1: length(tandas_CN) 
    tanda = string( tandas_CN(i) ) ;
    disp ('-------------------');
    disp(tanda);
    FilesFun = dir( strcat(path,tanda,'/FunImg') );
    % files con las imagenes funcionales de todos los subjects
   
    for sub=1:length(FilesFun)
        ext = string(FilesFun(sub).name) ; % nombre del archivo
        if (ext~= '.') && (ext~= '..') % filtro archivos basura
            % decalro el path de la carpeta con los archivos DICOM 
            path_DICOM = strcat(path, tanda, '/FunRaw/', string(ext));
            
            % declaro el path donde se van a guardar las iamgebes nii
            path_nii = strcat(path, tanda, '/FunImg/', string(ext) );  
            
            % conversion DICOM a nii con funcion dicm2nii
            % el ultimo numero deterima el tipo de archivo (0 es .nii)
            varargout = dicm2nii(path_DICOM, path_nii, 0) ;     
        end   
    end
end


%% [2] CLASIFICAR SEGUN ESCANER
% el slice order es impares-pares para GE y Philips 
% el slice order es pares-impares para Siemens
% tengo que separar las Siemens, para procesarlas con orden par-impar


% PARA LAS IMAGENES AD
% loop sin comentar
% exactamente igual que los CN
for i = 1: length(tandas_AD)
    tanda = string( tandas_AD(i) ) ;
    FilesFun = dir( strcat(path,tanda,'/FunImg') );
    disp (tanda) ;
    
    for sub=1:length(FilesFun)
        ext = string(FilesFun(sub).name) ; 
        if (ext~= '.') && (ext~= '..') 
            subID = string(ext) ;
            subject = strcat(path, tanda, '/FunImg/', subID );
            data = load( strcat(subject , headers) );
            nombre_archivo = fieldnames(data.h) ;
            nombre_archivo = string(nombre_archivo);
            metadata = getfield(data.h, nombre_archivo);
            escaner = metadata.Manufacturer ;
            escaner = string (escaner);
            disp(escaner)
            
            if (escaner == 'SIEMENS')
                cd (AD) ;
                carpeta_tanda = strcat (AD , '/',tanda,'/FunImg') ;
                carpeta_subject = strcat (carpeta_tanda,'/', subID) ;
                if (isfolder(tanda)==0)
                    mkdir (tanda) ;
                end
                mkdir( carpeta_tanda, subID) ;
                copyfile ( subject , carpeta_subject) ;
            end
        end   
    end 
    disp ('-------------------');
end


% PARA LAS IMAGENES CN
% loop comentado y explicado paso por paso
for i = 1: length(tandas_CN)
    tanda = string( tandas_CN(i) ) ;
    % abro el archivo con las imagenes fMRI ORIGINALES
    FilesFun = dir( strcat(path,tanda,'/FunImg') );
    disp (tanda) ;
    
    for sub=1:length(FilesFun)
        ext = string(FilesFun(sub).name) ; % nombre del archivo
        if (ext~= '.') && (ext~= '..') % filtro archivos basura
            subID = string(ext) ;
            % path de la carpeta ORIGINAL del subject con su fMRI y headers
            subject = strcat(path, tanda, '/FunImg/', subID );
        
            % cargo la data guardada en los headers
            data = load( strcat(subject , headers) );
            nombre_archivo = fieldnames(data.h) ;
            nombre_archivo = string(nombre_archivo);
            metadata = getfield(data.h, nombre_archivo);
            
            % consigo el tipo de escaner guardado en la metadata
            escaner = metadata.Manufacturer ;
            escaner = string (escaner);
            disp(escaner)
            
            % si es SIEMENS, separo esos subjects en otra carpeta
            % estructura carpetas: 
            % SIEMENS--> CN--> tanda--> FunImg--> subjectID--> fMRI.nii
            
            if (escaner == 'SIEMENS')
                cd (CN) ;
                % creo las carpetas necesarias
                carpeta_tanda = strcat (CN , '/',tanda, '/FunImg') ;
                carpeta_subject = strcat (carpeta_tanda,'/', subID) ;
                % control 
                if (isfolder(tanda)==0)
                    mkdir (tanda) ;
                end
                % creo una carpeta para ese subject y guardo su fMRI ahi
                mkdir( carpeta_tanda, subID) ;
                copyfile ( subject , carpeta_subject) ;              
            end
           
        end   
    end
    disp ('-------------------');
end

disp ('-------------------');
disp ('-------------------');
disp ('-------------------');
disp('SE CREARON TODAS LAS CARPETAS');
disp ('-------------------');
disp ('-------------------');
disp ('-------------------');


%% [3] PROCESAR EN DPARFS

% reproceso una por una las tandas SIEMENS
% cada una tiene un total de TimePoints distintos
% asi que no puedo procesarlas juntas

tanda = 'CN_tanda6' ; % lo voy cambiando a mano y compilo esta seccion para cada tanda
working_dir = strcat(CN,'/',tanda) ;
cd (working_dir) ; % asi DPARFS ya se abre en el working dir y no lo tengo que cambiar a mano

DPARSF


%% [4] EXTRAER SEÑALES BOLD

% lo hice a mano con una opcion en DPABI--> Utilities --> DPABI_ROISignalExtracter

% 0) agrego los paths de DPABI

% 1) creo un directorio donde guardar las ROI Signals 
%   --> yo los llame 'output_AplicarAtlas_AD' y '_CN'

% 2) tipeo DPABI_ROISignalExtracter en la Comand Window

% 3) aparece un display, para pasarle parametros al ROISignalExtracter
%   --> 'define ROI'--> 'sphere' le paso las coordenadas DSI_enhanced.txt (MNI)
%   --> 'add image' agrego la imagen nifi (de la carpeta FunImg ARWSD)
%       (no encontre como hacer el proceso para muchas imagenes al mismo
%       tiempo, siempre me daba error, asi que lo hice una por una)
%   --> 'ouput dir' completo .../output_AplicarAtlas_AD

% 4) en output_AplicarAtlas quedan, por cada subject, estos 7 archivos
%   --> ROISignals_ROISignal (.txt) y (.mat)
%   --> ROICorrelation_FisherZ (.txt) y (.mat)
%   --> ROI_OrderKey_ROISignal (.tsv)
%   --> ROICorrelation_FisherZ_ROISignal (.txt) y (.mat)

% 5) muevo a otra carpeta los 'ROISignals_ROISignal' de cada subject
%   --> yo las llame 'ROISignals_AD' y '_CN'



% [ SCRIPT PARA MOVER LOS ARCHIVOS ] (paso 5)

% [ CN ] comentado
path_output ='/home/triniibar/Desktop/Corregidos/998ROI/CN/output_AplicarAtlas_CN/';
cd (path_output);
carpeta= '/home/triniibar/Desktop/Corregidos/998ROI/CN/';
destination = strcat(carpeta,'ROISignals_CN'); % aca van a ir a parar los archivos

files = dir ;
% empiezo en 3 porque los primeros dos archivos son '.' y '..'
for i= 3: length(files)
    nombre = string(files(i).name);
    
    % solo muevo los archivos ROISignals_ROISignal.txt
    if ( contains(nombre,'ROISignals_ROISignal')==1) && (contains(nombre,'.txt')==1)
        source=strcat(carpeta,'output_AplicarAtlas_CN/',nombre);
        copyfile(source,destination);
    end  
    
    % solo muevo los archivos ROISignals_ROISignal.mat
    if ( contains(nombre,'ROISignals_ROISignal')==1) && (contains(nombre,'.mat')==1)
        source=strcat(carpeta,'output_AplicarAtlas_CN/',nombre);
        copyfile(source,destination);
    end
    
end

% [ AD ] (lo mismo que CN)
clear path_output
clear carpeta
clear source
clear destination
path_output ='/home/triniibar/Desktop/Corregidos/998ROI/AD/output_AplicarAtlas_AD/';
cd (path_output);
carpeta= '/home/triniibar/Desktop/Corregidos/998ROI/AD/';
destination = strcat(carpeta,'ROISignals_AD'); % aca van a ir a parar los archivos
files = dir ;
for i= 3: length(files)
    nombre = string(files(i).name);
    if ( contains(nombre,'ROISignals_ROISignal')==1) && (contains(nombre,'.txt')==1)
        source=strcat(carpeta,'output_AplicarAtlas_AD/',nombre);
        copyfile(source,destination);
    end  
    if ( contains(nombre,'ROISignals_ROISignal')==1) && (contains(nombre,'.mat')==1)
        source=strcat(carpeta,'output_AplicarAtlas_AD/',nombre);
        copyfile(source,destination);
    end
end



%% [5] CONCATENAR SEÑALES BOLD 

clc 
clear all 
close all 

% CN (comentado)

carpeta= '/home/triniibar/Desktop/Corregidos/998ROI/CN/';
path_ROISignals = strcat(carpeta,'ROISignals_CN');
archivos = dir (path_ROISignals);
ROISignals_CN = [] *998;

% empiezo en 3 porque los primeros dos archivos son '.' y '..'
for i=3:length(archivos)
    nombre = string(archivos(i).name);
    
    if (contains(nombre,'mat'))
        % extraigo las señalers BOLD de un subject
        Signals = load (strcat(path_ROISignals,'/',nombre)) ; % struct
        Signals = Signals.ROISignals ; % variable dentro del struct
        
        % Signals es una matriz de dimension TimePoints x 998
        % entonces se tiene que concatenar en vertical 
        ROISignals_CN = vertcat(ROISignals_CN,Signals);
        Signals=[]*[] ; 
        % vacio la variable al final para que no haya problemas al
        % concatenarla en el siguiente ciclo
    end
end

% necesito una matriz de dimension 998 x TimePoints 
ROISignals_CN = reshape (ROISignals_CN,998,[]);
% la normalizo para que "se peguen bien las señales"
% que coincidan los principios y finales de las señales de subjects consecutivos
normROISignals_CN = (ROISignals_CN - mean(ROISignals_CN)) ./std(ROISignals_CN) ;

% guardo las variables en archivos separados
cd /home/triniibar/Desktop/Corregidos/998ROI/
save('ROISignals_CN.mat','ROISignals_CN')
save('normROISignals_CN.mat','normROISignals_CN')



% AD (lo mismo que CN)
clear carpeta
clear Signals
clear archivos
carpeta= '/home/triniibar/Desktop/Corregidos/998ROI/AD/';
path_ROISignals = strcat(carpeta,'ROISignals_AD');
archivos = dir (path_ROISignals);
ROISignals_AD = [] *998;

for i=3:length(archivos)
    nombre = string(archivos(i).name);
    if (contains(nombre,'mat'))
        Signals = load (strcat(path_ROISignals,'/',nombre)) ;
        Signals = Signals.ROISignals ; 
        ROISignals_AD = vertcat(ROISignals_AD,Signals);
        Signals=[]*[] ; 
    end
end

ROISignals_AD = reshape (ROISignals_AD,998,[]);
normROISignals_AD = (ROISignals_AD - mean(ROISignals_AD)) ./std(ROISignals_AD) ;
cd /home/triniibar/Desktop/Corregidos/998ROI/
save('ROISignals_AD.mat','ROISignals_AD')
save('normROISignals_AD.mat','normROISignals_AD')



%% [6] CALCULAR ACF

% este codigo parte de los archivos guardados en el paso anterior
% no usa ninguna de las variables anteriores, solo archivos
clear all 
close all 
clc


% cargo los archivos
cd /home/triniibar/Desktop/Corregidos/998ROI/ % directorio donde estan los archivos
AD_signals = load('normROISignals_AD.mat'); % struct
AD_signals = AD_signals.normROISignals_AD ; % variable en el struct 
CN_signals = load('normROISignals_CN.mat'); % struct
CN_signals = CN_signals.normROISignals_CN ; % variable en el struct 


% calculo de la funcion autocorrelacion 
lags= 50 ;
roi= randi(998) ; % elijo un roi aleatorio
acf_AD= autocorr(AD_signals(roi,:),'NumLags',lags) ;
acf_CN= autocorr(CN_signals(roi,:),'NumLags',lags) ;


% plot de ambas acf
figure(1)
subplot(2,1,1)
plot(acf_AD,'r');
xlim([0,50]);
ylim([-0.5,1]);
texto = strcat('señales AD','--->','roi:',string(roi));
title(texto) ;
texto = strcat('acf','--> lags ', string(lags));
legend(texto) ;
subplot(2,1,2)
plot(acf_CN,'m');
xlim([0,50]);
ylim([-0.5,1]);
texto = strcat('señales CN','--->','roi:',string(roi));
title(texto) ;
texto = strcat('acf','--> lags ', string(lags));
legend(texto) ;
hold off


% invertir las señales
AD_ROI_flip = flip(AD_signals);
CN_ROI_flip = flip(CN_signals);
acf_flipAD= autocorr(AD_ROI_flip(roi,:),'NumLags' , lags) ;
acf_flipCN= autocorr(CN_ROI_flip(roi,:),'NumLags' , lags) ;


% plot acf señal original vs invertida 
figure(2)
subplot(2,1,1)
plot(acf_AD,'r')
hold on ;
plot(acf_flipAD)
xlim([0,50]);
ylim([-0.5,1]);
texto = strcat('señales AD','--->','roi:',string(roi));
title(texto) ;
legend('acf señal original','acf señal invertida');
hold off
subplot(2,1,2)
plot(acf_CN,'m')
hold on 
plot(acf_flipCN)
xlim([0,50]);
ylim([-0.5,1]);
texto = strcat('señales CN','--->','roi:',string(roi));
title(texto) ;
legend('acf señal original','acf señal invertida');
hold off


% plot de la acf para random roi 
figure(3)
roi_i=[] ;
for i=1:5
    roi= randi(998) ;
    acf_AD= autocorr(AD_signals(roi,:),'NumLags' , lags) ;
    acf_CN= autocorr(CN_signals(roi,:),'NumLags' , lags) ;
    
    subplot(2,1,1)
    plot(acf_AD) ; 
    xlim([0,50]);
    ylim([-0.5,1]);
    title ('acf para 5 ROI aleatorios (AD)');
    hold all ; 
    
    subplot(2,1,2)
    plot(acf_CN) ; 
    xlim([0,50]);
    ylim([-0.5,1]);
    title ('acf para 5 ROI aleatorios (CN)');
    hold all ; 
    
    roi_i= [roi_i, roi];
end 
legend( string(roi_i(1)) , string(roi_i(2)) ,string(roi_i(3)) ,string(roi_i(4)),string(roi_i(5)) );







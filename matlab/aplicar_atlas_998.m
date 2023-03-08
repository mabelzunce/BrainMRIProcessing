%% APLICAR ATLAS 
clear all 
close all 
clc

%% TYPE OF PROCESSING
firstCovThenStats = 1; % if 1, first the COV matrices are computed and then the mean matrices. 
%0, mean signals and the CC matrices.


%% carga la informacion del atlas

cd /home/triniibar/Documents/codigos_para_armar_la_database/
ATLAS = load('DSI_enhanced.mat');
path = '/home/triniibar/Desktop/Corregidos/998ROI/matrices/' ;

%% Declaro las variables utiles

ROIlabels = 1:998 ;

% las coords x y z de cada ROI, puestas una atras de la otra
coordsROI = ATLAS.roi_xyz_avg ;

% las coords x y z de cada ROI, puestas en columnas
ROIcentroids = transpose(coordsROI) ;


%% Distancia entre ROI

% defino la matriz de los centroides
M_ROIcentroids = reshape(ROIcentroids, size(ROIcentroids,1), 1, []);
M_ROIcentroids = repmat(M_ROIcentroids, 1, size(ROIcentroids,1), 1);

% fomrula para la distancia 
M_distROI = sqrt(sum((M_ROIcentroids-permute(M_ROIcentroids, [2 1 3])).^2,3));


%% Procesar las Señales AD
% empieza a partir de los archivos matrices de todos los subjects

% long ROI Signals
Signals = load(strcat(path,'M_longRoiSignals_AD.mat')) ;
Signals = Signals.M_longRoiSignals_AD;
for i=1:size(Signals,3)
    sub = Signals(:,:,i) ;
    sub = transpose(sub) ;
    longMcorr(:,:,i) = corr(sub) ;
end

% short ROI Signals
Signals = load(strcat(path,'M_shortRoiSignals_AD.mat')) ;
Signals = Signals.M_shortRoiSignals_AD;
for i=1:size(Signals,3)
    sub = Signals(:,:,i) ;
    sub = transpose(sub) ;
    shortMcorr(:,:,i) = corr(sub) ;
end

Mcorr_AD = zeros(998,998,55) ;
Mcorr_AD(:,:,(1:19)) = shortMcorr ;
Mcorr_AD(:,:,(20:55)) = longMcorr ;


%% Procesar las Señales CN

clear longMcorr
clear shortMcorr
% long ROI Signals
Signals = load(strcat(path,'M_longRoiSignals_CN.mat')) ;
Signals = Signals.M_longRoiSignals_CN;
for i=1:size(Signals,3)
    sub = Signals(:,:,i) ;
    sub = transpose(sub) ;
    longMcorr(:,:,i) = corr(sub) ;
end

% short ROI Signals
Signals = load(strcat(path,'M_shortRoiSignals_CN.mat')) ;
Signals = Signals.M_shortRoiSignals_CN;
for i=1:size(Signals,3)
    sub = Signals(:,:,i) ;
    sub = transpose(sub) ;
    shortMcorr(:,:,i) = corr(sub) ;
end

Mcorr_CN= zeros(998,998,52) ;
Mcorr_CN(:,:,(1:10)) = shortMcorr ;
Mcorr_CN(:,:,(11:52)) = longMcorr ;
whos longMcorr
whos shortMcorr



%% ESTADISTICAS

% AD
% el valor medio de la correlacion entre ROI1 y ROI2 
% ( y asi para las 116x116 combinaciones entre ROI)
% promediando entre los 55 subjects
AD_M_media = mean(Mcorr_AD,3);
AD_M_std = std(Mcorr_AD,[],3);
AD_M_mediana = median(Mcorr_AD,3);
adIqrCC = quantile(Mcorr_AD,[0.25 0.75], 3);

% CN
% lo mismo que AD, para los subjects CN
CN_M_media = mean(Mcorr_CN,3);
CN_M_std = std(Mcorr_CN,[],3);
CN_M_mediana = median(Mcorr_CN,3);
cnIqrCC = quantile(Mcorr_CN,[0.25 0.75], 3);


%% CORELACION EN FUNCION DE LA DISTANCIA 

% M distancia entre ROI tiene datos repetidos 
% (porq es simetrica CRA la diagonal)
% aplico una mascara al triangulo de arriba
% asi me quedo solo con un datp de cada distancia 
% (y saco los duplicados, q estan en el triangulo simetrico)

MascaraTriangulo = logical(tril(ones(size(AD_M_media))));

% al aplicar mascara triangulo
% guarda los valores de la matriz uno atras de otro, en un vector
% donde se aplico la mascara, guarda un 1
% donde no estaba la mascara, guarad el valor original
AD_Vec_media = AD_M_media(MascaraTriangulo);
CN_Vec_media = CN_M_media(MascaraTriangulo);
Vec_distROI = M_distROI(MascaraTriangulo);

% Ordeno las distancias
% a cada distancia le corresponde una posicion en el vector
% hay q guardar esa posicion
% para despues encontrar a qué par de ROI correspondia esa ditancia
[Vec_distROI, posicion_en_vector] = sort(Vec_distROI);
AD_Vec_media = AD_Vec_media(posicion_en_vector);
CN_Vec_media = CN_Vec_media(posicion_en_vector);



%% FILTRO 
filterSize = 80;
coeffFilter = ones(1, filterSize)/filterSize;

AD_Vec_media_filt = filter(coeffFilter, 1, AD_Vec_media);
CN_Vec_media_filt = filter(coeffFilter, 1, CN_Vec_media);
eje_distancias = Vec_distROI ;


%%  PLOT

figure;
eje_distancias = Vec_distROI ;

% Sin filtro
subplot(2,1,1)
plot(eje_distancias, AD_Vec_media,'r');
hold on 
plot(eje_distancias, CN_Vec_media, 'b') ;
title ('Cross Correlation 998 ROI') ;
xlabel('Distance between ROI');
ylabel('Mean cross correlation');
legend('AD', 'CN')
hold off

% Filtrado
eje_distancias = Vec_distROI ;
subplot(2,1,2)
plot(eje_distancias, AD_Vec_media_filt,'r');
hold on 
plot(eje_distancias, CN_Vec_media_filt, 'b') ;
title ('Cross Correlation 998 ROI (filtered)') ;
xlabel('Distance between ROI');
ylabel('Mean cross correlation (filtered)');
legend('AD', 'CN')
hold off


%% encontrar L R 

L_indices = [] ;
R_indices = [] ;
% casos particualres 
Extras_indices = [] ;
 
for roi = 1: 998
    cual_region = ATLAS.roi_lbls(roi) ;    
    anat= ATLAS.anat_lbls(cual_region,:);
    % si esta etiquetado como Left
    if (contains(anat,'l'))
        L_indices = [L_indices , roi] ;
        texto = strcat(string(roi),'--->',anat,'--->','left H');
        disp(texto);
    % si esta etiquetado como Right
    elseif (contains(anat,'r'))
        R_indices = [R_indices , roi] ;
        texto = strcat(string(roi),'--->',anat,'--->','rigth H');
        disp(texto);
    % si no esta ni en L ni en R
    else
         disp('algo salio mal');
    end

end 


%% Matrices ipisilateral y contralateral

% consegui los indices de los ROI LH y ROI RH
% con esos indices busque la distancia entre ROI
% tambien busque la correlacion media

% ipsilitareal L: entre un ROI L y todos los otros ROI L
% (y asi para cada uno de los ROI L)

% ipsilateral R: lo mismo, entre ROI R

% contralteral: entre un ROI L y todos los otros ROI R
% (y asi para cada uno de los ROI L)


% inicializar matrices
M_ipsilateral_L = [] ;
M_ipsilateral_R = [] ;
M_contralateral = [] ;
AD_corr_ipsilateral_R = [] ;
AD_corr_ipsilateral_L = [] ;
AD_corr_contralateral = [] ;
CN_corr_ipsilateral_R = [] ;
CN_corr_ipsilateral_L = [] ;
CN_corr_contralateral = [] ;

% LH
for m = 1: length (L_indices)
    for n = 1: length (L_indices)
        roi1 = L_indices(m) ;
        roi2 = L_indices(n) ;
        M_ipsilateral_L(m,n)=M_distROI (roi1,roi2) ;
        AD_corr_ipsilateral_L(m,n) = AD_M_media(roi1,roi2) ;
        CN_corr_ipsilateral_L(m,n) = CN_M_media(roi1,roi2) ;
    end
end

% RH
for m = 1: length (R_indices)
    for n = 1: length (R_indices)
        roi1 = R_indices(m) ;
        roi2 = R_indices(n) ;
        M_ipsilateral_R(m,n)=M_distROI (roi1,roi2) ;
        AD_corr_ipsilateral_R(m,n) = AD_M_media(roi1,roi2) ;
        CN_corr_ipsilateral_R(m,n) = CN_M_media(roi1,roi2) ;
    end
end

% L-R
for m = 1: length (L_indices)
    for n = 1: length (R_indices)
        roi1 = L_indices(m) ;
        roi2 = R_indices(n) ;
        % filas ROI L
        % columnas ROI R
        M_contralateral(m,n)=M_distROI (roi1,roi2) ;
        AD_corr_contralateral(m,n)=AD_M_media(roi1,roi2) ;
        CN_corr_contralateral(m,n)=CN_M_media(roi1,roi2) ;
    end
end


disp ('arme todas Las amtrices') ;


%% Vectores ipisilateral

% mismo proceso q para el cerebro entero

% ipsilateral 
MascaraTriangulo = logical(tril(ones(size(AD_corr_ipsilateral_L))));
AD_Vec_ipsilateral_L = AD_corr_ipsilateral_L(MascaraTriangulo);
CN_Vec_ipsilateral_L = CN_corr_ipsilateral_L(MascaraTriangulo);
Vec_dist_ipsilateral_L = M_ipsilateral_L(MascaraTriangulo);
[Vec_dist_ipsilateral_L, posicion_en_vector] = sort(Vec_dist_ipsilateral_L);
AD_Vec_ipsilateral_L = AD_Vec_ipsilateral_L(posicion_en_vector);
CN_Vec_ipsilateral_L = CN_Vec_ipsilateral_L(posicion_en_vector);

% ipsilateral R
MascaraTriangulo = logical(tril(ones(size(AD_corr_ipsilateral_R))));
AD_Vec_ipsilateral_R = AD_corr_ipsilateral_R(MascaraTriangulo);
CN_Vec_ipsilateral_R = CN_corr_ipsilateral_R(MascaraTriangulo);
Vec_dist_ipsilateral_R = M_ipsilateral_R(MascaraTriangulo);
[Vec_dist_ipsilateral_R, posicion_en_vector] = sort(Vec_dist_ipsilateral_R);
AD_Vec_ipsilateral_R = AD_Vec_ipsilateral_R(posicion_en_vector);
CN_Vec_ipsilateral_R = CN_Vec_ipsilateral_R(posicion_en_vector);


%% Vectores contralteral 

% la matriz no es smetrica, las distancias no se repiten
% --> no aplico mascara 

AD_Vec_contralateral = reshape( AD_corr_contralateral , [],1);
CN_Vec_contralateral = reshape( CN_corr_contralateral , [],1);
Vec_dist_contralateral = reshape( M_contralateral , [],1);
[Vec_dist_contralateral, posicion_en_vector] = sort(Vec_dist_contralateral);

AD_Vec_contralateral = AD_Vec_contralateral(posicion_en_vector);
CN_Vec_contralateral = CN_Vec_contralateral(posicion_en_vector);



%% filtrado para un mejor plot

%filterSize = 80;
coeffFilter = ones(1, filterSize)/filterSize;

% LH
AD_ipsilateral_L_filt = filter(coeffFilter, 1, AD_Vec_ipsilateral_L);
CN_ipsilateral_L_filt = filter(coeffFilter, 1, CN_Vec_ipsilateral_L);

% RH
AD_ipsilateral_R_filt = filter(coeffFilter, 1, AD_Vec_ipsilateral_R);
CN_ipsilateral_R_filt = filter(coeffFilter, 1, CN_Vec_ipsilateral_R);

% LH-RH
AD_contralateral_filt = filter(coeffFilter, 1, AD_Vec_contralateral);
CN_contralateral_filt = filter(coeffFilter, 1, CN_Vec_contralateral);


%% PLOT IPSILATERAL

figure;


% LH
eje_distancias = Vec_dist_ipsilateral_L ;

% sin filtrar
subplot(2,3,1)
plot(eje_distancias, AD_Vec_ipsilateral_L,'r');
hold on 
plot(eje_distancias, CN_Vec_ipsilateral_L, 'b') ;
title ('Ipsilateral L-H') ;
xlabel('Distance between ROI');
ylabel('Mean cross correlation');
legend('AD', 'CN')
hold off

% Filtrado
subplot(2,3,4)
plot(eje_distancias, AD_ipsilateral_L_filt,'r');
hold on 
plot(eje_distancias, CN_ipsilateral_L_filt, 'b') ;
title ('Ipsilateral L-H (filtered)') ;
xlabel('Distance between ROI');
ylabel('Mean cross correlation (filtered)');
legend('AD', 'CN')


% RH
eje_distancias = Vec_dist_ipsilateral_R ;

% sin filtrar
subplot(2,3,2)
plot(eje_distancias, AD_Vec_ipsilateral_R,'r');
hold on 
plot(eje_distancias, CN_Vec_ipsilateral_R, 'b') ;
title ('Ipsilateral R-H') ;
xlabel('Distance between ROI');
ylabel('Mean cross correlation');
legend('AD', 'CN')
hold off

% Filtrado
subplot(2,3,5)
plot(eje_distancias, AD_ipsilateral_R_filt,'r');
hold on 
plot(eje_distancias, CN_ipsilateral_R_filt, 'b') ;
title ('Ipsilateral R-H (filtered)') ;
xlabel('Distance between ROI');
ylabel('Mean cross correlation (filtered)');
legend('AD', 'CN')


% LH - RH

eje_distancias = Vec_dist_contralateral ;

% sin filtrar
subplot(2,3,3)
plot(eje_distancias, AD_Vec_contralateral,'r');
hold on 
plot(eje_distancias, CN_Vec_contralateral, 'b') ;
title ('Contralateral') ;
xlabel('Distance between ROI');
ylabel('Mean cross correlation');
legend('AD', 'CN')
hold off

% Filtrado
subplot(2,3,6)
plot(eje_distancias, AD_contralateral_filt ,'r');
hold on 
plot(eje_distancias, CN_contralateral_filt, 'b') ;
title ('Contralateral (filtered)') ;
xlabel('Distance between ROI');
ylabel('Mean cross correlation (filtered)');
legend('AD', 'CN')




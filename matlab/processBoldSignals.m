clear all
%close all
dataPartitionPath = '/data/'; %'D:/'
adniPartitionPath = '/data_imaging/'; %'F:/'
%% ADD PATHS
dpabiPath = [dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/'];
addpath([dataPartitionPath 'UNSAM/Brain/dicm2nii/'])
addpath(genpath(dpabiPath))
addpath([dataPartitionPath 'UNSAM/Brain/spm12/spm12/'])
addpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/DPARSF/'])
%% TYPE OF PROCESSING
normaliseBoldSignals = 1;
firstCovThenStats = 1; % if 1, first the COV matrices are computed and then the mean matrices. 
%0, mean signals and the CC matrices.
%% AAL ATLAS ROIs
aalRois = load([dpabiPath '/Templates/aal_Labels.mat']);
% Compute the distance between ROIs:
roisIds = [aalRois.Reference{:,2}];
roisCentroids = reshape([aalRois.Reference{:,3}],3,[])';
% Should be in order, but we play safe,so index 1 is label 1, and remove 0.
roisIds(roisIds==0) = [];
[roisIds, indices] = sort(roisIds);
roisCentroids = roisCentroids(indices,:);
% Euclidean distance between ROIs:
roisCentroidsMatrix = reshape(roisCentroids, size(roisCentroids,1), 1, []);
roisCentroidsMatrix = repmat(roisCentroidsMatrix, 1, size(roisCentroids,1), 1);
roisDistanceMatrix = sqrt(sum((roisCentroidsMatrix-permute(roisCentroidsMatrix, [2 1 3])).^2,3));
%% PROCESS ROI SIGNALS
boldPath = '/home/martin/data/UNSAM/CEMSC3/ProcesamientoADNI/atlas AAL-20230309T142347Z-001/atlas AAL/';
%% AD
adPath = [boldPath 'ROISignalsAD(individuales)/'];
roiSignalsFilename = dir([adPath '*.mat']);
j = 1; % index to count total processed subjects.
adTotalSubjects = numel(roiSignalsFilename);
for i = 1 : adTotalSubjects
    roiSignalsSubject = load([adPath roiSignalsFilename(i).name]);
    roiSignalsSubject = roiSignalsSubject.ROISignals;
    if normaliseBoldSignals
        meanValueRoiSignals = mean(roiSignalsSubject, 1);
        stdValueRoiSignals = std(roiSignalsSubject, 1);
        roiSignalsSubject = (roiSignalsSubject-meanValueRoiSignals)./stdValueRoiSignals;
    end
    adLengthBoldSignals(i) = size(roiSignalsSubject,1);
    adBoldSignals(1:adLengthBoldSignals(i),:,i) = roiSignalsSubject;
    adCrossCorr(:,:,i) = corr(roiSignalsSubject);
    for j = 1 : numel(roisIds)
        ac = autocorr(roiSignalsSubject(:,j));
        adAutocorrLag1(j,i) = ac(2);
    end

end
adMeanCC = mean(adCrossCorr,3);
adStdCC = std(adCrossCorr,[],3);
adMedianCC = median(adCrossCorr,3);
adIqrCC = quantile(adCrossCorr,[0.25 0.75], 3);

%% CN
cnPath = [boldPath 'ROISignalsCN(individuales)/'];
cnRoiSignalsFilename = dir([cnPath '*.mat']);
j = 1; % index to count total processed subjects.
cnTotalSubjects = numel(cnRoiSignalsFilename);
for i = 1 : cnTotalSubjects
    roiSignalsSubject = load([cnPath cnRoiSignalsFilename(i).name]);
    roiSignalsSubject = roiSignalsSubject.ROISignals;
    if normaliseBoldSignals
        meanValueRoiSignals = mean(roiSignalsSubject, 1);
        stdValueRoiSignals = std(roiSignalsSubject, 1);
        roiSignalsSubject = (roiSignalsSubject-meanValueRoiSignals)./stdValueRoiSignals;
    end
    cnLengthBoldSignals(i) = size(roiSignalsSubject,1);
    cnBoldSignals(1:cnLengthBoldSignals(i),:,i) = roiSignalsSubject;
    cnCrossCorr(:,:,i) = corr(roiSignalsSubject);
    for j = 1 : numel(roisIds)
        ac = autocorr(roiSignalsSubject(:,j));
        cnAutocorrLag1(j,i) = ac(2);
    end

end
cnMeanCC = mean(cnCrossCorr,3);
cnStdCC = std(cnCrossCorr,[],3);
cnMedianCC = median(cnCrossCorr,3);
cnIqrCC = quantile(cnCrossCorr,[0.25 0.75], 3);
%%
figure; plot(cnAutocorrLag1')

%% CC AS FUNCTION OF THE DISTANCE
% GET THE CC AND DISTANCE IN A VECTOR.
% First the lower traingular part of the matrix as they are symmetric.
lowerDiagonalMask = logical(tril(ones(size(adMeanCC))));
adVecMeanCC = adMeanCC(lowerDiagonalMask);
cnVecMeanCC = cnMeanCC(lowerDiagonalMask);
% vecStdCC = stdCC(lowerDiagonalMask);
% vecMedianCC = medianCC(lowerDiagonalMask);
vecDistance = roisDistanceMatrix(lowerDiagonalMask);
% Sort the distances:
[vecDistance, indicesSortDistance] = sort(vecDistance);
adVecMeanCC = adVecMeanCC(indicesSortDistance);
cnVecMeanCC = cnVecMeanCC(indicesSortDistance);
% PLOT
figure;
plot(vecDistance, adVecMeanCC, vecDistance, cnVecMeanCC)
legend('AD', 'CN')
%% FILTER DATA
filterSize = 50;
coeffFilter = ones(1, filterSize)/filterSize;

adVecMeanCCfilt = filter(coeffFilter, 1, adVecMeanCC);
cnVecMeanCCfilt = filter(coeffFilter, 1, cnVecMeanCC);
figure;
plot(vecDistance, adVecMeanCCfilt, vecDistance, cnVecMeanCCfilt)
legend('AD', 'CN')
%% LEFT AND RIGHT
labelNames = {aalRois.Reference{:,1}};
indicesRoisR = find(endsWith(labelNames, '_R')>0);
indicesRoisL = find(endsWith(labelNames, '_L')>0);

for m = 1: length (indicesRoisL)
    for n = 1: length (indicesRoisL)
        roi1 = indicesRoisL(m) ;
        roi2 = indicesRoisL(n) ;
        M_ipsilateral_L(m,n)=roisDistanceMatrix (roi1,roi2) ;
        AD_corr_ipsilateral_L(m,n) = adMeanCC(roi1,roi2) ;
        CN_corr_ipsilateral_L(m,n) = cnMeanCC(roi1,roi2) ;
    end
end

for m = 1: length (indicesRoisR)
    for n = 1: length (indicesRoisR)
        roi1 = indicesRoisR(m) ;
        roi2 = indicesRoisR(n) ;
        M_ipsilateral_R(m,n)=roisDistanceMatrix (roi1,roi2) ;
        AD_corr_ipsilateral_R(m,n) = adMeanCC(roi1,roi2) ;
        CN_corr_ipsilateral_R(m,n) = cnMeanCC(roi1,roi2) ;
    end
end

% L-R
for m = 1: length (indicesRoisL)
    for n = 1: length (indicesRoisR)
        roi1 = indicesRoisL(m) ;
        roi2 = indicesRoisR(n) ;
        % filas ROI L
        % columnas ROI R
        M_contralateral(m,n)=roisDistanceMatrix (roi1,roi2) ;
        AD_corr_contralateral(m,n)=adMeanCC(roi1,roi2) ;
        CN_corr_contralateral(m,n)=cnMeanCC(roi1,roi2) ;
    end
end
%%
%% Vectores ipisilateral

% mismo proceso q para el cerebro entero
MascaraTriangulo = lowerDiagonalMask;
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
%% Ipsilateral Right
%indicesRoisRMatrix = reshape(indicesRoisR, size(indicesRoisR,2), 2, []);
%indicesRoisRMatrix = repmat(indicesRoisRMatrix, 1, size(indicesRoisR,1), 1);
%% COMPUTE AUTOCORRELATIONS
% Crop signals to the shortest time series to concatenate them:
adLengthShortestBold = min(adLengthBoldSignals);
cnLengthShortestBold = min(cnLengthBoldSignals);
adBoldSignals_cropped = permute(adBoldSignals, [2 1 3]);
cnBoldSignals_cropped = permute(cnBoldSignals, [2 1 3]);
adBoldSignals_cropped = adBoldSignals_cropped(:,1:adLengthShortestBold,:);
cnBoldSignals_cropped = cnBoldSignals_cropped(:,1:cnLengthShortestBold,:);
adConcatenatedSignals = reshape(adBoldSignals_cropped, size(adBoldSignals_cropped,1),[]);
cnConcatenatedSignals = reshape(cnBoldSignals_cropped, size(cnBoldSignals_cropped,1),[]);

% Autocorrelation for each ROI:
numLags = 100;
figure;
for i = 1 : size(adConcatenatedSignals,1)
    adAutoCorr = autocorr(adConcatenatedSignals(i,:), numLags);
    plot(adAutoCorr);
    hold on;
end
figure;
for i = 1 : size(cnConcatenatedSignals,1)
    cnAutoCorr = autocorr(cnConcatenatedSignals(i,:), numLags);
    plot(cnAutoCorr);
    hold on;
end
%% SAVE DATA
cnData.boldSignals = cnBoldSignals;
cnData.boldRoisCorr = cnCrossCorr;
cnData.boldLength = cnLengthBoldSignals;
cnData.boldNumSubjects = cnTotalSubjects;
save([boldPath 'cnData.mat'], 'cnData');

adData.boldSignals = adBoldSignals;
adData.boldRoisCorr = adCrossCorr;
adData.boldLength = adLengthBoldSignals;
adData.boldNumSubjects = adTotalSubjects;
save([boldPath 'adData.mat'], 'adData');
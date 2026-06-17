
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
addpath('./DPARSF/')
%%
adBoldSignals = load('/home/martin/data/UNSAM/Estudiantes/JulianFurios/DatosfMRI/ADNI_EA.mat');
cnBoldSignals = load('/home/martin/data/UNSAM/Estudiantes/JulianFurios/DatosfMRI/ADNI_CN.mat');
hcBoldSignals = load('/home/martin/data/UNSAM/Estudiantes/JulianFurios/DatosfMRI/HCP.mat');

adBoldSignals = adBoldSignals.ADNI_EA;
cnBoldSignals = cnBoldSignals.ADNI_CN;
hcpBoldSignals = hcBoldSignals.HCP;
%%
adBoldSignals = adBoldSignals(1:140,:,1:4);
cnBoldSignals = cnBoldSignals(1:140,:,:);
% %% CHECK THE IMAGES
% for i = 1 : 10%size(adBoldSignals,3)
%     fig = check_fMRI_bold_signals(adBoldSignals(:,:,i));
% end
% %% CHECK THE IMAGES
% for i = 1 : 10%size(cnBoldSignals,3)
%     fig = check_fMRI_bold_signals(cnBoldSignals(:,:,i));
% end
% %% CHECK THE IMAGES
% for i = 1 : 5%size(cnBoldSignals,3)
%     fig = check_fMRI_bold_signals(hcpBoldSignals(:,:,i));
% end
%%
for i = 1 : size(adBoldSignals,3)
    adCrossCorr(:,:,i) = corr(adBoldSignals(:,:,i));
end
%%
for i = 1 : size(cnBoldSignals,3)
    cnCrossCorr(:,:,i) = corr(cnBoldSignals(:,:,i));
end
%%
for i = 1 : size(hcpBoldSignals,3)
    hcpCrossCorr(:,:,i) = corr(hcpBoldSignals(:,:,i));
end

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
%%
adMeanCC = mean(adCrossCorr,3);
adStdCC = std(adCrossCorr,[],3);
adMedianCC = median(adCrossCorr,3);
adIqrCC = quantile(adCrossCorr,[0.25 0.75], 3);

cnMeanCC = mean(cnCrossCorr,3);
cnStdCC = std(cnCrossCorr,[],3);
cnMedianCC = median(cnCrossCorr,3);
cnIqrCC = quantile(cnCrossCorr,[0.25 0.75], 3);

hcpMeanCC = mean(hcpCrossCorr,3);
hcpStdCC = std(hcpCrossCorr,[],3);
hcpMedianCC = median(hcpCrossCorr,3);
hcpIqrCC = quantile(hcpCrossCorr,[0.25 0.75], 3);

%%
lowerDiagonalMask = logical(tril(ones(size(adMeanCC))));
adVecMeanCC = adMeanCC(lowerDiagonalMask);
cnVecMeanCC = cnMeanCC(lowerDiagonalMask);
hcpVecMeanCC = hcpMeanCC(lowerDiagonalMask);

adVecStdCC = adStdCC(lowerDiagonalMask);
cnVecStdCC = cnStdCC(lowerDiagonalMask);
hcpVecStdCC = hcpStdCC(lowerDiagonalMask);
% vecMedianCC = medianCC(lowerDiagonalMask);
vecDistance = roisDistanceMatrix(lowerDiagonalMask);
% Sort the distances:
[vecDistance, indicesSortDistance] = sort(vecDistance);
adVecMeanCC = adVecMeanCC(indicesSortDistance);
cnVecMeanCC = cnVecMeanCC(indicesSortDistance);
hcpVecMeanCC = hcpVecMeanCC(indicesSortDistance);
% PLOT
figure;
plot(vecDistance, adVecMeanCC, vecDistance, cnVecMeanCC, vecDistance, hcpVecMeanCC)
legend('AD', 'CN')
%%
%% FILTER DATA
filterSize = 20;
coeffFilter = ones(1, filterSize)/filterSize;

adVecMeanCCfilt = filter(coeffFilter, 1, adVecMeanCC);
cnVecMeanCCfilt = filter(coeffFilter, 1, cnVecMeanCC);
hcpVecMeanCCfilt = filter(coeffFilter, 1, hcpVecMeanCC);
figure;
plot(vecDistance, adVecMeanCCfilt, vecDistance, cnVecMeanCCfilt,...
    vecDistance, hcpVecMeanCCfilt)
legend('AD', 'CN')

%%
figure;
errorbar(vecDistance,adVecMeanCCfilt,adVecStdCC);
hold on;
errorbar(vecDistance,cnVecMeanCCfilt,cnVecStdCC);
errorbar(vecDistance,hcpVecMeanCCfilt,hcpVecStdCC);
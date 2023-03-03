clear all
%close all

%% ADD PATHS
dpabiPath = 'D:\UNSAM\Brain\DPABI_V6.2_220915\';
addpath('D:\UNSAM\Brain\dicm2nii\')
addpath(genpath(dpabiPath))
addpath('D:\UNSAM\Brain\spm12\spm12\')
addpath('D:\UNSAM\Brain\DPABI_V6.2_220915\DPARSF')
%% TYPE OF PROCESSING
firstCovThenStats = 1; % if 1, first the COV matrices are computed and then the mean matrices. 
%0, mean signals and the CC matrices.
%% AAL ATLAS ROIs
aalRois = load([dpabiPath '\Templates\aal_Labels.mat']);
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
boldPath = 'D:\UNSAM\CEMSC3\ProcesamientoADNI\Señales BOLD\';
%% AD
adBatches = 7;
j = 1; % index to count total processed subjects.
adTotalSubjects = 0;
for i = 1 : adBatches
    adBatch = load([boldPath sprintf('AD_tanda%d.mat', i)]);
    cantSubjects = adBatch.cant_subjects;
    adTotalSubjects = adTotalSubjects + cantSubjects;
    % Process each subject:
    for s = 1 : cantSubjects
        adLengthBoldSignals(j) = size(adBatch.S(:,:,s),1);
        adBoldSignals(1:adLengthBoldSignals(j),:,j) = adBatch.S(:,:,s);
        adCrossCorr(:,:,j) = corr(adBatch.S(:,:,s));
        j = j + 1;
    end
    %adBoldSignals(j:j+cantSubjects-1) = adBatch.S;
    %j = j + cantSubjects;
end
% if firstCovThenStats
%     % Process each subject:
%     for s = 1 : cantSubjects
%         crossCorr(:,:,j) = corr(adBatch.S(:,:,s));
%         j = j + 1;
%     end
% else
%     % Process each subject:
%     meanBoldRois = mean(adBatch.S, 3);
%     meanCC = corr(meanBoldRois);
%     %stdCC = std(crossCorr,[],3);
%     %medianCC = median(crossCorr,3);
%     %iqrCC = quantile(crossCorr,[0.25 0.75], 3);
% end
adMeanCC = mean(adCrossCorr,3);
adStdCC = std(adCrossCorr,[],3);
adMedianCC = median(adCrossCorr,3);
adIqrCC = quantile(adCrossCorr,[0.25 0.75], 3);
%% CN
cnBatches = 6;
j = 1; % index to count total processed subjects.
cnTotalSubjects = 0;
for i = 1 : cnBatches
    cnBatch = load([boldPath sprintf('CN_tanda%d.mat', i)]);
    cantSubjects = cnBatch.cant_subjects;
    cnTotalSubjects = cnTotalSubjects + cantSubjects;
    % Process each subject:
    for s = 1 : cantSubjects
        cnLengthBoldSignals(j) = size(cnBatch.S(:,:,s),1);
        cnBoldSignals(1:cnLengthBoldSignals(j),:,j) = cnBatch.S(:,:,s);
        cnCrossCorr(:,:,j) = corr(cnBatch.S(:,:,s));
        j = j + 1;
    end
    %adBoldSignals(j:j+cantSubjects-1) = adBatch.S;
    %j = j + cantSubjects;
end
% if firstCovThenStats
%     % Process each subject:
%     for s = 1 : cantSubjects
%         crossCorr(:,:,j) = corr(adBatch.S(:,:,s));
%         j = j + 1;
%     end
% else
%     % Process each subject:
%     meanBoldRois = mean(adBatch.S, 3);
%     meanCC = corr(meanBoldRois);
%     %stdCC = std(crossCorr,[],3);
%     %medianCC = median(crossCorr,3);
%     %iqrCC = quantile(crossCorr,[0.25 0.75], 3);
% end
cnMeanCC = mean(cnCrossCorr,3);
cnStdCC = std(cnCrossCorr,[],3);
cnMedianCC = median(cnCrossCorr,3);
cnIqrCC = quantile(cnCrossCorr,[0.25 0.75], 3);
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
filterSize = 20;
coeffFilter = ones(1, filterSize)/filterSize;

adVecMeanCCfilt = filter(coeffFilter, 1, adVecMeanCC);
cnVecMeanCCfilt = filter(coeffFilter, 1, cnVecMeanCC);
figure;
plot(vecDistance, adVecMeanCCfilt, vecDistance, cnVecMeanCCfilt)
legend('AD', 'CN')
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
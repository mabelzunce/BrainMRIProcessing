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
%% DATA PATHS
adniPath = [adniPartitionPath '/ADNIdata/'];
processedDataPath = [adniPath '/Processed2/'];
filename998Rois = '/data/UNSAM/Brain/DSI_enhanced.mat';
outputPath = [adniPath '/fMRI_Analysis/'];
%% FOLDERS FOR DPARSF
t1NameNifti = 'T1Img';
fmriNameNifti = 'FunImg';
fmriSliceCorrectedNameNifti = 'FunImgA';
fmriFullyPreprocessedNameNifti = 'FunImgARWSD'; %Filtered_4DVolume.nii
fmriFullyPreprocessedBoldSignals = '/Results/ROISignals_FunImgARWSDCF/';
fieldmapBaseNameNifti = 'FieldMap';
fieldmapPhaseNameNifti = 'PhaseDiffImg';
fieldmapMagNameNifti = 'Maginute1Img';
%% LOAD DATA FROM prepareADNIDataFromDicom
load([processedDataPath 'workspaceProcessingScript']);
%% TYPE OF PROCESSING
firstCovThenStats = 1; % if 1, first the COV matrices are computed and then the mean matrices. 
%0, mean signals and the CC matrices.
%% AAL ATLAS ROIs
aalImage = niftiread([dpabiPath '/Templates/aal.nii']);
aalImage_info = niftiinfo([dpabiPath '/Templates/aal.nii']);
aalImage_3mm = imresize3(aalImage, ceil(size(aalImage).*aalImage_info.PixelDimensions./[3 3 3]), 'nearest');
aalRois = load([dpabiPath '/Templates/aal_Labels.mat']);
% Resample image
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
%% 998 ATLAS ROIs
atlas998Rois = load(filename998Rois);
% las coords x y z de cada ROI, puestas una atras de la otra
coords998Rois = atlas998Rois.roi_xyz_avg ;

% las coords x y z de cada ROI, puestas en columnas
rois998Centroids = transpose(coords998Rois) ;
% Compute the distance between ROIs:
rois998Ids = atlas998Rois.roi_lbls;
% Sort the ROI indices:
[rois998Ids, indices] = sort(rois998Ids);
rois998Centroids = rois998Centroids(indices,:);
% Euclidean distance between ROIs:
rois998CentroidsMatrix = reshape(rois998Centroids, size(rois998Centroids,1), 1, []);
rois998CentroidsMatrix = repmat(rois998CentroidsMatrix, 1, size(rois998Centroids,1), 1);
rois998DistanceMatrix = sqrt(sum((rois998CentroidsMatrix-permute(rois998CentroidsMatrix, [2 1 3])).^2,3));
%% READ ALL THE BOLD SIGNALS FOR EACH GROUP AND SCANNER
% group them by clinical group:
indexSubject = 1;
for s = 1 : numel(scanners)
    for g = 1 : numel(groups)
        dparsfThisGroupBasePath = [dparsfPerGroupPath groups{g} '_' scanners{s} '/'];
        dparsfThisGroupFunPath = [dparsfThisGroupBasePath fmriNameNifti '/'];
        dparsfThisGroupFunSliceCorrectedPath = [dparsfThisGroupBasePath fmriSliceCorrectedNameNifti '/'];
        dparsfThisGroupFunFullyPrecorrectedPath = [dparsfThisGroupBasePath fmriFullyPreprocessedNameNifti '/'];        
        listDir = dir(dparsfThisGroupFunPath);
        % Get all the subject names:
        subjectNamesThisGroup = listDir(3:end); % Remove ., ..
        for i = 1 : numel(subjectNamesThisGroup)
            subjectNames{indexSubject} = subjectNamesThisGroup(i).name;
            outputPathThisSubject = [outputPath subjectNames{indexSubject} '/'];
            if ~isdir(outputPathThisSubject)
                mkdir(outputPathThisSubject);
            end
            % Read images:
            filenameFunImg = [dparsfThisGroupFunPath subjectNames{indexSubject} '/'...
                dir([dparsfThisGroupFunPath subjectNames{indexSubject} '/*.nii']).name];
            flienameSliceCorrectedImg = [dparsfThisGroupFunSliceCorrectedPath subjectNames{indexSubject} '/'...
                dir([dparsfThisGroupFunSliceCorrectedPath subjectNames{indexSubject} '/*.nii']).name];
            filenameFullyPreprocessedImg = [dparsfThisGroupFunFullyPrecorrectedPath subjectNames{indexSubject} '/'...
                dir([dparsfThisGroupFunFullyPrecorrectedPath subjectNames{indexSubject} '/*.nii']).name];
            fmriRawImage = niftiread(filenameFunImg);
            fmriSliceTimingCorrectedImage = niftiread(flienameSliceCorrectedImg);
            fmriFullyPreprocessedImage = niftiread(filenameFullyPreprocessedImg);
            headerinfo = niftiinfo(filenameFullyPreprocessedImg);
            % Get ROIs:
            boldLength(indexSubject) = size(fmriFullyPreprocessedImage,4);
            for j = 1 : max(aalImage_3mm(:))
                maskRoi = aalImage_3mm == j;
                boldSignalsPerRoi(1:boldLength(indexSubject), j, indexSubject) = squeeze(sum(fmriFullyPreprocessedImage.*maskRoi, [1 2 3]))./...
                    repmat(sum(maskRoi(:)), [size(fmriFullyPreprocessedImage,4) 1]);
                 % Compute autocorr:
                    
                 cnAutocorrLag1(j,i) = ac(2);
            end
            % Save image to verify:
            SaveOverlayImageWithRois(fmriFullyPreprocessedImage, aalImage_3mm, [outputPathThisSubject 'meanfMriAAL.gif']);
            % Compute cross correlation:
            boldCrossCorr(:,:,indexSubject) = corr(boldSignalsPerRoi(:,:,indexSubject));
            % Read the output from DPARSF:
            boldSignalsPerRoiDparsfThisSubject = load([dparsfThisGroupBasePath fmriFullyPreprocessedBoldSignals 'ROISignals_' subjectNames{indexSubject} '.mat']);
            boldCorrDparsfThisSubject = load([dparsfThisGroupBasePath fmriFullyPreprocessedBoldSignals 'ROICorrelation_' subjectNames{indexSubject} '.mat']);
            boldCorrZDparsfThisSubject = load([dparsfThisGroupBasePath fmriFullyPreprocessedBoldSignals 'ROICorrelation_FisherZ_' subjectNames{indexSubject} '.mat']);
            roisCentreDparsfThisSubject = load([dparsfThisGroupBasePath fmriFullyPreprocessedBoldSignals 'ROI_CenterOfMass_' subjectNames{indexSubject} '.mat']);
            boldSignalsDparsf(1:boldLength(indexSubject),:,indexSubject) = boldSignalsPerRoiDparsfThisSubject.ROISignals;
            boldCorrMatrixDparsf(:,:,indexSubject) = boldCorrDparsfThisSubject.ROICorrelation;
            boldCorrZMatrixDparsf(:,:,indexSubject) = boldCorrZDparsfThisSubject.ROICorrelation_FisherZ;
            indexSubject = indexSubject + 1;
        end
    end
end
fullProcessedResults.boldSignals = boldSignalsPerRoi;
fullProcessedResults.boldLength = boldLength;
fullProcessedResults.boldSignalsDparsf = boldSignalsDparsf;
fullProcessedResults.boldCorrMatrix = boldCrossCorr;
fullProcessedResults.boldCorrMatrixDparsf = boldCorrMatrixDparsf;
save([boldPath 'fullProcessedResults.mat'], 'fullProcessedResults');
%% CHECK CORRELATION MATRIX
figure;
for i = 1 : size(boldCrossCorr,3)
    subplot(1,2,1);
    imagesc(boldCrossCorr(:,:,i));
    subplot(1,2,2);
    imagesc(boldCorrMatrixDparsf(:,:,i));
    pause(0.5);
end
%% CHECK CORRELATION MATRIX
% First the lower traingular part of the matrix as they are symmetric.
lowerDiagonalMask = logical(tril(ones(size(boldCrossCorr(:,:,1)))));
figure;
for i = 1 : size(boldCrossCorr,3)
    matrixCC = boldCrossCorr(:,:,i);
    vecCC = matrixCC(lowerDiagonalMask);
    subplot(1,2,1);
    hist(vecCC);
    matrixCC = boldCorrMatrixDparsf(:,:,i);
    vecCC = matrixCC(lowerDiagonalMask);
    subplot(1,2,2);
    hist(vecCC);
    pause(0.5);
end
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
clear all
close all

dataPartitionPath = '/data/'; %'D:/'
adniPartitionPath = '/data_imaging/'; %'F:/'
currentPath = pwd;
% Load data:
loadData = 1;

%% ADD PATHS
addpath([dataPartitionPath 'UNSAM/Brain/dicm2nii/'])
addpath(genpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/']))
addpath([dataPartitionPath 'UNSAM/Brain/spm12/spm12/'])
addpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/DPARSF/'])
%% PREPROCESSED DATA PATH
dataPath = '/home/martin/data_imaging/ADNIdata/fMRI_first_batch_to_process_2023_04_27/Procesadas_2mm_highpass/';
dataPath = '/home/martin/data_imaging/CovidProject/Estudio2/PreprocessedMRI/DPARSF/';
preprocessedImageFilenames = 'Filtered_4DVolume.nii';
preprocessedFolder = 'FunImgARWSDCF';
outputNormalizedImageSubdir = '/FunImgARWSDCFN/';
preprocessedDataPath = fullfile(dataPath, preprocessedFolder);
%% OUTPUTPATH
suffixROIFilenames = 'ROISignals_';
outputPathImages = fullfile(dataPath, outputNormalizedImageSubdir);
outputPathSignals = fullfile(dataPath, 'Results', ['ROISignals_' outputNormalizedImageSubdir]);

if ~isfolder(outputPathImages)
    mkdir(outputPathImages);
end

if ~isfolder(outputPathSignals)
    mkdir(outputPathSignals);
end
%% MASKS
masksDparsfSubfolder = '/Masks/';
masksDparsfPath = [dataPath '/' masksDparsfSubfolder '/'];
brainMaskFilename = 'AllResampled_BrainMask_05_91x109x91.nii';
gmMaskFilename = 'AllResampled_GreyMask_02_91x109x91.nii';
wmMaskFilename = 'AllResampled_WhiteMask_09_91x109x91.nii';
csfMaskFilename = 'AllResampled_CsfMask_07_91x109x91.nii';
%% GRAND MEAN NORMALIZATION, CHECK COVERAGE AND SPATIAL ALGINMENT
saveNormalizedImages = 0;
createNormalizedRoiSignals = 0;
scaleMeanValue = 10000;
maxTimePoints = 200; % FOR THE STANDARD PROTOCOL.
% read masks:
brain_mask = logical(niftiread(fullfile(masksDparsfPath,brainMaskFilename)));
brain_mask_4d = repmat(brain_mask,1,1,1,maxTimePoints);
gm_mask = logical(niftiread(fullfile(masksDparsfPath,gmMaskFilename)));
gm_mask_4d = repmat(gm_mask,1,1,1,maxTimePoints);
wm_mask = logical(niftiread(fullfile(masksDparsfPath,wmMaskFilename)));
wm_mask_4d = repmat(wm_mask,1,1,1,maxTimePoints);
csf_mask = logical(niftiread(fullfile(masksDparsfPath,csfMaskFilename)));
csf_mask_4d = repmat(csf_mask,1,1,1,maxTimePoints);
% Non brain mask
non_brain_mask = ~brain_mask;
se = strel('cube',10);
non_brain_mask = imerode(non_brain_mask,se);
non_brain_mask_4d = repmat(non_brain_mask,1,1,1,maxTimePoints);

% Get a list of preprocessed fMRI files in the directory
fmriSubjects = dir(preprocessedDataPath);
fmriSubjects = fmriSubjects([fmriSubjects(:).isdir]);
% Remove the first two (.,..)
fmriSubjects = {fmriSubjects(3:end).name};
% MAtrix to store all the data:
SignalsAllSubjects = {};
% Iterate over each preprocessed fMRI file
for i = 1 : numel(fmriSubjects)
    % Read the fMRI data
    niftiFilename = dir(fullfile(preprocessedDataPath, fmriSubjects{i}, '*.nii'));
    outputFilename = fullfile(dataPath, outputNormalizedImageSubdir, fmriSubjects{i}, preprocessedImageFilenames);
    fmriVolume = niftiread(fullfile(preprocessedDataPath, fmriSubjects{i}, niftiFilename.name));
    fmriInfo = niftiinfo(fullfile(preprocessedDataPath, fmriSubjects{i}, niftiFilename.name));
    
    % Get the voxel dimensions
    voxelSize_mm3 = fmriInfo.PixelDimensions(1:3);
    tR_sec = fmriInfo.PixelDimensions(4);
    imageSize_voxels = fmriInfo.ImageSize(1:3);
    timePoints = fmriInfo.ImageSize(4);   
    grandMeanValue(i) = mean(fmriVolume(brain_mask_4d(:,:,:,size(fmriVolume,4))));
    gmMeanValue(i) = mean(fmriVolume(gm_mask_4d(:,:,:,size(fmriVolume,4))));
    stdMeanValue(i) = std(fmriVolume(gm_mask_4d(:,:,:,size(fmriVolume,4))));
    wmMeanValue(i) = mean(fmriVolume(wm_mask_4d(:,:,:,size(fmriVolume,4))));
    csfMeanValue(i) = mean(fmriVolume(csf_mask_4d(:,:,:,size(fmriVolume,4))));
    fmriVolume_norm = fmriVolume ./ grandMeanValue(i) .* scaleMeanValue; % Grand mean intensity normalization.
    if saveNormalizedImages && ~exist(outputFilename) % if exists dont'w write. %TODO
        if ~exist(fullfile(dataPath, outputNormalizedImageSubdir, fmriSubjects{i}))
            mkdir(fullfile(dataPath, outputNormalizedImageSubdir, fmriSubjects{i}))
        end
        niftiwrite(fmriVolume_norm, outputFilename, fmriInfo,  'Compressed', 1);
    end
    grandMeanValuePost(i) = mean(fmriVolume_norm(brain_mask_4d(:,:,:,size(fmriVolume,4))));
    gmMeanValuePost(i) = mean(fmriVolume_norm(gm_mask_4d(:,:,:,size(fmriVolume,4))));
    stdMeanValuePost(i) = std(fmriVolume_norm(gm_mask_4d(:,:,:,size(fmriVolume,4))));
    wmMeanValuePost(i) = mean(fmriVolume_norm(wm_mask_4d(:,:,:,size(fmriVolume,4))));
    csfMeanValuePost(i) = mean(fmriVolume_norm(csf_mask_4d(:,:,:,size(fmriVolume,4))));
    % Compute iamge quality metrics: SNR, tSNR and CNR.
    meanBrainTimeCourses = mean(fmriVolume_norm,4);
    stdBrainTimeCourses = std(fmriVolume_norm,0,4);
    tSNR(i) = mean(meanBrainTimeCourses(brain_mask)./stdBrainTimeCourses(brain_mask));
    tSNR_gm(i) = mean(meanBrainTimeCourses(gm_mask)./stdBrainTimeCourses(gm_mask));
    SNR(i) = grandMeanValuePost(i)./std(fmriVolume_norm(non_brain_mask_4d(:,:,:,size(fmriVolume,4))));
    SNR_gm(i) = gmMeanValuePost(i)./std(fmriVolume_norm(non_brain_mask_4d(:,:,:,size(fmriVolume,4))));
    CNR(i) = (gmMeanValuePost(i)-wmMeanValuePost(i))./std(fmriVolume_norm(non_brain_mask_4d(:,:,:,size(fmriVolume,4))));
    %tSNR = mean(fmriVolume_norm(gm_mask_4d(:,:,:,size(fmriVolume,4))));
end

%% PLOT MEAN VALUES

% PLOTS
figure;
plot(grandMeanValue);
hold on;
plot(gmMeanValue);
plot(wmMeanValue);
plot(csfMeanValue);
figure; plot(stdMeanValue)


% PLOTS
figure;
plot(grandMeanValuePost);
hold on;
plot(gmMeanValuePost);
plot(wmMeanValuePost);
plot(csfMeanValuePost);
figure; plot(stdMeanValuePost)
%% PLOT SNR METRICS
outliersSNR = find(isoutlier(SNR)>0);
outliersSNR_gm = find(isoutlier(SNR_gm)>0);
outlierstSNR_gm = find(isoutlier(tSNR_gm)>0);
figure;
subplot(1,3,1);
plot(SNR);
hold on
plot(outliersSNR, SNR(outliersSNR),'x');
subplot(1,3,2);
plot(SNR_gm);
hold on
plot(outliersSNR_gm, SNR_gm(outliersSNR_gm),'x');
subplot(1,3,3);
plot(tSNR_gm);
hold on
plot(outlierstSNR_gm, tSNR_gm(outlierstSNR_gm),'x');
%% GET NAMES TO CHECK OUTLIERS
% Get outliers with lowe SNR:
outliersSNR = find(isoutlier(SNR) & (SNR < median(SNR, "omitnan")));
outliersSNR_gm = find(isoutlier(SNR_gm) & (SNR_gm < median(SNR_gm, "omitnan")));
outlierstSNR_gm = find(isoutlier(tSNR_gm) & (tSNR_gm < median(tSNR_gm, "omitnan")));
% Print subjects names:
fprintf('Subjects with low SNR: %s\n', fmriSubjects{outliersSNR_gm})
fprintf('Subjects with low tSNR: %s\n', fmriSubjects{outlierstSNR_gm})
% Prints subjects with NaN:
fprintf('Subjects with NaN SNR: %s\n', fmriSubjects{isnan(SNR_gm)})
fprintf('Subjects with NaN tSNR: %s\n', fmriSubjects{isnan(tSNR_gm)})
%% CREATE A TABLE WITH THESE METRICS
tableImageQuality = table();
tableImageQuality.subjectsNames = fmriSubjects';
tableImageQuality.grandMeanValue = grandMeanValue';
tableImageQuality.grandMeanValueNorm = grandMeanValuePost';
tableImageQuality.SNR = SNR';
tableImageQuality.SNR_gm = SNR_gm';
tableImageQuality.tSNR_gm = tSNR_gm';
tableImageQuality.CNR = CNR';
%%
save([dataPath 'tableImageQuality.mat'], "tableImageQuality")
clear all
close all

dataPartitionPath = '/data/'; %'D:/'
imagingPartitionPath = '/data_imaging/'; %'F:/'
currentPath = pwd;
% Load data:
loadData = 1;

%% ADD PATHS
addpath([dataPartitionPath 'UNSAM/Brain/dicm2nii/'])
addpath(genpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/']))
addpath([dataPartitionPath 'UNSAM/Brain/spm12/spm12/'])
addpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/DPARSF/'])
addpath('../DPARSF/')
%% DATA PATHs
dataPath = [imagingPartitionPath '/CovidProject/Estudio/PreprocessedMRI/'];
%dataPath = [imagingPartitionPath '/CovidProject/Estudio2/PreprocessedMRI/'];
niftiPath = [dataPath '/Nifti/'];
dparsfDataPath = [imagingPartitionPath '/CovidProject/Estudio/PreprocessedMRI/DPARSF/'];
%dparsfDataPath = [imagingPartitionPath '/CovidProject/Estudio2/PreprocessedMRI/DPARSF/'];
outputCheckPath = [dparsfDataPath '/check/'];
if ~isdir(outputCheckPath)
    mkdir(outputCheckPath)
end
indexScanner = 1; % Siemens=1, GE=2, Philips=3.
% Filename with the available MRI data:
filenameMriInfo = [dataPath 'mriInfoAndProcessing_2024_06_03.mat'];
%filenameMriInfo = [dataPath 'mriInfoAndProcessing_2024_09_16.mat'];
%% DATA INFO
mriInfo = load(filenameMriInfo);
%% CONFIG
bandPassFilter = 0; %If not band pass, only high pass

%% FOLDERS FOR DPARSF IMAGES
t1NameNifti = 'T1Img';
fmriNameNifti = 'FunImg';
fmriSliceCorrectedNameNifti = 'FunImgA';
fmriFullyPreprocessedNameNifti = 'FunImgARWSD';
fieldmapBaseNameNifti = 'FieldMap';
fieldmapPhaseNameNifti = 'PhaseDiffImg';
fieldmapMagNameNifti = 'Maginute1Img';

t1DparsfPath = [dparsfDataPath '/' t1NameNifti '/'];
fmriDparsfPath = [dparsfDataPath '/' fmriNameNifti '/'];

roiSignalsDparsfSubdir = '/Results/ROISignals_FunImgARWSDCF/';
roiSignalsDparsfPath = [dparsfDataPath '/' roiSignalsDparsfSubdir '/'];
suffixROIFilenames = 'ROISignals_';

masksDparsfSubfolder = '/Masks/';
masksDparsfPath = [dparsfDataPath '/' masksDparsfSubfolder '/'];
brainMaskFilename = 'AllResampled_BrainMask_05_91x109x91.nii';
gmMaskFilename = 'AllResampled_GreyMask_02_91x109x91.nii';
wmMaskFilename = 'AllResampled_WhiteMask_09_91x109x91.nii';
csfMaskFilename = 'AllResampled_CsfMask_07_91x109x91.nii';
overwriteNifti = 1; % If 0 won't process existing Nfiti data in the DPARSF folder
format = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';

%% CHECK DATA
listDir = dir(fmriDparsfPath);
% Get all subjects that are similar:
i = 1;
for k = 1 : numel(listDir)-2
    tempName = listDir(k+2).name;
    % Check in MriInfo the data:
    indexThisSubject = find(strcmp(tempName, mriInfo.casesToProcess));
    if ~isempty(indexThisSubject)
        subjectNames{i} = tempName;
        imageSize_voxels(i,:) = mriInfo.fMRI_imageSize_voxels(indexThisSubject,:);
        timePoints(i) = imageSize_voxels(i,4);
        numSlices(i) = imageSize_voxels(i,3);
        i = i + 1;
    end
end

% %% CHECK DATA
% for i = 1 : numel(subjectNames)
%     if exist([roiSignalsDparsfPath suffixROIFilenames subjectNames{i} '.mat'])
%         signals = load([roiSignalsDparsfPath suffixROIFilenames subjectNames{i} '.mat']);
%         fig = check_fMRI_bold_signals(signals.ROISignals);
%         saveas(gca, fullfile(outputCheckPath, subjectNames{i}), 'tif');
%         close all
%     end
% end
%% GRAND MEAN NORMALIZATION, CHECK COVERAGE AND SPATIAL ALGINMENT
saveNormalizedImages = 1;
scaleMeanValue = 10000;
preprocessedImageSubdir = '/FunImgARWSDCF/';
outputNormalizedImageSubdir = '/FunImgARWSDCFN/';
preprocessedImageFilenames = 'Filtered_4DVolume.nii';
preprocessedImagesDparsfPath = [dparsfDataPath '/' preprocessedImageSubdir '/'];
% read masks:
maxTimePoints = max(mriInfo.fMRI_imageSize_voxels(:,4));
brain_mask = logical(niftiread(fullfile(masksDparsfPath,brainMaskFilename)));
brain_mask_4d = repmat(brain_mask,1,1,1,maxTimePoints);
gm_mask = logical(niftiread(fullfile(masksDparsfPath,gmMaskFilename)));
gm_mask_4d = repmat(gm_mask,1,1,1,maxTimePoints);
wm_mask = logical(niftiread(fullfile(masksDparsfPath,wmMaskFilename)));
wm_mask_4d = repmat(wm_mask,1,1,1,maxTimePoints);
csf_mask = logical(niftiread(fullfile(masksDparsfPath,csfMaskFilename)));
csf_mask_4d = repmat(csf_mask,1,1,1,maxTimePoints);
j=1;
for i = 1 : numel(subjectNames)
    if exist(fullfile(preprocessedImagesDparsfPath, subjectNames{i}, preprocessedImageFilenames))
        image = niftiread(fullfile(preprocessedImagesDparsfPath, subjectNames{i}, preprocessedImageFilenames));
        info = niftiinfo(fullfile(preprocessedImagesDparsfPath, subjectNames{i}, preprocessedImageFilenames));       
        image_z = (image - mean(image, 4)) ./ std(image, [], 4);
        grandMeanValue(i) = mean(image(brain_mask_4d(:,:,:,size(image,4))));
        gmMeanValue(i) = mean(image(gm_mask_4d(:,:,:,size(image,4))));
        stdMeanValue(i) = std(image(gm_mask_4d(:,:,:,size(image,4))));
        wmMeanValue(i) = mean(image(wm_mask_4d(:,:,:,size(image,4))));
        csfMeanValue(i) = mean(image(csf_mask_4d(:,:,:,size(image,4))));
        image_norm = image ./ grandMeanValue(i) .* scaleMeanValue; % Grand mean intensity normalization.
        if saveNormalizedImages
            if ~exist(fullfile(dparsfDataPath, outputNormalizedImageSubdir, subjectNames{i}))
                mkdir(fullfile(dparsfDataPath, outputNormalizedImageSubdir, subjectNames{i}))
            end
            niftiwrite(image_norm,fullfile(dparsfDataPath, outputNormalizedImageSubdir, subjectNames{i}, [preprocessedImageFilenames '.gz']), info)
        end
        gmVoxels(j,:) = image(gm_mask_4d(:,:,:,size(image,4)));
        preprocessedImages(:,:,:,j) = mean(image_norm, 4);
        %meanBrainSignal = 
        j = j + 1;
    end
end
% Get mean and std images:
info = niftiinfo(fullfile(preprocessedImagesDparsfPath, subjectNames{i}, preprocessedImageFilenames));
info.ImageSize = [91 109 91];
info.PixelDimensions = [2 2 2];
meanImage = mean(preprocessedImages,4);
stdImage = std(preprocessedImages,[],4);
niftiwrite(meanImage, fullfile(outputCheckPath,'mean_image.nii.gz'), info);
niftiwrite(stdImage, fullfile(outputCheckPath, 'std_image.nii.gz'), info);

% PLOTS
figure;
plot(grandMeanValue);
hold on;
plot(gmMeanValue);
plot(wmMeanValue);
plot(csfMeanValue);
figure; plot(stdMeanValue)
figure; imshow(gmVoxels,[])
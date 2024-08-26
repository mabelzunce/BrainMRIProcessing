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
addpath(genpath('../'))
%% PREPROCESSED DATA PATH
dataPath = '/home/martin/data_imaging/ADNIdata/fMRI_first_batch_to_process_2023_04_27/Procesadas_2mm_highpass/';
dataPath = '/home/martin/data_imaging/CovidProject/Estudio/PreprocessedMRI/DPARSF/';
preprocessedFolder = 'FunImgARWSDCF';
preprocessedDataPath = fullfile(dataPath, preprocessedFolder);
indexScanner = 1; % Siemens=1, GE=2, Philips=3.
%% OUTPUTPATH
outputPath = [dataPath 'SignalsSchaefer2018_1000Parcels_17Networks/'];
outputPathSignals = [outputPath '/ROISignals/'];
outputPathImages = [outputPath '/CheckImages/'];
suffixROIFilenames = 'ROISignals_';
if ~isfolder(outputPathSignals)
    mkdir(outputPathSignals);
end
if ~isfolder(outputPathImages)
    mkdir(outputPathImages);
end
%% Schaefer Parcelation
atlasPAth = '/home/martin/data/UNSAM/Brain/Atlases/';
filenameAtlasNifti =  fullfile(atlasPAth, 'Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.nii.gz');
filenameAtlasRois =  fullfile(atlasPAth, 'Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
filenameAtlasColormap = fullfile(atlasPAth, 'Schaefer2018_1000Parcels_Kong2022_17Networks_order.txt');
% read atlas volume:
atlasVolume = niftiread(filenameAtlasNifti);
atlasRois = readtable(filenameAtlasRois);
% las coords x y z de cada ROI en MNI space:
coordsRois = [atlasRois.R atlasRois.A atlasRois.S];
% Compute the ROIs ids:
roisIds = atlasRois.ROILabel;
numLabels = numel(roisIds);
% Colormaps
tableCmap = readtable(filenameAtlasColormap);
atlasColormap = [tableCmap.Var3 tableCmap.Var4 tableCmap.Var5];
atlasColormap = atlasColormap./255;
%% PROCESS ALL DATA
% Get a list of preprocessed fMRI files in the directory
fmriSubjects = dir(preprocessedDataPath);
fmriSubjects = fmriSubjects([fmriSubjects(:).isdir]);
% Remove the first two (.,..)
fmriSubjects = fmriSubjects(3:end);
% MAtrix to store all the data:
schaeferSignalsAllSubjects = {};
% Iterate over each preprocessed fMRI file
for i = 1 : numel(fmriSubjects)
    % Read the fMRI data
    niftiFilename = dir(fullfile(preprocessedDataPath, fmriSubjects(i).name, '*.nii'));
    fmriVolume = niftiread(fullfile(preprocessedDataPath, fmriSubjects(i).name, niftiFilename.name));
    fmriInfo = niftiinfo(fullfile(preprocessedDataPath, fmriSubjects(i).name, niftiFilename.name));
    
    % Get the voxel dimensions
    voxelSize_mm3 = fmriInfo.PixelDimensions(1:3);
    tR_sec = fmriInfo.PixelDimensions(4);
    imageSize_voxels = fmriInfo.ImageSize(1:3);
    timePoints = fmriInfo.ImageSize(4);
    
    % Create a matrix to store the signals for the Schaefer atlas
    schaeferSignals = zeros(timePoints, numLabels);
    
    % Iterate over each unique label in the Schaefer atlas
    for label = roisIds'
        % Get the indices corresponding to the current label in the Schaefer atlas
        labelMask = (atlasVolume == label);
        %indicesLabel = find(atlasVolume == label);
        %[indLabelX, indLabelY, indLabelZ] = ind2sub(size(atlasVolume), indicesLabel);
        % Convert mask into 4d:
        labelMask4d = repmat(labelMask, 1, 1, 1, size(fmriVolume,4));
        % Extract the mean signal for the current label
        labelSignal = squeeze(sum(fmriVolume.*labelMask4d, [1 2 3]))./sum(labelMask,'all');
        % Store the label signal in the Schaefer signals matrix
        schaeferSignals(:, label) = labelSignal;
    end
    % Save image to verify:
    SaveOverlayImageWithRois(fmriVolume, atlasVolume, atlasColormap, ...
        [outputPathImages 'meanfMRISchaefer_' fmriSubjects(i).name '.gif']);
            
    % Save them as a .mat and csv:
    save([outputPathSignals suffixROIFilenames fmriSubjects(i).name '.mat'], 'schaeferSignals');
    writematrix(schaeferSignals, [outputPathSignals suffixROIFilenames fmriSubjects(i).name '.txt'],...
        "Delimiter", ',');
end
schaeferSignalsAllSubjects{i} = schaeferSignals;

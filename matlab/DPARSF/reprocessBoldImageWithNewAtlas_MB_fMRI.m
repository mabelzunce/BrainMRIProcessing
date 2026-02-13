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
dataPath = '/home/martin/data_imaging/ADNIdata/ADNI3_Advanced_fMRI/ADNI3_AdvancedFmri/Processed/DPARSF/';
%dataPath = '//home/martin/ADNI3_AdvancedFmri/Processed/DPARSF/';
preprocessedFolder = 'FunImgRWSDCF';
preprocessedDataPath = fullfile(dataPath, preprocessedFolder);
indexScanner = 3; % Siemens=1, GE=2, Philips=3.
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
    % If exists don't process it.
    if ~exist([outputPathSignals suffixROIFilenames fmriSubjects(i).name '.mat'])
        fprintf('Resampling subject %d of %d: %s.\n', i, numel(fmriSubjects), fmriSubjects(i).name);
        tic
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
        schaeferSignalsAllSubjects{i} = schaeferSignals;
        toc
    end
end
%% CREATE A CSV WITH THE MERGE DATA OF THE  DATA
adniCollectionFullFilename = '/data_imaging/ADNIdata/ADNI3_Advanced_fMRI/ADNI3_AdvancedFmri/ADNI3_Advanced_fMRI_6_27_2023.csv';
% Load adni dicom collection csv:
adniCollectionData = readtable(adniCollectionFullFilename);
% Add new fields to the table:
adniCollectionData.PTEDUCAT = NaN(size(adniCollectionData,1), 1);
adniCollectionData.CDRSB = NaN(size(adniCollectionData,1), 1);
adniCollectionData.MMSE = NaN(size(adniCollectionData,1), 1);
adniCollectionData.DIGITSCOR = NaN(size(adniCollectionData,1), 1);
adniCollectionData.MOCA = NaN(size(adniCollectionData,1), 1);
adniCollectionData.Ventricles = NaN(size(adniCollectionData,1), 1);
adniCollectionData.Hippocampus = NaN(size(adniCollectionData,1), 1);
adniCollectionData.WholeBrain = NaN(size(adniCollectionData,1), 1);
adniCollectionData.MidTemp = NaN(size(adniCollectionData,1), 1);
adniCollectionData.ABETA = NaN(size(adniCollectionData,1), 1);
adniCollectionData.TAU = NaN(size(adniCollectionData,1), 1);
adniCollectionData.PTAU = NaN(size(adniCollectionData,1), 1);

% Load the adni merge csv with a summary of all the data
adniMergeFullFilename = '/home/martin/data_imaging/ADNIdata/StudyInfo/Study_Info/ADNIMERGE_12Jul2023.csv';
adniMergeData = readtable(adniMergeFullFilename);
for i = 1 : numel(fmriSubjects)
    % It should appear only once in this data collection.
    indexSubjectInCollection = find(strcmp(fmriSubjects(i).name, adniCollectionData.SubjectID)>0);
    % If more than one, use the last one (assumes that they are all from
    % the same visit):
    indexSubjectInCollection = indexSubjectInCollection(end);
    infoThisImage = adniCollectionData(indexSubjectInCollection,:);
    % If not in merge, it's because has not been updated in the MERGE file,
    % leave it empty
    indexMerge = find(strcmp(fmriSubjects(i).name, adniMergeData.PTID)>0);
    if ~isempty(indexMerge)
        % There can be more than 1, check the closer date.
        % Convert format of mergedate into the csv collection format.
        dateTimeMerge = datetime(adniMergeData(indexMerge,:).EXAMDATE, 'InputFormat', 'yyyy-MM-dd', 'Format', 'dd/MM/yyyy');
        if ~isempty(indexMerge)
            % Add fields for cognitive tests and other metrics.
            % Not all the fields have data for every date, look for the closer
            % date for each field.
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).PTEDUCAT));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            infoThisImage.PTEDUCAT = adniMergeData(indexMerge(indexMin),:).PTEDUCAT;
    
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).CDRSB));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            infoThisImage.CDRSB = adniMergeData(indexMerge(indexMin),:).CDRSB;
    
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).MMSE));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            infoThisImage.MMSE = adniMergeData(indexMerge(indexMin),:).MMSE;
    
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).DIGITSCOR));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            if ~isempty(indexMin)
                infoThisImage.DIGITSCOR = adniMergeData(indexMerge(indexMin),:).DIGITSCOR;
            else
                infoThisImage.DIGITSCOR = NaN;
            end
    
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).MOCA));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            if ~isempty(indexMin)
                infoThisImage.MOCA = adniMergeData(indexMerge(indexMin),:).MOCA;
            else
                infoThisImage.MOCA = NaN;
            end
    
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).Ventricles));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            if ~isempty(indexMin)
                infoThisImage.Ventricles = adniMergeData(indexMerge(indexMin),:).Ventricles;
            else
                infoThisImage.Ventricles = NaN;
            end
    
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).Hippocampus));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            if ~isempty(indexMin)
                infoThisImage.Hippocampus = adniMergeData(indexMerge(indexMin),:).Hippocampus;
            else
                infoThisImage.Hippocampus = NaN;
            end
    
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).WholeBrain));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            if ~isempty(indexMin)
                infoThisImage.WholeBrain = adniMergeData(indexMerge(indexMin),:).WholeBrain;
            else
                infoThisImage.WholeBrain = NaN;
            end
    
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).MidTemp));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            if ~isempty(indexMin)
                infoThisImage.MidTemp = adniMergeData(indexMerge(indexMin),:).MidTemp;
            else
                infoThisImage.MidTemp = NaN;
            end       
            
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).ABETA));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            if ~isempty(indexMin)
                infoThisImage.ABETA = adniMergeData(indexMerge(indexMin),:).ABETA;
            else
                infoThisImage.ABETA = NaN;
            end
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).TAU));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            if ~isempty(indexMin)
                infoThisImage.TAU = adniMergeData(indexMerge(indexMin),:).TAU;
            else
                infoThisImage.TAU = NaN;
            end
    
            indicesNotNan = find(~isnan(adniMergeData(indexMerge,:).PTAU));
            [minValue, indexMin] = min(abs(infoThisImage.StudyDate - dateTimeMerge(indicesNotNan)));
            if ~isempty(indexMin)
                infoThisImage.PTAU = adniMergeData(indexMerge(indexMin),:).PTAU;
            else
                infoThisImage.PTAU = NaN;
            end
        end
    else
        warning(sprintf('Subject %s not found in the ADNIMERGE file.', fmriSubjects(i).name));
        infoThisImage.PTEDUCAT = NaN;
        infoThisImage.CDRSB = NaN;
        infoThisImage.MMSE = NaN;
        infoThisImage.DIGITSCOR = NaN;
        infoThisImage.MOCA = NaN;
        infoThisImage.Ventricles = NaN;
        infoThisImage.Hippocampus = NaN;
        infoThisImage.WholeBrain = NaN;
        infoThisImage.MidTemp = NaN;
        infoThisImage.ABETA = NaN;
        infoThisImage.TAU = NaN;
        infoThisImage.PTAU = NaN;
    end
    infoProcessedData(i,:) = infoThisImage;
    adniCollectionData(indexSubjectInCollection,:) = infoThisImage;
end
writetable(infoProcessedData, fullfile(dataPath, 'SubjctsDataAndTests.csv'))
[path, filenameCollection, ext] = fileparts(adniCollectionFullFilename);
fullFilenameCollectionExtended = fullfile(path, [filenameCollection '_extended' ext]);
writetable(adniCollectionData, fullFilenameCollectionExtended);
%% STATS OF EACH CASE
[values, channels] = hist(categorical(infoProcessedData.ResearchGroup));
[valuesCollection, channelsCollection] = hist(categorical(adniCollectionData.ResearchGroup));
%% CHECK THE IMAGES
for i = 1 : numel(fmriSubjects)
    signals = load([outputPathSignals suffixROIFilenames fmriSubjects(i).name '.mat']);
    fig = check_fMRI_bold_signals(signals.schaeferSignals);
    close all
end
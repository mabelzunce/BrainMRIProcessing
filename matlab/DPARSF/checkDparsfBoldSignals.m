clear all
%close all
dataPartitionPath = '/data/'; %'D:/'
imagingPartitionPath = '/data_imaging/'; %'F:/'
%% ADD PATHS
dpabiPath = [dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/'];
addpath([dataPartitionPath 'UNSAM/Brain/dicm2nii/'])
addpath(genpath(dpabiPath))
addpath([dataPartitionPath 'UNSAM/Brain/spm12/spm12/'])
addpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/DPARSF/'])
%% CONFIG
createPlots = 0;
%% DATA PATHS
dataPath = [dataPartitionPath '/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/DataBaseDrive-20231031T195909Z-001/DataBase/ADNI_fMRI_screening/AAL/'];
dataPath = [imagingPartitionPath '/CovidProject/Estudio/PreprocessedMRI/DPARSF/'];
dataPath = '/home/martin/data_imaging/ADNIdata/fMRI_ADNI2_ADNI3_initial_visit/';
adniCollectionFullFilename = '/home/martin/data/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/fMRI_screening_initalvisit/fMRI_screening_initalvisit_2023_04_27.csv';
adniMergeFullFilename = '/home/martin/data_imaging/ADNIdata/StudyInfo/Study_Info/ADNIMERGE_14Oct2024.csv';
adniImageQualityAdni2Filename = '/home/martin/data_imaging/ADNIdata/StudyInfo/MR_Image_Quality/MAYOADIRL_MRI_IMAGEQC_12_08_15_14Oct2024.csv';
adniImageQualityAdni3Filename = '/home/martin/data_imaging/ADNIdata/StudyInfo/MR_Image_Quality/MAYOADIRL_MRI_QUALITY_ADNI3_14Oct2024.csv';
roiSignalsPath = [dataPath '/Results/ROISignals_FunImgARWSDCF/'];
roiSignalsSubdir = fullfile('ResultsAAL', 'ROISignals_NiftiPreprocessedAllBatchesNorm/');
roiSignalsPath = fullfile(dataPath, 'ResultsAAL', 'ROISignals_NiftiPreprocessedAllBatchesNorm/');

preprocessedFolder = 'NiftiPreprocessedAllBatchesNorm';%'FunImgARWSDCFN';
preprocessedDataPath = fullfile(dataPath, preprocessedFolder);

%roiSignalsPath = '/home/martin/data/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/ADNI3_Advanced_MB_fMRI/AAL/';
csvFilename = [dataPath 'SubjectsDataAndTests.csv'];
outputPath = fullfile(dataPath, 'check/');
if ~isdir(outputPath)
    mkdir(outputPath)
end
%% READ CSV
subjectsData = readtable(csvFilename);
mergeData = readtable(adniMergeFullFilename);
iqAdni2 = readtable(adniImageQualityAdni2Filename);
iqAdni3 = readtable(adniImageQualityAdni3Filename);
%% READ ALL SIGNALS
% Get a list of preprocessed fMRI files in the directory
fmriSubjects = dir(fullfile(roiSignalsPath, '*.mat'));
for i = 1 : numel(fmriSubjects)
    aux = strsplit(fmriSubjects(i).name, '.mat');
    aux = strsplit(aux{1}, 'ROISignals_');
    subjectNames{i} = aux{2};
end
%% CREATE PLOTS
if createPlots
    for i = 1 : numel(subjectNames)
        resultSignals = load(fullfile(roiSignalsPath, fmriSubjects(i).name));
        roiSignals = resultSignals.signals;
        fig = check_fMRI_bold_signals(roiSignals);
        saveas(gcf, fullfile(outputPath, [subjectNames{i}, '.png']));
        close all
    end
end
%% CREATE A CSV WITH THE MERGE DATA OF THE  DATA
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
adniCollectionData.ScannerManufacturer = repmat("Empty",size(adniCollectionData,1), 1);
for i = 1 : size(adniCollectionData, 1)
    substrings = strsplit(adniCollectionData.ImagingProtocol{i},';');
    % Find amnufacturer:
    indexManuf = find(contains(substrings, 'Manufacturer')>0);
    if contains(substrings{indexManuf}, 'Philips','IgnoreCase',true)
        adniCollectionData.ScannerManufacturer(i) = 'Philips';
    elseif contains(substrings{indexManuf}, 'Siemens','IgnoreCase',true)
        adniCollectionData.ScannerManufacturer(i) = 'Siemens';
    elseif contains(substrings{indexManuf}, 'GE','IgnoreCase',true)
        adniCollectionData.ScannerManufacturer(i) = 'GE';
    end
end
% Load the adni merge csv with a summary of all the data
adniMergeData = readtable(adniMergeFullFilename);
adni2IQ = readtable(adniImageQualityAdni2Filename);
adni3IQ = readtable(adniImageQualityAdni3Filename);
% Get images IDs for IQ 2 and 3
for i = 1 : size(adni2IQ,1)
    if ~isempty(adni2IQ.loni_image{i})
        aux = strsplit(adni2IQ.loni_image{i}, 'I');
        imageIdIQ2(i) = str2num(aux{2});
    end
end

% Merge data from different CSVs
for i = 1 : numel(subjectNames)
    % It should appear only once in this data collection.
    indexSubjectInCollection = find(strcmp(subjectNames(i), adniCollectionData.SubjectID)>0);
    % If more than one, use the last one (assumes that they are all from
    % the same visit):
    indexSubjectInCollection = indexSubjectInCollection(end);
    infoThisImage = adniCollectionData(indexSubjectInCollection,:);
    % If not in merge, it's because has not been updated in the MERGE file,
    % leave it empty
    indexMerge = find(strcmp(subjectNames(i), adniMergeData.PTID)>0);

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
        warning(sprintf('Subject %s not found in the ADNIMERGE file.', subjectNames{i}));
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
    adniCollectionData(indexSubjectInCollection,:) = infoThisImage; % The extended without image quality
    % Image quality
    indexIQ = find((infoThisImage.ImageID==imageIdIQ2)>0);
    if ~isempty(indexIQ)
        infoThisImage.series_quality = adni2IQ.series_quality(indexIQ);
        infoThisImage.study_overallpass = adni2IQ.study_overallpass(indexIQ);
        infoThisImage.study_comments = adni2IQ.study_comments(indexIQ);
        infoThisImage.fMRI_motion_display = adni2IQ.fMRI_motion_display(indexIQ);
        infoThisImage.study_medical_abnormalities = adni2IQ.study_medical_abnormalities(indexIQ);
        infoThisImage.fMRI_mean_snr = adni2IQ.fMRI_mean_snr(indexIQ);
        infoThisImage.fMRI_yaw = adni2IQ.fMRI_yaw(indexIQ);
        infoThisImage.fMRI_roll = adni2IQ.fMRI_roll(indexIQ);
        infoThisImage.fMRI_pitch = adni2IQ.fMRI_pitch(indexIQ);
        infoThisImage.fMRI_mmz = adni2IQ.fMRI_mmz(indexIQ);
        infoThisImage.fMRI_mmy = adni2IQ.fMRI_mmy(indexIQ);
        infoThisImage.fMRI_mmx = adni2IQ.fMRI_mmx(indexIQ);
        infoThisImage.fMRI_temporal_slice_order = adni2IQ.fMRI_temporal_slice_order(indexIQ);
    else
        indexIQ = find((infoThisImage.ImageID==adni3IQ.LONI_IMAGE)>0);
        if numel(indexIQ) > 1
            indexIQ = indexIQ(1);
        end
        infoThisImage.series_quality = adni3IQ.SERIES_QUALITY(indexIQ);
        infoThisImage.study_overallpass = adni3IQ.STUDY_OVERALLPASS(indexIQ);
        infoThisImage.study_comments = adni3IQ.STUDY_COMMENTS(indexIQ);
        infoThisImage.study_medical_abnormalities = adni3IQ.STUDY_MEDICAL_ABNORMALITIES(indexIQ);
        infoThisImage.fMRI_mean_snr = NaN;
        infoThisImage.fMRI_motion_display = NaN;
        infoThisImage.fMRI_yaw = NaN;
        infoThisImage.fMRI_roll = NaN;
        infoThisImage.fMRI_pitch = NaN;
        infoThisImage.fMRI_mmz = NaN;
        infoThisImage.fMRI_mmy = NaN;
        infoThisImage.fMRI_mmx = NaN;
        infoThisImage.fMRI_temporal_slice_order = {NaN};
    end
    infoProcessedData(i,:) = infoThisImage;

end
writetable(infoProcessedData, fullfile(dataPath, ['SubjectsDataAndTests' char(datetime('today','Format','y_MM_dd')) '.csv']))
save(fullfile(dataPath, ['SubjctsDataAndTests' char(datetime('today','Format','y_MM_dd')) '.mat']), 'infoProcessedData')
[path, filenameCollection, ext] = fileparts(adniCollectionFullFilename);
fullFilenameCollectionExtended = fullfile(path, [filenameCollection '_extended_' char(datetime('today','Format','y_MM_dd')) ext]);
writetable(adniCollectionData, fullFilenameCollectionExtended);

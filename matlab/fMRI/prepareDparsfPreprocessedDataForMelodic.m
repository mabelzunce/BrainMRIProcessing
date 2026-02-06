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
dataPath = '/home/martin/data_imaging/ADNIdata/Trini/Batch4/2mm/';
dataPath = '/home/martin/data_imaging/ADNIdata/Trini/'
adniCollectionFullFilename = '/home/martin/data/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/fMRI_screening_initalvisit/fMRI_screening_initalvisit_2023_04_27.csv';
adniMergeFullFilename = '/home/martin/data_imaging/ADNIdata/StudyInfo/Study_Info/ADNIMERGE_12Jul2023.csv';
preprocessedSubjectsDataFilename = fullfile(dataPath, 'SubjctsDataAndTests2024_10_17.csv');
preprocessedFolder = 'NiftiPreprocessedAllBatchesNorm';%'FunImgARWSDCFN';
preprocessedDataPath = fullfile(dataPath, preprocessedFolder);
indexScanner = 1; % Siemens=1, GE=2, Philips=3.

%% OUTPUTPATH
outputPath = fullfile(dataPath, 'Melodic');
if ~isfolder(outputPath)
    mkdir(outputPath);
end
filenameSubjectsMelodic = fullfile(outputPath, 'input_files.txt');
%% GET THE LIST OF ALL SCANS
subjectsToExclude = {'114_S_6039', '035_S_6953', '128_S_2002', '031_S_4021', '130_S_5231'};
% Get a list of preprocessed fMRI files in the directory
fmriSubjects = dir(preprocessedDataPath);
fmriSubjects = fmriSubjects([fmriSubjects(:).isdir]);
% Remove the first two (.,..)
fmriSubjects = fmriSubjects(3:end);
% MAtrix to store all the data:
schaeferSignalsAllSubjects = {};
% Iterate over each preprocessed fMRI file
for i = 1 : numel(fmriSubjects)
    % Get names and filenames:
    if sum(strcmp(fmriSubjects(i).name, subjectsToExclude)) == 0
        subjectNames{i} = fmriSubjects(i).name;
        filename = dir(fullfile(preprocessedDataPath, fmriSubjects(i).name, '*.nii.gz'));
        preprocessedImageFilenames{i} = fullfile(preprocessedDataPath, fmriSubjects(i).name, filename.name);
    end
end
%% WRITE A LIST OF FILENAMES TO PROCESS WITH MELODIC
%Choose based on image quality:
subjectsDataInfo = readtable(preprocessedSubjectsDataFilename);
% Get subjects with higher SNR
indicesSelectedSubjects = find(subjectsDataInfo.fMRI_mean_snr > 45);


% Write the input filenames
fid = fopen(filenameSubjectsMelodic,'w');
% First lines to prepare the data:
for i = 1 : numel(indicesSelectedSubjects)
    fprintf(fid, '%s\n', preprocessedImageFilenames{indicesSelectedSubjects(i)});
end
fclose(fid);
%% CREATE MELODIC FILES
%melodic -i input_files.txt -o GroupICA_25 --tr=3.0 --nobet -a concat -m $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask.nii.gz --report --Oall -d 25
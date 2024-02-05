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
addpath('../DPARSF/')
%% DATA PATHs
dataPath = '/home/martin/data/UNSAM/CovidProject/Estudio/PreprocessedMRI/';
niftiPath = '/home/martin/data/UNSAM/CovidProject/Estudio/PreprocessedMRI/Nifti/';
dparsfDataPath = '/home/martin/data/UNSAM/CovidProject/Estudio/PreprocessedMRI/DPARSF/';
outputCheckPath = [dparsfDataPath '/check/'];
if ~isdir(outputCheckPath)
    mkdir(outputCheckPath)
end
indexScanner = 1; % Siemens=1, GE=2, Philips=3.
% Filename with the available MRI data:
filenameMriInfo = [dataPath 'mriInfo_2023_09_05.mat'];
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

overwriteNifti = 1; % If 0 won't process existing Nfiti data in the DPARSF folder
format = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';

%% CHECK DATA
listDir = dir(fmriDparsfPath);
% Get all subjects that are similar:
for i = 1 : numel(listDir)-2
    subjectNames{i} = listDir(i+2).name;
    % Check in MriInfo the data:
    indexThisSubject = find(strcmp(subjectNames{i}, mriInfo.casesToProcess));
    imageSize_voxels(i,:) = mriInfo.fMRI_imageSize_voxels(indexThisSubject,:);
    timePoints(i) = imageSize_voxels(i,4);
    numSlices(i) = imageSize_voxels(i,3);
end

%% CHECK DATA
for i = 1 : numel(subjectNames)
    if exist([roiSignalsDparsfPath suffixROIFilenames subjectNames{i} '.mat'])
        signals = load([roiSignalsDparsfPath suffixROIFilenames subjectNames{i} '.mat']);
        fig = check_fMRI_bold_signals(signals.ROISignals);
        saveas(gca, fullfile(outputCheckPath, subjectNames{i}), 'tif');
        close all
    end
end
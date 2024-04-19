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
%% DEFAULT CONFIG FILE
[ProgramPath, fileN, extn] = fileparts(which('DPARSF.m'));
load([ProgramPath,filesep,'Jobmats',filesep,'Template_CalculateInMNISpace_TraditionalOrder.mat']);


%% DATA PATHs
dataPath = '/home/martin/data/UNSAM/CovidProject/Estudio/PreprocessedMRI/';
niftiPath = '/home/martin/data/UNSAM/CovidProject/Estudio/PreprocessedMRI/Nifti/';
dparsfDataPath = '/home/martin/data/UNSAM/CovidProject/Estudio/PreprocessedMRI/DPARSF_bandpass/';
indexScanner = 1; % Siemens=1, GE=2, Philips=3.
% Filename with the available MRI data:
filenameMriInfo = [dataPath 'mriInfo_2024_03_20.mat'];
%% DATA INFO
mriInfo = load(filenameMriInfo);
%% CONFIG
bandPassFilter = 1; %If not band pass, only high pass
overwritePreprocessedData = 0; % if 0, if the preprocessed singals exist, it doesn't re run it.
%% FOLDERS FOR DPARSF
t1NameNifti = 'T1Img';
fmriNameNifti = 'FunImg';
fmriSliceCorrectedNameNifti = 'FunImgA';
fmriFullyPreprocessedNameNifti = 'FunImgARWSD';
fieldmapBaseNameNifti = 'FieldMap';
fieldmapPhaseNameNifti = 'PhaseDiffImg';
fieldmapMagNameNifti = 'Maginute1Img';

% t1DparsfPath = [dparsfPath '/' t1NameNifti '/'];

fmriDparsfPath = [dparsfDataPath '/' fmriNameNifti '/'];
fmriPreprocessedSignalsDparsfPath = [dparsfDataPath '/Results/ROISignals_FunImgARWSDCF/'];
%% DPARSF CONFIG
overwriteNifti = 0; % If 0 won't process existing Nfiti data in the DPARSF folder
format = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';

% Load config files:
configThisManufacturer = load([dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARSF/DPARSF_config_siemens.mat']);
cfg = configThisManufacturer.Cfg;
listDir = dir(fmriDparsfPath);

% Remove the high cut off of the band pass filter (high pass filter now):
cfg.Filter.ASamplePeriod = cfg.TR;
if bandPassFilter ~= 1
    cfg.Filter.ALowPass_HighCutoff = 0; % 0 mean no cut off.
end

% Voxel size:
cfg.Normalize.VoxSize = [2 2 2];
% Remove first time points:
cfg.IsRemoveFirstTimePoints = 0;
cfg.RemoveFirstTimePoints = 0;

% Restore mean value after regression:
cfg.Covremove.IsAddMeanBack=1;
cfg.Covremove.IsHeadMotionScrubbingRegressors = 0;
auxCSF = cfg.Covremove.CSF;
auxWholeBrain = cfg.Covremove.WholeBrain;
cfg.Covremove = rmfield(cfg.Covremove, 'CSF');
cfg.Covremove.WM = Cfg.Covremove.WM;
cfg.Covremove.CSF = Cfg.Covremove.CSF;
cfg.Covremove.WholeBrain = Cfg.Covremove.WholeBrain;
cfg.Covremove.CSF.IsRemove = auxCSF;
cfg.Covremove.WholeBrain.IsRemove = auxWholeBrain;
cfg.Covremove.WM.IsRemove = cfg.Covremove.WhiteMatter; %YAN Chao-Gan, 20160415. Fixed the bug. Cfg.Covremove.WM.IsRemov = CfgBasic.Covremove.WhiteMatter;

cfg.IsNeedReorientFunImgInteractively = 0; % Do not reorient manually.
cfg.IsNeedReorientT1ImgInteractively = 0;
cfg.Normalize.Timing = 'OnFunctionalData';
cfg.Covremove.Timing = 'AfterNormalize';
cfg.Filter.Timing='AfterNormalize';
cfg.Smooth.Timing = 'OnFunctionalData';
switch cfg.IsNormalize
    case 1
        cfg.IsNeedT1CoregisterToFun=0;
        cfg.IsSegment=0;
        cfg.IsDARTEL=0;
    case 2
        cfg.IsNeedT1CoregisterToFun=1;
        cfg.IsSegment=1;
        cfg.IsDARTEL=0;
    case 3
        cfg.IsNeedT1CoregisterToFun=1;
        cfg.IsSegment=2;
        cfg.IsDARTEL=1;
end

if cfg.IsNormalize==0
    cfg.IsNeedT1CoregisterToFun=0;
    cfg.IsSegment=0;
    cfg.IsDARTEL=0;
    
    cfg.IsNeedReorientFunImgInteractively=0; 
    cfg.IsNeedReorientT1ImgInteractively=0; 
    cfg.IsBet=0; 
    cfg.IsAutoMask=0;  
end
%Turn off VMHC
cfg.IsNormalizeToSymmetricGroupT1Mean=0; %YAN Chao-Gan, 150706
cfg.IsCalVMHC=0;
%cfg.IsAllowGUI = 0;
% File to the AAL atlas:
cfg.CalFC.ROIDef = {[dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/Templates/aal.nii']};
%% CHECK DATA
j = 0;
% Get all subjects:
for i = 1 : numel(listDir)-2
    % Check if preprocessed signals already exists:
    signalsThisSubjectFilname = [fmriPreprocessedSignalsDparsfPath 'ROISignals_' listDir(i+2).name '.mat'];
    if ~exist(signalsThisSubjectFilname) || overwritePreprocessedData
        j = j + 1;
        subjectNames{j} = listDir(i+2).name;
        % Check in MriInfo the data:
        indexThisSubject = find(strcmp(subjectNames{j}, mriInfo.casesToProcess));
        imageSize_voxels(j,:) = mriInfo.fMRI_imageSize_voxels(indexThisSubject,:);
        timePoints(j) = imageSize_voxels(j,4);
        numSlices(j) = imageSize_voxels(j,3);
    end
end

% Group all the subjects with the same numSlices and timePoints:
indicesSameNumSlices = numSlices == mode(numSlices);
if sum(indicesSameNumSlices) ~= numel(subjectNames)
    
    modeNumSlices = mode(numSlices);
    subjectNamesToProcess = subjectNames(indicesSameNumSlices);
    warning(['Not all the images have the same number of slices.This subjects won''t be processed:' subjectNames{~indicesSameNumSlices}]);
end
indicesTimePoints = mode(timePoints);
if sum(indicesTimePoints) ~= numel(subjectNames)
    warning('Not all the images have the same number of time points.')
end
% Set the slice order:
if rem(modeNumSlices,2) == 0
    % Even
    cfg.SliceTiming.SliceOrder = [2:2:modeNumSlices 1:2:modeNumSlices];
else
    % Odd
    cfg.SliceTiming.SliceOrder = [1:2:modeNumSlices 2:2:modeNumSlices];
end
% subjects and paths:
cfg.SubjectID = subjectNamesToProcess;
cfg.SubjectNum = numel(subjectNamesToProcess);
cfg.WorkingDir = dparsfDataPath;
cfg.DataProcessDir = dparsfDataPath;
cfg.TimePoints = 0;
tic
[Error, cfg]=DPARSFA_run(cfg);
toc

clear all
close all

dataPartitionPath = '/data/'; %'D:/'
adniPartitionPath = '/data_imaging/'; %'F:/'
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
%% NIFTI DATA PATH
dataPath = '/home/martin/data_imaging/ADNIdata/ADNI3_Advanced_fMRI/ADNI3_AdvancedFmri/Processed/DPARSFBandPass/';
subjectsInfo = '/home/martin/data_imaging/ADNIdata/ADNI3_Advanced_fMRI/ADNI3_AdvancedFmri/ADNI3_Advanced_fMRI_6_27_2023.csv';
indexScanner = 3; % Siemens=1, GE=2, Philips=3.
%% FOLDERS FOR DPARSF
t1NameNifti = 'T1Img';
fmriNameNifti = 'FunImg';
fmriSliceCorrectedNameNifti = 'FunImgA';
fmriFullyPreprocessedNameNifti = 'FunImgARWSD';
fieldmapBaseNameNifti = 'FieldMap';
fieldmapPhaseNameNifti = 'PhaseDiffImg';
fieldmapMagNameNifti = 'Maginute1Img';

% t1DparsfPath = [dparsfPath '/' t1NameNifti '/'];

fmriDparsfPath = [dataPath '/' fmriNameNifti '/'];

%% PATHS AND FILENAMES
overwriteNifti = 1; % If 0 won't process existing Nfiti data in the DPARSF folder
format = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';

dparsfConfigFilenamesPerManufacturer = {[dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARFS/DPARSF_config_siemens.mat'], ...
    [dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARFS/DPARSF_config_ge.mat'],...
    [dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARFS/DPARSF_config_philips.mat']};

% Load config files:
configThisManufacturer = load(dparsfConfigFilenamesPerManufacturer{indexScanner});
% Update TR for the advanced:
cfg = configThisManufacturer.Cfg;
cfg.TR = 0.607;
% No slice timing for short TRs
cfg.IsSliceTiming = 0;
listDir = dir(fmriDparsfPath);
% Remove the high cut off of the band pass filter (high pass filter now):
cfg.Filter.ASamplePeriod = cfg.TR;
%cfg.Filter.ALowPass_HighCutoff = 0; % 0 mean no cut off.

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
%%
% Get all subjects that are similar:
for i = 1 : 2 %numel(listDir)-2
    subjectNames{i} = listDir(i+2).name;
    niftiFilesThisSubejct = dir([fmriDparsfPath subjectNames{i} '/*.nii*']);
    image = niftiread([fmriDparsfPath subjectNames{i} '/' niftiFilesThisSubejct(1).name]);
    imageSize_voxels(i,:) = size(image);
    timePoints(i) = imageSize_voxels(i,4);
    numSlices(i) = imageSize_voxels(i,3);
end

% Group all the subjects with the same numSlices and timePoints:
indicesSameNumSlices = numSlices==numSlices(1);
if sum(indicesSameNumSlices) ~= numel(subjectNames)
    warning('Not all the images have the same number of slices.')
end
indicesTimePoints = timePoints==timePoints(1);
if sum(indicesTimePoints) ~= numel(subjectNames)
    warning('Not all the images have the same number of time points.')
end
% subjects and paths:
cfg.SubjectID = subjectNames;
cfg.SubjectNum = numel(subjectNames);
cfg.WorkingDir = dataPath;
cfg.DataProcessDir = dataPath;
cfg.TimePoints = 0;
tic
[Error, cfg]=DPARSFA_run(cfg);
toc

%% RERUN
cfg.StartingDirName = 'FunImgRWSD';
cfg.IsRealign = 0;
cfg.IsNormalize = 0;
cfg.IsSmooth = 0;
cfg.IsDetrend = 0;
[Error, cfg]=DPARSFA_run(cfg);
%% CHECK THE IMAGES

tableSubjects = readtable(subjectsInfo);
for i = 1 : numel(subjectNames)
    resultSignals = load([dataPath '/Results/ROISignals_FunImgRWSDCF/ROISignals_' subjectNames{i} '.mat']);
    roiSignals = resultSignals.ROISignals;
    fig = check_fMRI_bold_signals(roiSignals);
end
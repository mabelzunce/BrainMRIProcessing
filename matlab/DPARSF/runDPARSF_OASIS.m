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


%% NIFTI DATA PATH
dataPath = '/home/martin/data_imaging/OASIS/second_batch/OASIS_next_60CN_60AD_2026_05_26_RAW_DATA_120of120-001/OASIS_next_60CN_60AD_2026_05_26_processed/DPARSF/';
indexScanner = 1; % Siemens=1, GE=2, Philips=3. Oasis is Siemens.
%% CONFIG
bandPassFilter = 1; %If not band pass, only high pass

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
startingDir = fmriNameNifti;

%% PATHS AND FILENAMES
overwriteNifti = 1; % If 0 won't process existing Nfiti data in the DPARSF folder
format = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';
% Get data to process:
subjectNames = dir(fmriDparsfPath); subjectNames = {subjectNames(3:end).name};

dparsfConfigFilenamesPerManufacturer = {[dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARSF/DPARSF_config_siemens.mat'], ...
    [dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARSF/DPARSF_config_ge.mat'],...
    [dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARSF/DPARSF_config_philips.mat']};

% Load config files:
configThisManufacturer = load(dparsfConfigFilenamesPerManufacturer{indexScanner});
cfg = configThisManufacturer.Cfg;
% Updaet for Oasis config (hardcoded at the moment):
cfg.TR = 2.2;
cfg.TimePoints = 164;
cfg.SliceTiming.SliceNumber = 36;
cfg.SliceTiming.TR = 2.2;
cfg.SliceTiming.TR = cfg.SliceTiming.TR - cfg.SliceTiming.TR/cfg.SliceTiming.SliceNumber;
cfg.SliceTiming.ReferenceSlice = 18;
cfg.SliceTiming.SliceOrder = [2:2:cfg.SliceTiming.SliceNumber 1:2:cfg.SliceTiming.SliceNumber];
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



%%
% subjects and paths:
batchSize = 10;
numBatches = ceil(numel(subjectNames)./batchSize);
% Get all subjects that are similar:
k=0; % counts the number of subjecdts processed

for i = 1 : numBatches%numel(listDir)-2
    clear subjectNamesThisBatch
    clear imageSize_voxels
    clear timePoints
    clear numSlices
    for j = 1 : min(batchSize, (numel(subjectNames)-k))
        indexDir = (i-1)*batchSize + j;
        subjectNamesThisBatch{j} = subjectNames{indexDir};
    end

    % subjects and paths:
    cfg.SubjectID = subjectNamesThisBatch;
    cfg.SubjectNum = numel(subjectNamesThisBatch);
    cfg.WorkingDir = dataPath;
    cfg.DataProcessDir = dataPath;
    cfg.TimePoints = 0;
    cfg.StartingDirName = startingDir;
    tic
    [Error, cfg_out{indexDir}]=DPARSFA_run(cfg);
    toc
    % Rename realign parameter summary files:
    movefile([dataPath '/RealignParameter/ExcludeSubjectsAccordingToMaxHeadMotion.txt'], [dataPath sprintf('/RealignParameter/ExcludeSubjectsAccordingToMaxHeadMotion_batch%d.txt',i)])
    movefile([dataPath '/RealignParameter/HeadMotion.mat'], [dataPath sprintf('/RealignParameter/HeadMotion_batch%d.mat',i)])
    movefile([dataPath '/RealignParameter/HeadMotion.tsv'], [dataPath sprintf('/RealignParameter/HeadMotion_batch%d.tsv',i)])
    % Free space:
    rmdir([dataPath '/FunImgA/'],'s')
    rmdir([dataPath '/FunImgAR/'],'s')
    rmdir([dataPath '/FunImgARW/'],'s')
    rmdir([dataPath '/FunImgARWS/'],'s')
    rmdir([dataPath '/FunImgARWSD/'],'s')
    rmdir([dataPath '/FunImgARWSDC/'],'s')
    cd(currentPath)
    k = k + min(batchSize, (numel(subjectNames)-k));
end

%% CREATE PLOTS
%% CHECK THE IMAGES
checkPath = fullfile(dataPath, "Check");
mkdir(checkPath);
for i = 1 : numel(subjectNames)
    resultSignals = load([dataPath '/Results/ROISignals_FunImgARWSDCF/ROISignals_' subjectNames{i} '.mat']);
    roiSignals = resultSignals.ROISignals;
    fig = check_fMRI_bold_signals(roiSignals);
    saveas(fig, fullfile(checkPath, [subjectNames{i}, '.png']))
    close all
end
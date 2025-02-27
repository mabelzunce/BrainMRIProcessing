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
dataPath = '/home/martin/data_imaging/ADNIdata/fMRI_first_batch_to_process_2023_04_27/Procesadas_2mm_highpass/';
indexScanner = 3; % Siemens=1, GE=2, Philips=3.
%% CONFIG
bandPassFilter = 0; %If not band pass, only high pass

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

dparsfConfigFilenamesPerManufacturer = {[dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARSF/DPARSF_config_siemens.mat'], ...
    [dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARSF/DPARSF_config_ge.mat'],...
    [dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARSF/DPARSF_config_philips.mat']};

% Load config files:
configThisManufacturer = load(dparsfConfigFilenamesPerManufacturer{indexScanner});
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
%%
% Get all subjects that are similar:
for i = 1 : numel(listDir)-2
    subjectNames{i} = listDir(i+2).name;
    dcmTags = load([fmriDparsfPath subjectNames{i} '/dcmHeaders.mat']);
    niftiFilesThisSubejct = fieldnames(dcmTags.h); % Get all the nifti files converted, it should be only one.
    dcmTagsRsFmri{i} = getfield(dcmTags.h, niftiFilesThisSubejct{1});
    image = niftiread([fmriDparsfPath subjectNames{i} '/' niftiFilesThisSubejct{1}]);
    imageSize_voxels(i,:) = size(image);
    timePoints(i) = imageSize_voxels(i,4);
    numSlices(i) = imageSize_voxels(i,3);
    fMRI_tR(k) = dcmTagsRsFmri{k}.RepetitionTime;
    fMRI_tE(k) = dcmTagsRsFmri{k}.EchoTime;
    if isfield(dcmTagsRsFmri{k}, 'MosaicRefAcqTimes')
        fMRI_sliceAcqTimes{k} = dcmTagsRsFmri{k}.MosaicRefAcqTimes;
        fMRI_fslSliceOrcer{k} = dcmTagsRsFmri{k}.SliceTiming;
        [times, fMRI_sliceOrder{k}] = sort(dcmTagsRsFmri{k}.MosaicRefAcqTimes);
    elseif isfield(dcmTagsRsFmri{k}, 'SliceTiming')
        fMRI_sliceAcqTimes{k} = (0.5 - dcmTagsRsFmri{k}.SliceTiming) * dcmTagsRsFmri{k}.RepetitionTime;
        fMRI_fslSliceOrcer{k} = dcmTagsRsFmri{k}.SliceTiming;
        [times, fMRI_sliceOrder{k}] = sort(fMRI_sliceAcqTimes{k});
    end
    %image = niftiread([niftifMriFilenames{k}]);
    info = niftiinfo(niftiFilename);
    fMRI_imageSize_voxels(k,:) = info.ImageSize;
    fMRI_voxelSize_mm(k,:) = info.PixelDimensions;
    fMRI_inPlanePhaseEncodingDirection{k} = dcmTagsRsFmri{k}.InPlanePhaseEncodingDirection;
    fMRI_unwarpDirection(k,:) = dcmTagsRsFmri{k}.UnwarpDirection;
    if isfield(dcmTagsRsFmri{k}, 'EffectiveEPIEchoSpacing')
        fMRI_effectiveEPIEchoSpacing(k) = dcmTagsRsFmri{k}.EffectiveEPIEchoSpacing;
    end
    manufacturer{k} = dcmTagsRsFmri{k}.Manufacturer;
    model{k} = dcmTagsRsFmri{k}.ManufacturerModelName;
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
%[Error]=DPARSF_run(cfg);
%% CREATE PLOTS
cd(currentPath)
subjectsInfo = '/home/martin/data_imaging/ADNIdata/fMRI_first_batch_to_process_2023_04_27/fMRI_first_batch_to_process.csv';
tableSubjects = readtable(subjectsInfo);
for i = 1 : numel(subjectNames)
    fig = check_fMRI_DPARSF (dataPath ,  tableSubjects , subjectNames{i}, '/Results/ROISignals_FunImgARWSDCF/' , [dataPath 'check/'], 1 );
end
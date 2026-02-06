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
dataPath = '/home/martin/data_imaging/ADNIdata/DPARSF_initialvisit/';
subjectsInfo = [dataPath 'SubjctsDataAndTests2024_10_21.csv'];
subjectsMat = [dataPath 'SubjctsDataAndTests2024_10_21.mat'];

dataSubjects = readtable(subjectsInfo);
dataSubjects = load(subjectsMat);
dataSubjects = dataSubjects.infoProcessedData;
indexScanner = 3; % Siemens=1, GE=2, Philips=3.
%% CONFIG
bandPassFilter = 1; %If not band pass, only high pass
removeTimePoints = 3;
%% FOLDERS FOR DPARSF
t1NameNifti = 'T1Img';
fmriNameNifti = 'Nifti';
fmriNameDparsfNifti = 'FunImg';
fmriSliceCorrectedNameNifti = 'FunImgA';
fmriFullyPreprocessedNameNifti = 'FunImgARWSD';
fieldmapBaseNameNifti = 'FieldMap';
fieldmapPhaseNameNifti = 'PhaseDiffImg';
fieldmapMagNameNifti = 'Maginute1Img';

% t1DparsfPath = [dparsfPath '/' t1NameNifti '/'];
fmriNiftiPath = [dataPath '/' fmriNameNifti '/'];
fmriDparsfPath = [dataPath '/' fmriNameDparsfNifti '/'];

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
listDir = dir(fmriNiftiPath);
for i = 1 : numel(listDir)-2
    subjectNames{i} = listDir(i+2).name;
    niftiFilename = dir(fullfile(fmriNiftiPath, subjectNames{i} , '*.nii.gz'));
    if isempty(niftiFilename)
        niftiFilename = dir(fullfile(fmriNiftiPath, subjectNames{i} , '*.nii'));
    end
    niftiFilename = niftiFilename.name;
    dcmTags = load([fmriNiftiPath subjectNames{i} '/dcmHeaders.mat']);
    niftiFilesThisSubejct = fieldnames(dcmTags.h); % Get all the nifti files converted, it should be only one.
    % Check if nifti filename matches the conversion, won't match if the
    % nifti image has been renamed.
    [~, filenameNoExt, ext] = fileparts(niftiFilename);
    if contains(filenameNoExt, '.') % For compressed .nii.gz
        [~, filenameNoExt, ext] = fileparts(filenameNoExt);
    end
    
    if ~strcmp(filenameNoExt, niftiFilesThisSubejct)
        warning('Nifti filename does not match dcmHeaders from  dicom to nifti converter')
    end
    dcmTagsRsFmri{i} = getfield(dcmTags.h, niftiFilesThisSubejct{1});
    fMRI_tR(i) = dcmTagsRsFmri{i}.RepetitionTime;
    fMRI_tE(i) = dcmTagsRsFmri{i}.EchoTime;
    if isfield(dcmTagsRsFmri{i}, 'MosaicRefAcqTimes')
        fMRI_sliceAcqTimes{i} = dcmTagsRsFmri{i}.MosaicRefAcqTimes;
        fMRI_fslSliceOrcer{i} = dcmTagsRsFmri{i}.SliceTiming;
        [times, fMRI_sliceOrder{i}] = sort(dcmTagsRsFmri{i}.MosaicRefAcqTimes);
    elseif isfield(dcmTagsRsFmri{i}, 'SliceTiming')
        fMRI_sliceAcqTimes{i} = (0.5 - dcmTagsRsFmri{i}.SliceTiming) * dcmTagsRsFmri{i}.RepetitionTime;
        fMRI_fslSliceOrcer{i} = dcmTagsRsFmri{i}.SliceTiming;
        [times, fMRI_sliceOrder{i}] = sort(fMRI_sliceAcqTimes{i});
    end
    %image = niftiread([niftifMriFilenames{i}]);
    info = niftiinfo(fullfile(fmriNiftiPath, subjectNames{i}, niftiFilename));
    fMRI_imageSize_voxels(i,:) = info.ImageSize;
    timePoints(i) = fMRI_imageSize_voxels(i,4);
    numSlices(i) = fMRI_imageSize_voxels(i,3);
    fMRI_voxelSize_mm(i,:) = info.PixelDimensions;
    fMRI_inPlanePhaseEncodingDirection{i} = dcmTagsRsFmri{i}.InPlanePhaseEncodingDirection;
    fMRI_unwarpDirection(i,:) = dcmTagsRsFmri{i}.UnwarpDirection;
    if isfield(dcmTagsRsFmri{i}, 'EffectiveEPIEchoSpacing')
        fMRI_effectiveEPIEchoSpacing(i) = dcmTagsRsFmri{i}.EffectiveEPIEchoSpacing;
    end
    manufacturer{i} = dcmTagsRsFmri{i}.Manufacturer;
    model{i} = dcmTagsRsFmri{i}.ManufacturerModelName;

    % Now get info from the merged CSV with subject data nad image quality
    % metrics:
    indexSubjectsData  = find(strcmp(subjectNames{i}, dataSubjects.SubjectID)>0);
    fmri_sliceOrderAdniIq{i} = dataSubjects.fMRI_temporal_slice_order(indexSubjectsData);
end
%% CHECK SLICE TIMINGSS
indexPhilips = find(strncmpi(manufacturer, 'Philips', numel('Philips'))>0);
indexSiemens = find(strncmpi(manufacturer, 'siemens', numel('siemens'))>0);
indexGE = find(strncmpi(manufacturer, 'GE', numel('GE'))>0);
if (numel(indexGE) + numel(indexSiemens) + numel(indexPhilips)) ~= numel(manufacturer)
    error('Not all scanners classified')
end
indicesUnmatched = [];
indicesDflt = [];
%fMRI_sliceOrderMerged % to emrge both data

% Philips should be empty based on the dixom data:
for i = 1 : numel(indexPhilips)
    if isempty(fMRI_sliceOrder{indexPhilips(i)})
        if ~isnan(fmri_sliceOrderAdniIq{indexPhilips(i)}{1}(1)) && (numel(str2num(fmri_sliceOrderAdniIq{indexPhilips(i)}{1})) == numSlices(indexPhilips(i)))
            fMRI_sliceOrderMerged{indexPhilips(i)} = str2num(fmri_sliceOrderAdniIq{indexPhilips(i)}{1});
        else
            % No data from any of them, GE by default is ascending,
            % odd first:
            fMRI_sliceOrderMerged{indexPhilips(i)}  = [1:2:numSlices(indexPhilips(i)) 2:2:numSlices(indexPhilips(i))];
            indicesDflt = [indicesDflt indexPhilips(i)];
        end
    else
        fMRI_sliceOrderMerged{indexPhilips(i)} = fMRI_sliceOrder{indexPhilips(i)};
        % Check if they are the same:
        if ~isnan(fmri_sliceOrderAdniIq{indexPhilips(i)}{1}) % Is is nan, just leave the on from Dicom
            sliceOrder = str2num(fmri_sliceOrderAdniIq{indexPhilips(i)}{1});
            if  sum(fMRI_sliceOrder{indexPhilips(i)} == sliceOrder) ~= numSlices(indexPhilips(i))
                indicesUnmatched = [indicesUnmatched indexPhilips(i)];
                warning(sprintf('Slicer oder for subject %d %s, manufcaturer %s dont match', indexPhilips(i), subjectNames{indexPhilips(i)}, manufacturer{indexPhilips(i)}))
                
            else
                disp('Matched');
            end
        end
    end
end
% Now with GE, idem strategy
for i = 1 : numel(indexGE)
    if isempty(fMRI_sliceOrder{indexGE(i)})
        if ~isnan(fmri_sliceOrderAdniIq{indexGE(i)}{1})  && (numel(str2num(fmri_sliceOrderAdniIq{indexGE(i)}{1})) == numSlices(indexGE(i)))
            fMRI_sliceOrderMerged{indexGE(i)} = str2num(fmri_sliceOrderAdniIq{indexGE(i)}{1});
        else
            % No data from any of them, default GE
            fMRI_sliceOrderMerged{indexGE(i)}  = [1:2:numSlices(indexGE(i)) 2:2:numSlices(indexGE(i))];
            indicesDflt = [indicesDflt indexGE(i)];
        end
    else
        fMRI_sliceOrderMerged{indexGE(i)} = fMRI_sliceOrder{indexGE(i)};
        % Check if they are the same:
        if ~isnan(fmri_sliceOrderAdniIq{indexGE(i)}{1}) % Is is nan, just leave the on from Dicom
            sliceOrder = str2num(fmri_sliceOrderAdniIq{indexGE(i)}{1});
            if  sum(fMRI_sliceOrder{indexGE(i)} == sliceOrder) ~= numSlices(indexGE(i))
                indicesUnmatched = [indicesUnmatched indexGE(i)];
                warning(sprintf('Slicer oder for subject %d %s, manufcaturer %s dont match', indexGE(i), subjectNames{indexGE(i)}, manufacturer{indexGE(i)}))
            else
                disp('Matched');
            end
        end
    end
end

% Now with Siemens, different default
for i = 1 : numel(indexSiemens)
    if isempty(fMRI_sliceOrder{indexSiemens(i)})
        if ~isnan(fmri_sliceOrderAdniIq{indexSiemens(i)}{1}) && (numel(str2num(fmri_sliceOrderAdniIq{indexSiemens(i)}{1})) == numSlices(indexSiemens(i)))
            fMRI_sliceOrderMerged{indexSiemens(i)} = str2num(fmri_sliceOrderAdniIq{indexSiemens(i)}{1});
        else
            % No data from any of them, Siemens is even first for even
            % slices or odd first for odd slices
            if mod(numSlices(1),2) == 0
                fMRI_sliceOrderMerged{indexSiemens(i)}  = [2:2:numSlices(indexSiemens(i)) 1:2:numSlices(indexSiemens(i))];
            else
                fMRI_sliceOrderMerged{indexSiemens(i)}  = [1:2:numSlices(indexSiemens(i)) 2:2:numSlices(indexSiemens(i))];
            end
            indicesDflt = [indicesDflt indexSiemens(i)];
        end
    else
        % Check if they are the same:
        fMRI_sliceOrderMerged{indexSiemens(i)} = fMRI_sliceOrder{indexSiemens(i)};
        if ~isempty(fmri_sliceOrderAdniIq{indexSiemens(i)})
            if ~isnan(fmri_sliceOrderAdniIq{indexSiemens(i)}{1}) % Is is nan, just leave the on from Dicom
                sliceOrder = str2num(fmri_sliceOrderAdniIq{indexSiemens(i)}{1});
                if  sum(fMRI_sliceOrder{indexSiemens(i)} == sliceOrder) ~= numSlices(indexSiemens(i))
                    indicesUnmatched = [indicesUnmatched indexSiemens(i)];
                    warning(sprintf('Slicer oder for subject %d %s, manufcaturer %s dont match', indexSiemens(i), subjectNames{indexSiemens(i)}, manufacturer{indexSiemens(i)}))
                else
                    disp('Matched');
                end
            end
        end
    end
end

%% RUN DPARSF FOR THE SUBJECTS THAT ARE SIMILAR
% Run DPARSF for each set of patients that have the same maufacturer, slice
% timing and number of time points:
indicesProcessed = [];
indicesToProcess = [1 : numel(subjectNames)];
indicesAll = [1 : numel(subjectNames)];
batch = 0;
%%
while ~isempty(indicesToProcess)
    batch = batch + 1;
    dataPathBatch = fullfile(dataPath, sprintf('batch%d', batch));
    fmriDparsfPath = fullfile(dataPathBatch, fmriNameDparsfNifti);
    indexReference = indicesToProcess(1);
    indexMatched = find(numSlices==numSlices(indexReference) & timePoints==timePoints(indexReference));
    % Check if slice order is matched:
    indicesToRemove = [];
    for i = 1 : numel(indexMatched)
        if sum(fMRI_sliceOrderMerged{indexMatched(i)} == fMRI_sliceOrderMerged{indexReference}) ~= numSlices(i)
            indicesToRemove = [indicesToRemove i];
        end
    end
    % Remove the indices that don't go.
    indexMatched(indicesToRemove) = [];
    % Copy images to process to FunImg
    if ~isdir(fmriDparsfPath)
        mkdir(fmriDparsfPath)
    end
    for i = 1 : numel(indexMatched)
        copyfile(fullfile(fmriNiftiPath, subjectNames{indexMatched(i)}), ...
            fullfile(fmriDparsfPath, subjectNames{indexMatched(i)}))
    end
    % subjects and paths:
    % Reconfigure cfg:
    % Load config files:
    configThisManufacturer = load(dparsfConfigFilenamesPerManufacturer{indexScanner});
    cfg = configThisManufacturer.Cfg;
    
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
    cfg.SliceTiming.SliceOrder = vec(fMRI_sliceOrderMerged{indexReference});
    cfg.SubjectID = subjectNames(indexMatched);
    cfg.SubjectNum = numel(cfg.SubjectID);
    cfg.WorkingDir = dataPathBatch;
    cfg.DataProcessDir = dataPathBatch;
    cfg.TimePoints = timePoints(indexReference);
    cfg.IsRemoveFirstTimePoints = 1;
    cfg.RemoveFirstTimePoints = removeTimePoints;
    cfg.IsDARTEL = 0;
    tic
    [Error, cfg] =DPARSFA_run(cfg);
    toc
    for i = 1 : numel(indexMatched)
        indicesToProcess(indicesToProcess==indexMatched(i)) = [];
        indicesProcessed = [indicesProcessed indexMatched(i)];
    end
    % %[Error]=DPARSF_run(cfg);
    % %% CREATE PLOTS
    % cd(currentPath)
    % tableSubjects = readtable(subjectsInfo);
    % for i = 1 : numel(subjectNames)
    %     fig = check_fMRI_DPARSF (dataPath ,  tableSubjects , subjectNames{i}, '/Results/ROISignals_FunImgARWSDCF/' , [dataPath 'check/'], 1 );
% end
end



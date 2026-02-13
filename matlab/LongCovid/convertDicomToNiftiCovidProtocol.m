clear all
close all

dataPartitionPath = '/data/'; %'D:/'
dataImagingPartitionPath = '/data_imaging/'; %'F:/'
% % Load data:
% loadData = 0;
% Overwrite Nifti:
overwriteNifti = 0; % 0: if nifti conversion exists, just reads the images and headers, 1: covnerts and overwrites again
%% ADD PATHS
addpath([dataPartitionPath 'UNSAM/Brain/dicm2nii/'])
addpath(genpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/']))
addpath([dataPartitionPath 'UNSAM/Brain/spm12/spm12/'])
addpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/DPARSF/'])
addpath('../FSL/')
freesurferSubjectsPath =  [dataPartitionPath '/freesurfer_subjects/'];

%% PATHS AND FILENAMES
niftiExtension = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';
% First study
dicomDataPath = [dataPartitionPath '/UNSAM/CovidProject/Estudio/MRI/'];
preprocessedDataPath = [dataPartitionPath '/UNSAM/CovidProject/Estudio/PreprocessedMRI/'];
preprocessedDataPath = [dataImagingPartitionPath '/CovidProject/Estudio/PreprocessedMRI/'];
% Second study:
dicomDataPath = [dataPartitionPath '/UNSAM/CovidProject2/Imaging/MRI/'];
preprocessedDataPath = [dataImagingPartitionPath '/CovidProject/Estudio2/PreprocessedMRI/'];
niftiDataPath = [preprocessedDataPath '/Nifti/'];
dparsfPath = [preprocessedDataPath '/DPARSF/'];
if ~isdir(niftiDataPath)
    mkdir(niftiDataPath)
end
if ~isdir(dparsfPath)
    mkdir(dparsfPath)
end
%% NAME SEQUENCES
dfltNameT1 = 't1_mprage_1x1x1';
nameFmri = 'funcional_3s_30te_2_4_iso';
nameFieldmappingMag1 = 'gre_field_mapping_2mm_e1';
nameFieldmappingMag2 = 'gre_field_mapping_2mm_e2';
nameFieldmappingPhase = 'gre_field_mapping_2mm_phase';
%% CASES TO PROCESS
casesToProcess = [];%{'CP0002', 'CP0006', 'CP0007', 'CP0008', 'CP0009', 'CP0010'};
%% FOLDERS FOR FREEFURSFER
filenameFreesurferScript = [preprocessedDataPath 'freesurferScript.sh'];
if ~isdir(freesurferSubjectsPath)
    mkdir(freesurferSubjectsPath);
end
freesurferLines = [];
%% FOLDERS FOR DPARSF
t1NameNifti = 'T1Img';
fmriNameNifti = 'FunImg';
fieldmapBaseNameNifti = 'FieldMap';
fieldmapPhaseNameNifti = 'PhaseDiffImg';
fieldmapMag1NameNifti = 'Maginute1Img';
fieldmapMag2NameNifti = 'Maginute2Img';

t1DparsfPath = [dparsfPath '/' t1NameNifti '/'];
if ~isdir(t1DparsfPath)
    mkdir(t1DparsfPath);
end

fmriDparsfPath = [dparsfPath '/' fmriNameNifti '/'];
if ~isdir(fmriDparsfPath)
    mkdir(fmriDparsfPath);
end

fieldmapPhaseDparsfPath = [dparsfPath '/' fieldmapBaseNameNifti '/' fieldmapPhaseNameNifti '/'];
if ~isdir(fieldmapPhaseDparsfPath)
    mkdir(fieldmapPhaseDparsfPath);
end

fieldmapMag1DparsfPath = [dparsfPath '/' fieldmapBaseNameNifti '/' fieldmapMag1NameNifti '/'];
if ~isdir(fieldmapMag1DparsfPath)
    mkdir(fieldmapMag1DparsfPath);
end

fieldmapMag2DparsfPath = [dparsfPath '/' fieldmapBaseNameNifti '/' fieldmapMag2NameNifti '/'];
if ~isdir(fieldmapMag2DparsfPath)
    mkdir(fieldmapMag2DparsfPath);
end
%% FOLDER FOR FSL
fslPreprocessedDataPath = [preprocessedDataPath '/FSL/'];
if ~isdir(fslPreprocessedDataPath)
    mkdir(fslPreprocessedDataPath);
end

%% FOLDER FOR FREESURFER
% This where the freesurfer processed data will be stored, not the
% freesurfer_subject folder. It is used mainly to check if the freesurfer
% data is already available.
freesurferPreprocessedDataPath = [preprocessedDataPath '/Freesurfer/'];
if ~isdir(freesurferPreprocessedDataPath)
    mkdir(freesurferPreprocessedDataPath);
end
freesurferAparcAsegFilename = 'mri/aparc.a2009s+aseg.mgz';
% Counter for subjects to process with freesurfer:
indexFreesurfer = 0;
%% PROCESS EACH CASE
if isempty(casesToProcess)
    dirPaths = dir(dicomDataPath);
    casesToProcess = {dirPaths(3:end).name};
end
%casesToProcess = {'CP0024'};
for i = 1 : numel(casesToProcess)
    %% CONVERT DICOM TO NIFTI
    niftiPathThisSubject = [niftiDataPath casesToProcess{i} '/'];
    dicomPathThisSubject = [dicomDataPath casesToProcess{i}];
    if ~exist(niftiPathThisSubject) || overwriteNifti
        converterOut = dicm2nii(dicomPathThisSubject, niftiPathThisSubject, niftiExtension);
    end    
    dcmTags{i} = load([niftiPathThisSubject dcmHeadersFilename]);
%     %% CREATE FOLDERS FOR MELODIC
%     melodicPath = [pathThisSubject '/Nifti/'];
    % Images for this subject:
    sequencesPerSubject{i} = fieldnames(dcmTags{i}.h);
    %% GET ALL PARAMETERS FOR EACH SEQUENCE
    %% T1
    indexT1 = find(strncmp(sequencesPerSubject{i}, dfltNameT1, numel(dfltNameT1)) > 0);
    % If there are more than 1 t1, use the last one (assumes that one was
    % repeated).
    if numel(indexT1) > 1 % There are two fmris, use the second one.
        % Use only the last one, as this means this sequence was repeated.
        indexT1 = indexT1(end);
    end
    nameT1 = sequencesPerSubject{i}{indexT1};
    dcmTagsT1{i} = getfield(dcmTags{i}.h,nameT1);
    niftiT1Filenames{i} = [niftiPathThisSubject nameT1 '.nii.gz'];
    %fMRI_voxelSize_mm(i,:) = [dcmTagsRsFmri{i}.PixelSpacing' dcmTagsRsFmri{i}.SliceThickness];
    %fMRI_matrixSize_voxels(i,:) = [dcmTagsRsFmri{i}.AcquisitionMatrix(1) dcmTagsRsFmri{i}.AcquisitionMatrix(4) dcmTagsRsFmri{i}.LocationsInAcquisition]; 
    t1_tR(i) = dcmTagsT1{i}.RepetitionTime;
    t1_tE(i) = dcmTagsT1{i}.EchoTime;
    info = niftiinfo([niftiT1Filenames{i}]);
    t1_imageSize_voxels(i,:) = info.ImageSize;
    t1_voxelSize_mm(i,:) = info.PixelDimensions;
    % Get age and sex
    if isfield(dcmTagsT1{i}, 'PatientAge')
        age_years(i) = str2num(dcmTagsT1{i}.PatientAge(1:end-1));
    else
        age_years(i) = 0;
    end
    if isfield(dcmTagsT1{i}, 'PatientSex')
        sex(i) = dcmTagsT1{i}.PatientSex;
    else
        sex(i) = 'N';
    end
    %% fMRI
    indexfMri = find(strncmp(sequencesPerSubject{i}, 'funcional', numel('funcional'))>0);
    if numel(indexfMri) > 1 % There are two fmris, use the second one.
        indexfMriNoMoco = []; % Esclude MoCoSeries.
        for j = 1 : numel(indexfMri)
            namefMri = sequencesPerSubject{i}{indexfMri(j)};
            auxDcmTagsRsFmri = getfield(dcmTags{i}.h,namefMri);
            if ~strcmp(auxDcmTagsRsFmri.SeriesDescription, 'MoCoSeries')
                % If there are Motion Corrected Series
                % Use only the last one, as this means this sequence was repeated.
                indexfMriNoMoco = [indexfMriNoMoco indexfMri(j)];
            end
        end
        % Use the last one
        indexfMri = indexfMriNoMoco(end);
    end
    if ~isempty(indexfMri)
        namefMri = sequencesPerSubject{i}{indexfMri};
        dcmTagsRsFmri{i} = getfield(dcmTags{i}.h,namefMri);
        niftifMriFilenames{i} = [niftiPathThisSubject namefMri '.nii.gz'];
        %fMRI_voxelSize_mm(i,:) = [dcmTagsRsFmri{i}.PixelSpacing' dcmTagsRsFmri{i}.SliceThickness];
        %fMRI_matrixSize_voxels(i,:) = [dcmTagsRsFmri{i}.AcquisitionMatrix(1) dcmTagsRsFmri{i}.AcquisitionMatrix(4) dcmTagsRsFmri{i}.LocationsInAcquisition]; 
        fMRI_tR(i) = dcmTagsRsFmri{i}.RepetitionTime;
        fMRI_tE(i) = dcmTagsRsFmri{i}.EchoTime;
        if isfield(dcmTagsRsFmri{i}, 'MosaicRefAcqTimes')
            fMRI_sliceAcqTimes{i} = dcmTagsRsFmri{i}.MosaicRefAcqTimes;
            fMRI_fslSliceOrcer{i} = dcmTagsRsFmri{i}.SliceTiming;
            [times, fMRI_sliceOrder{i}] = sort(dcmTagsRsFmri{i}.MosaicRefAcqTimes);
        else
            fMRI_sliceAcqTimes{i} = (0.5 - dcmTagsRsFmri{i}.SliceTiming) * dcmTagsRsFmri{i}.RepetitionTime;
            fMRI_fslSliceOrcer{i} = dcmTagsRsFmri{i}.SliceTiming;
            [times, fMRI_sliceOrder{i}] = sort(fMRI_sliceAcqTimes{i});
        end
        %image = niftiread([niftifMriFilenames{i}]);
        info = niftiinfo([niftifMriFilenames{i}]);
        fMRI_imageSize_voxels(i,:) = info.ImageSize;
        fMRI_voxelSize_mm(i,:) = info.PixelDimensions;
        fMRI_inPlanePhaseEncodingDirection(i,:) = dcmTagsRsFmri{i}.InPlanePhaseEncodingDirection;
        fMRI_unwarpDirection(i,:) = dcmTagsRsFmri{i}.UnwarpDirection;
        fMRI_effectiveEPIEchoSpacing(i) = dcmTagsRsFmri{i}.EffectiveEPIEchoSpacing;
    else
        fMRI_tR(i) = NaN;
        fMRI_tE(i) = NaN;
        fMRI_imageSize_voxels(i,:) = NaN;
        fMRI_voxelSize_mm(i,:) = NaN;
        fMRI_inPlanePhaseEncodingDirection(i,:) = NaN;
        fMRI_unwarpDirection(i,:) = NaN;
        fMRI_effectiveEPIEchoSpacing(i) = NaN;
        fMRI_sliceOrder{i} = [];
        fMRI_sliceAcqTimes{i} = [];
        fMRI_fslSliceOrcer{i} = [];
        dcmTagsRsFmri{i} = [];
        niftifMriFilenames{i} = [];
    end

    %% FIELD MAPPING
    indexField = contains(sequencesPerSubject{i}, 'mapping');    
    namesFieldMapping = sequencesPerSubject{i}(indexField);
    if numel(namesFieldMapping) < 3
        warning('Field mapping images of %s are incomplete', casesToProcess{i});
        fieldMapping_unwarpDirection{i}{1} = [];
        fieldMapping_tR(i,1) = NaN;
        fieldMapping_tE(i,1) = NaN;
        fieldMapping_deltaTE(i,1) = NaN;
        fieldMapping_effectiveEPIEchoSpacing(i,1) = NaN;
        fieldMapping_imageSize_voxels(i,1,:) = NaN;
        fieldMapping_voxelSize_mm(i,1,:) = NaN;
    end
    for j = 1 : numel(namesFieldMapping)
        dcmTagsFieldMapping{i}{j} = getfield(dcmTags{i}.h,namesFieldMapping{j});
        niftiFieldMappingFilenames{i}{j} = [niftiPathThisSubject namesFieldMapping{j} '.nii.gz'];
        if contains(namesFieldMapping{j}, 'e1')
            nameFieldmappingMag1 = namesFieldMapping{j};
            fieldMapping_indexMag1(i) = j;
        elseif contains(namesFieldMapping{j}, 'e2')
            nameFieldmappingMag2 = namesFieldMapping{j};
            fieldMapping_indexMag2(i) = j;
        elseif contains(namesFieldMapping{j}, 'phase')
            nameFieldmappingPhase = namesFieldMapping{j};
            fieldMapping_indexPhase(i) = j;
        end
        if contains(namesFieldMapping{j}, 'e1') || contains(namesFieldMapping{j}, 'e2') || contains(namesFieldMapping{j}, 'phase')
            info = niftiinfo([niftiFieldMappingFilenames{i}{j}]);
            fieldMapping_unwarpDirection{i}{j} = dcmTagsFieldMapping{i}{j}.UnwarpDirection;
            fieldMapping_tR(i,j) = dcmTagsFieldMapping{i}{j}.RepetitionTime;
            fieldMapping_tE(i,j) = dcmTagsFieldMapping{i}{j}.EchoTime;
            fieldMapping_deltaTE(i,j) = dcmTagsFieldMapping{i}{j}.deltaTE;
            fieldMapping_effectiveEPIEchoSpacing(i,j) = dcmTagsFieldMapping{i}{j}.EffectiveEPIEchoSpacing;
            fieldMapping_imageSize_voxels(i,j,:) = info.ImageSize;
            fieldMapping_voxelSize_mm(i,j,:) = info.PixelDimensions;            
        end
    end
    %% DTI
    indexDti = contains(sequencesPerSubject{i}, 'TENSOR');    
    namesDti = sequencesPerSubject{i}(indexDti);
    for j = 1 : numel(namesDti)
        dcmTagsDti{i}{j} = getfield(dcmTags{i}.h,namesDti{j});
        niftiDtiFilenames{i}{j} = [niftiPathThisSubject namesDti{j} '.nii.gz'];        
        try
            info = niftiinfo([niftiDtiFilenames{i}{j}]);
            dti_imageSize_voxels(i,j,1:numel(info.ImageSize)) = info.ImageSize; % It can be 3D and 4D.
            dti_voxelSize_mm(i,j,1:numel(info.PixelDimensions)) = info.PixelDimensions;
            dti_unwarpDirection{i}{j} = dcmTagsDti{i}{j}.UnwarpDirection;
            dti_inPlanePhaseEncoding{i}{j} = dcmTagsDti{i}{j}.InPlanePhaseEncodingDirection;
            dti_inPlanePhaseEncodingPositive{i}(j) = dcmTagsDti{i}{j}.CSAImageHeaderInfo.PhaseEncodingDirectionPositive;
            dti_diffusionGradientDirection{i}{j} = dcmTagsDti{i}{j}.CSAImageHeaderInfo.DiffusionGradientDirection;
            dti_tR{i}{j} = dcmTagsDti{i}{j}.RepetitionTime;
            dti_tE1{i}{j} = dcmTagsDti{i}{j}.EchoTime;
        catch exc
            warning('Error when reading the %s image.', namesDti{j});
        end
    end
    %% ASL
    indexAsl = find(contains(sequencesPerSubject{i}, 'asl'));
    if ~isempty(indexAsl)
        nameAsl = sequencesPerSubject{i}{indexAsl};
        dcmTagsASL{i} = getfield(dcmTags{i}.h,nameAsl);
        niftiAslFilenames{i} = [niftiPathThisSubject nameAsl '.nii.gz'];
        info = niftiinfo([niftiAslFilenames{i}]);
        asl_tR(i) = dcmTagsASL{i}.RepetitionTime;
        asl_tE(i) = dcmTagsASL{i}.EchoTime;
        asl_imageSize_voxels(i,:) = info.ImageSize;
        asl_voxelSize_mm(i,:) = info.PixelDimensions;
    end
    %% CREATE FOLDERS FOR THIS SUBJECT DPARSF AND COPY FILES
    t1DparsfPathThisSubject = [t1DparsfPath '/' casesToProcess{i} '/'];
    fmriDparsfPathThisSubject = [fmriDparsfPath '/' casesToProcess{i} '/'];
    fieldmapMag1DparsfPathThisSubject = [fieldmapMag1DparsfPath '/' casesToProcess{i} '/'];
    fieldmapMag2DparsfPathThisSubject = [fieldmapMag2DparsfPath '/' casesToProcess{i} '/'];
    fieldmapPhaseDparsfPathThisSubject = [fieldmapPhaseDparsfPath '/' casesToProcess{i} '/'];
    if exist([niftiPathThisSubject namefMri niftiExtension])
        mkdir(t1DparsfPathThisSubject);
        copyfile([niftiPathThisSubject nameT1 niftiExtension], t1DparsfPathThisSubject);
        mkdir(fmriDparsfPathThisSubject);
        copyfile([niftiPathThisSubject namefMri niftiExtension], fmriDparsfPathThisSubject);
        if exist([niftiPathThisSubject nameFieldmappingMag1 niftiExtension])
            mkdir(fieldmapMag1DparsfPathThisSubject);
            copyfile([niftiPathThisSubject nameFieldmappingMag1 niftiExtension], fieldmapMag1DparsfPathThisSubject);
        end
    end
    if exist([niftiPathThisSubject nameFieldmappingMag2 niftiExtension])
        mkdir(fieldmapMag2DparsfPathThisSubject);
        copyfile([niftiPathThisSubject nameFieldmappingMag2 niftiExtension], fieldmapMag2DparsfPathThisSubject);
    end
    if exist([niftiPathThisSubject nameFieldmappingPhase niftiExtension])
        mkdir(fieldmapPhaseDparsfPathThisSubject);
        copyfile([niftiPathThisSubject nameFieldmappingPhase niftiExtension], fieldmapPhaseDparsfPathThisSubject);
    end
    %% FREESURFER
    % First check if freesurfer is already available.
    freesurferPathThisSubject = [freesurferPreprocessedDataPath '/' casesToProcess{i} '/'];
    if ~exist([freesurferPathThisSubject freesurferAparcAsegFilename])
        % If does not exist, add it later to the script
        indexFreesurfer = indexFreesurfer + 1;
        freesurferSubjectsToProcess{indexFreesurfer} = casesToProcess{i};
        freesurferT1FilenameToProcess{indexFreesurfer} = [niftiPathThisSubject nameT1 niftiExtension];        
    end
end
%% SAVE DATA
save(strcat([preprocessedDataPath 'mriInfo_' ], string(datetime('today','Format','y_MM_dd'))))
%% CHECK IF ALL THE SUBJECTS HAVE ALL THE IMAGES

%% PARAMTERS PER SEQUENCE
tableT1 = table(casesToProcess', t1_voxelSize_mm, t1_imageSize_voxels, t1_tR', t1_tE');
writetable(tableT1, [preprocessedDataPath 't1_parameters' ])
tablefMRI = table(casesToProcess', fMRI_voxelSize_mm, fMRI_imageSize_voxels, fMRI_tR', fMRI_tE', fMRI_effectiveEPIEchoSpacing', fMRI_unwarpDirection);
writetable(tablefMRI, [preprocessedDataPath 'fMRI_parameters' ])
tableFieldMapping = table(casesToProcess', permute(fieldMapping_voxelSize_mm(:,1,:), [1 3 2]), ... % TR and voxel sizes should be the same for all images
    permute(fieldMapping_imageSize_voxels(:,1,:), [1 3 2]), fieldMapping_tR(:,1,:), ...
    fieldMapping_tE(:,1:2,:), fieldMapping_deltaTE(:,1,:), fieldMapping_unwarpDirection');
writetable(tableFieldMapping, [preprocessedDataPath 'fieldMapping_parameters' ]);

%% CREATE FREESURFER SCRIPT
if indexFreesurfer > 0
    % If in windows, save the script and run it later in Ubuntu for windows:
    fid = fopen(filenameFreesurferScript,'w');
    fprintf(fid,'#!/bin/bash \n');
    % First lines to prepare the data:
    for i = 1 : numel(freesurferSubjectsToProcess)
        fprintf(fid, 'recon-all -subject %s -i %s\n', freesurferSubjectsToProcess{i}, freesurferT1FilenameToProcess{i});
    end
    % Second the recon-all lines
    for i = 1 : numel(freesurferSubjectsToProcess)
        fprintf(fid, 'recon-all -subject %s -all -openmp 12 \n', freesurferSubjectsToProcess{i});
    end
    fclose(fid);
else
    disp('Freesurfer is already available for all the subjects.')
end
%% PREPROCESSING fMRI WITH FSL
for i = 1 : numel(casesToProcess)
    if ~isempty(numel(casesToProcess))
        fslPreprocessedDataPathThisSubject = [fslPreprocessedDataPath '/' casesToProcess{i} '/fMRI/'];
        if ~isdir(fslPreprocessedDataPathThisSubject)
            mkdir(fslPreprocessedDataPathThisSubject);
        end
        % Check if fieldmapping available and correct:
        if (fieldMapping_indexMag1(i) ~= 0) && (fieldMapping_indexPhase(i)~=0)
            filenameMag1 = niftiFieldMappingFilenames{i}{fieldMapping_indexMag1(i)};
            filenamePhase = niftiFieldMappingFilenames{i}{fieldMapping_indexPhase(i)};
            outputPathFieldmap = fslPreprocessedDataPathThisSubject;
            dTE = fieldMapping_deltaTE(i,fieldMapping_indexMag1(i));
            outputFilename = [outputPathFieldmap 'fmap_rads.nii.gz']; % default output
            if ~exist(outputFilename)
                output = FslPrepareFieldmap(filenamePhase, filenameMag1, outputPathFieldmap, dTE);
            end
        end
        % Then do the preprocessing with Feat:
        % fMRI_sliceAcqTimes{i} = dcmTagsRsFmri{i}.MosaicRefAcqTimes;
        % fMRI_fslSliceOrcer{i} = dcmTagsRsFmri{i}.SliceTiming;
        % [times, fMRI_sliceOrder{i}] = sort(dcmTagsRsFmri{i}.MosaicRefAcqTimes);
    end
end
%% STRUCTURAL FOR EACH SUBJECT
% Calls robustfov and the bet for the T1_mprage
croppedT1Filename = 't1_robustfov.nii.gz';
betParameters = '-R -f 0.35 -g 0  -o -m'; % Using robustFov and 0.35, owrks better than the centroid that we were using before.
fastParameters = '-t 1 -n 3 -H 0.1 -I 4 -l 20.0 -b';
offset_mm = 50; % use 5 cm above the image centre to start with the centre of the brain.
for i = 1 : numel(casesToProcess)
    fslStructDataPathThisSubject = [fslPreprocessedDataPath '/' casesToProcess{i} '/Structural/'];
    if ~isdir(fslStructDataPathThisSubject)
        mkdir(fslStructDataPathThisSubject);
    end
    
    % Output robustfov:
    niftiT1RobustFovFilenames{i} = [fslStructDataPathThisSubject croppedT1Filename];
    if ~exist(niftiT1RobustFovFilenames{i})
        % 1) Call robustfov
        command = sprintf('robustfov -i %s -r %s', niftiT1Filenames{i}, ...
            niftiT1RobustFovFilenames{i});
        out=system(command);
    end

    % Get the name to check if the brain image is already available and not
    % do it:
    [filepath,name,ext] = fileparts(niftiT1RobustFovFilenames{i});
    if contains(name, '.')
        [p, name, ext2] = fileparts(name);
    else
        ext2 = '';
    end

    % Start above the centre in z, as the protocol includes the neck.
    %centre = round(t1_imageSize_voxels(i,:)./2);
    %centre(3) = centre(3) + offset_mm./t1_voxelSize_mm(3);
    % WE DONT USE THE CTRE PARAMETER NOW.
    fslT1BetFilenames{i} = [fslStructDataPathThisSubject name '_brain' ext2 ext];
    if ~exist(fslT1BetFilenames{i})
        fslT1BetFilenames{i} = FslBet(niftiT1RobustFovFilenames{i}, fslStructDataPathThisSubject, betParameters);
    end

    % Call fast for segmentation:
    fslT1FastFilenames{i} = [fslStructDataPathThisSubject name '_brain_seg' ext2 ext];
    if ~exist(fslT1FastFilenames{i})
        fslT1FastFilenames{i} = FslFast(fslT1BetFilenames{i}, fslStructDataPathThisSubject, fastParameters);
    end
end

%% SIENAX - SECOND METHOD - FOR EACH SUBJECT
clear brainVolumes
clear tableSienax
% In this approach instead of using bet with a custom center, first we
% crop the fov of the image and then we call a standard BET.
betParameters = '-R -f 0.35 -g 0  -o -m'; % Robust size and 0.35.
fastParameters = '';
croppedT1Filename = 't1_cropped.nii.gz';
outputSubdir = 'output';
offset_mm = 50; % use 5 cm above the image centre to start with the centre of the brain.
for i = 1 : numel(casesToProcess)
    fslStructDataPathThisSubject = [fslPreprocessedDataPath '/' casesToProcess{i} '/Structural/'];
    fslSienaxDataPathThisSubject = [fslStructDataPathThisSubject '/Sienax/'];
    if ~isdir(fslSienaxDataPathThisSubject)
        mkdir(fslSienaxDataPathThisSubject);
    end

    % We need to copy the T1 to this directory:
    copyfile(niftiT1RobustFovFilenames{i}, fslSienaxDataPathThisSubject);
    [filepath,name,ext] = fileparts(niftiT1RobustFovFilenames{i});
    if contains(name, '.')
        [p, name, ext2] = fileparts(name);
    else
        ext2 = '';
    end
    inputFilenameRobustFovT1 = [fslSienaxDataPathThisSubject name ext2 ext];
    
    sienaxReportFilename = [fslSienaxDataPathThisSubject '/' outputSubdir '/' 'report.sienax'];
    if ~exist(sienaxReportFilename)
        % 1) Call siena qith robustfovt1
        FslSienax(inputFilenameRobustFovT1, outputSubdir, betParameters, fastParameters);
    end
    [auxBrainVolumes, labels] = FslReadSienaxReport([fslSienaxDataPathThisSubject '/' outputSubdir '/' 'report.sienax']);
    brainVolumes(:,:,i) = auxBrainVolumes;
    delete(inputFilenameRobustFovT1);
end
tableSienax.SubjectNames = casesToProcess';
tableSienax.GreyMatterVolume = squeeze(brainVolumes(1,1,:));
tableSienax.WhiteMatterVolume = squeeze(brainVolumes(2,1,:));
tableSienax.BrainVolume = squeeze(brainVolumes(3,1,:));
tableSienax.GreyMatterUnNormVolume = squeeze(brainVolumes(1,2,:));
tableSienax.WhiteMatterUnNNormVolume = squeeze(brainVolumes(2,2,:));
tableSienax.BrainUnNNormVolume = squeeze(brainVolumes(3,2,:));
tableSienax = struct2table(tableSienax);
writetable(tableSienax, [fslPreprocessedDataPath '/SienaxResults.csv']);
%% PREPARE DATA FOR MELODIC

%% SAVE DATA
save(strcat([preprocessedDataPath 'mriInfoAndProcessing_' ], string(datetime('today','Format','y_MM_dd'))))

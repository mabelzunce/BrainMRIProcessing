clear all
close all

% Overwrite Nifti:
overwriteNifti = 0; % 0: if nifti conversion exists, just reads the images and headers, 1: covnerts and overwrites again

%% PATHS AND FILENAMES
niftiExtension = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';

% First study
%dicomDataPath = [dataPartitionPath '/UNSAM/CovidProject/Estudio/MRI/'];
preprocessedDataPath = ['/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/COVID/Prueba_d2n/'];
% Second study:
niftiDataPath = [preprocessedDataPath '/Nifti/'];
dparsfPath = [preprocessedDataPath '/DPARSF/'];

if ~isdir(niftiDataPath)
    mkdir(niftiDataPath)
end
if ~isdir(dparsfPath)
    mkdir(dparsfPath)
end
%% NAME SEQUENCES -----------
dfltNameT1 = 't1_mprage_1x1x1';
nameFmri = 'funcional_3s_30te_2_4_iso';
%nameFieldmappingMag1 = 'gre_field_mapping_2mm_e1'; % para Estudio1
nameFieldmappingMag1 = 'gre_field_mapping_2mm_estric_e1'; % para Estudio2
%nameFieldmappingMag2 = 'gre_field_mapping_2mm_e2';
nameFieldmappingMag2 = 'gre_field_mapping_2mm_estric_e2';
%nameFieldmappingPhase = 'gre_field_mapping_2mm_phase';
nameFieldmappingPhase = 'gre_field_mapping_2mm_estric_phase';
%% CASES TO PROCESS
casesToProcess = [];
%% FOLDERS FOR DPARSF
t1NameNifti = 'T1Img';
fmriNameNifti = 'FunImg';
fieldmapBaseNameNifti = 'FunFieldMap';
fieldmapPhaseNameNifti = 'PhaseDiff';
fieldmapMag1NameNifti = 'Magnitude1';
fieldmapMag2NameNifti = 'Magnitude2';

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

%% PROCESS EACH CASE
if isempty(casesToProcess)
    %dirPaths = dir(dicomDataPath);----------------
    dirPaths = dir(niftiDataPath);
    casesToProcess = {dirPaths(3:end).name};
end
%%
for i = 1 : numel(casesToProcess)
    niftiPathThisSubject = [niftiDataPath casesToProcess{i} '/'];
    if ~exist(niftiPathThisSubject) || overwriteNifti
        converterOut = dicm2nii(dicomPathThisSubject, niftiPathThisSubject, niftiExtension);
    end    
    dcmTags{i} = load([niftiPathThisSubject dcmHeadersFilename]);
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
        jsonStruct = struct();
    jsonStruct.PhaseEncodingDirection = dcmTagsRsFmri{i}.UnwarpDirection;
    jsonStruct.EffectiveEchoSpacing = dcmTagsRsFmri{i}.EffectiveEPIEchoSpacing/1000;
    if strcmp(dcmTagsRsFmri{i}.InPlanePhaseEncodingDirection, 'COL')
        ReconMatrixPE = dcmTagsRsFmri{i}.AcquisitionMatrix(4);
    elseif strcmp(dcmTagsRsFmri{i}.InPlanePhaseEncodingDirection, 'ROW')
        ReconMatrixPE = dcmTagsRsFmri{i}.AcquisitionMatrix(1);
    else
        ReconMatrixPE = NaN; 
    end
    jsonStruct.ReconMatrixPE = ReconMatrixPE;
    jsonPath = fullfile(niftiPathThisSubject, 'FunImg.json');
    jsonText = jsonencode(jsonStruct, 'PrettyPrint', true);
    fid = fopen(jsonPath, 'w');
    fwrite(fid, jsonText);
    fclose(fid);
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
        warning('Funcional json of %s was not created', casesToProcess{i});
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
            EchoTime1 = dcmTagsFieldMapping{i}{j}.EchoTime;
        elseif contains(namesFieldMapping{j}, 'e2')
            nameFieldmappingMag2 = namesFieldMapping{j};
            fieldMapping_indexMag2(i) = j;
            EchoTime2 = dcmTagsFieldMapping{i}{j}.EchoTime;
        elseif contains(namesFieldMapping{j}, 'phase')
            nameFieldmappingPhase = namesFieldMapping{j};
            fieldMapping_indexPhase(i) = j;
            EchoTime2 = dcmTagsFieldMapping{i}{j}.SecondEchoTime;
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
    if ~isnan(EchoTime1) || ~isnan(EchoTime2) 
        jsonStruct = struct();
        jsonStruct.EchoTime1 = EchoTime1/1000;
        jsonStruct.EchoTime2 = EchoTime2/1000;
        jsonText = jsonencode(jsonStruct, 'PrettyPrint', true);
        fid = fopen(fullfile(niftiPathThisSubject, 'PhaseDiff.json'), 'w');
        fwrite(fid, jsonText);
        fclose(fid);
    else
        warning('Field mapping json of %s was not created', casesToProcess{i});
    end
    EchoTime1 = NaN;
    EchoTime2 = NaN;
   
    %% CREATE FOLDERS FOR THIS SUBJECT DPARSF AND COPY FILES
    t1DparsfPathThisSubject = [t1DparsfPath '/' casesToProcess{i} '/'];
    fmriDparsfPathThisSubject = [fmriDparsfPath '/' casesToProcess{i} '/'];
    fieldmapMag1DparsfPathThisSubject = [fieldmapMag1DparsfPath '/' casesToProcess{i} '/'];
    fieldmapMag2DparsfPathThisSubject = [fieldmapMag2DparsfPath '/' casesToProcess{i} '/'];
    fieldmapPhaseDparsfPathThisSubject = [fieldmapPhaseDparsfPath '/' casesToProcess{i} '/'];
    newT1Name = ['co' nameT1 niftiExtension];
    if exist([niftiPathThisSubject namefMri niftiExtension])
        mkdir(t1DparsfPathThisSubject);
        copyfile([niftiPathThisSubject nameT1 niftiExtension], [t1DparsfPathThisSubject, newT1Name]);
        mkdir(fmriDparsfPathThisSubject);
        copyfile([niftiPathThisSubject namefMri niftiExtension], fmriDparsfPathThisSubject);
        gunzip([fmriDparsfPathThisSubject namefMri niftiExtension]);
        delete([fmriDparsfPathThisSubject namefMri niftiExtension]);
        if exist([niftiPathThisSubject 'FunImg.json'])
        copyfile([niftiPathThisSubject 'FunImg.json'], fmriDparsfPathThisSubject);
        delete([niftiPathThisSubject 'FunImg.json']);
        end
        if exist([niftiPathThisSubject nameFieldmappingMag1 niftiExtension])
            mkdir(fieldmapMag1DparsfPathThisSubject);
            copyfile([niftiPathThisSubject nameFieldmappingMag1 niftiExtension], fieldmapMag1DparsfPathThisSubject);
            gunzip([fieldmapMag1DparsfPathThisSubject nameFieldmappingMag1 niftiExtension]);
            delete([fieldmapMag1DparsfPathThisSubject nameFieldmappingMag1 niftiExtension]);
        end
    end
    if exist([niftiPathThisSubject nameFieldmappingMag2 niftiExtension])
        mkdir(fieldmapMag2DparsfPathThisSubject);
        copyfile([niftiPathThisSubject nameFieldmappingMag2 niftiExtension], fieldmapMag2DparsfPathThisSubject);
    end
    if exist([niftiPathThisSubject nameFieldmappingPhase niftiExtension])
        mkdir(fieldmapPhaseDparsfPathThisSubject);
        copyfile([niftiPathThisSubject nameFieldmappingPhase niftiExtension], fieldmapPhaseDparsfPathThisSubject);
        gunzip([fieldmapPhaseDparsfPathThisSubject nameFieldmappingPhase niftiExtension]);
        delete([fieldmapPhaseDparsfPathThisSubject nameFieldmappingPhase niftiExtension]);
        if exist([niftiPathThisSubject 'PhaseDiff.json'])
            copyfile([niftiPathThisSubject 'PhaseDiff.json'], fieldmapPhaseDparsfPathThisSubject);
            delete([niftiPathThisSubject 'PhaseDiff.json']);
        end
    end
end
%% SAVE DATA
save(strcat([preprocessedDataPath 'mriInfo_' ], string(datetime('today','Format','y_MM_dd'))))
%% SAVE DATA
save(strcat([preprocessedDataPath 'mriInfoAndProcessing_' ], string(datetime('today','Format','y_MM_dd'))))

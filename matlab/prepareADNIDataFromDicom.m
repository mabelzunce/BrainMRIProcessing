clear all
close all

dataPartitionPath = '/data/'; %'D:/'
adniPartitionPath = '/data_imaging/'; %'F:/'
% Load data:
loadData = 0;
%% ADD PATHS
addpath([dataPartitionPath 'UNSAM/Brain/dicm2nii/'])
addpath(genpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/']))
addpath([dataPartitionPath 'UNSAM/Brain/spm12/spm12/'])
addpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/DPARSF/'])

%% PATHS AND FILENAMES
overwriteNifti = 1; % If 0 won't process existing Nfiti data in the DPARSF folder
format = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';
adniPath = [adniPartitionPath '/ADNIdata/ADNI3_Advanced_fMRI/ADNI3_AdvancedFmri/'];
dataPath = [adniPath '/ADNI/'];

dparsfConfigFilenamesPerManufacturer = {[dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARFS/DPARSF_config_siemens.mat'], ...
    [dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARFS/DPARSF_config_ge.mat'],...
    [dataPartitionPath 'UNSAM/CEMSC3/ProcesamientoADNI/DPARFS/DPARSF_config_philips.mat']};

outputPath = [adniPath '/Processed/'];
dparsfPath = [outputPath '/DPARSF/'];
niftiPath = [outputPath '/Nifti/'];
if ~isdir(dparsfPath)
    mkdir(dparsfPath)
end
if ~isdir(niftiPath)
    mkdir(niftiPath)
end
%% ADNI DATABASE COLLECTION
filenameADNICollection = 'ADNI3_AdvancedFmri_6_27_2023.csv';
nameRsFmriInCollectionDataBase = {'Axial rsfMRI (Eyes Open)', 'Axial fcMRI', 'Axial - Advanced fMRI'}; %{'Axial rsfMRI (Eyes Open)', 'Axial fcMRI'};
nameRsFmriAdvancedInCollectionDataBase = {'Axial MB rsfMRI'};%'Axial MB rsfMRI (Eyes Open)';
% Load adni dicom collection csv:
adniCollectionData = readtable([adniPath filenameADNICollection]);
%% LOAD PREVIOUS RUN
% If didn't finish processing, load the workspace and to cntinue from the
% PROCESS EACH CASE cell changing the i index.
if loadData
    RunDparsf(dparsfConfigFilenamesPerManufacturer, [outputPath 'workspaceProcessingScript.mat'])
    exit
end
%% NAME SEQUENCES
nameT1 = 't1_mprage_1x1x1';
nameRsFmri = {'Axial_rsfMRI', 'Axial_fcMRI'};
nameRsFmriAdvanced = 'Axial_MB_rsfMRI';
% nameFieldmappingMag1 = 'gre_field_mapping_2mm_e1';
% nameFieldmappingMag2 = 'gre_field_mapping_2mm_e2';
% nameFieldmappingPhase = 'gre_field_mapping_2mm_phase';
%% CASES TO PROCESS
casesToProcess = [];
j = 1;
if isempty(casesToProcess)
    casesToProcess = {};
    casesToProcessAdvanced = {};
    pathRsFmri = {};
    pathRsFmriAdvanced = {};
    dirPath = dir(dataPath);
    for i = 3 : numel(dirPath)
        if(dirPath(i).isdir) 
            dirSubpath = dir([dataPath dirPath(i).name]);
            disp({dirSubpath(:).name});
            for j = 1 : numel(nameRsFmri)
                indicesRsFmri = strncmpi(nameRsFmri{j}, {dirSubpath(:).name}, numel(nameRsFmri{j}));
                if(sum(indicesRsFmri))
                    break;
                end
            end
            indiceRsFmriAdvanced = strncmp(nameRsFmriAdvanced, {dirSubpath(:).name}, numel(nameRsFmriAdvanced));
            if sum(indicesRsFmri) > 0
                casesToProcess{end+1} = dirPath(i).name;
                index = numel(casesToProcess);
                % indices if more than one subfolder was found:
                indices = find(indicesRsFmri > 0);
                subfolderScanDates = cell(0);
                indicesSubPaths = [];
                for j = 1 : numel(indices)
                    tempPath = [dataPath dirPath(i).name '/' dirSubpath(indices(j)).name '/'];
                    % If more than one date, get the last one:
                    tempDir = dir(tempPath);
                    subfolderScanDates = [subfolderScanDates {tempDir(:).name}];
                    indicesSubPaths = [indicesSubPaths repmat(indices(j), 1, numel(subfolderScanDates))];
                end
                [s, sortedIndices] = sort(subfolderScanDates);
                subfolderLatestScan = subfolderScanDates{sortedIndices(end)};
                % Create datetime struct:
                dateTimeThisScan = datetime(subfolderLatestScan(1:10),'Format','y-MM-dd');
                % Continue with subfolders
                tempPath = [dataPath dirPath(i).name '/' dirSubpath(indicesSubPaths(sortedIndices(end))).name '/' subfolderLatestScan '/'];
                tempDir = dir(tempPath);
                subfolderStudy = sort({tempDir(:).name});
                subfolderStudy = subfolderStudy{end}; % For resting state should be only one.
                pathRsFmri{end+1} = [tempPath '/' subfolderStudy '/'];
                % Look for the image in the collection database:
                %% MATCH IMAGE DATA WITH IMAGE COLLECTION DATA
                indicesMatchedName = find((strcmpi(casesToProcess{end}, adniCollectionData.Subject)) > 0); 
                % Search for fmri sequence:
                for j = 1 : numel(nameRsFmriInCollectionDataBase)
                    indicesMatchedNameAndSequence = find(strncmpi(nameRsFmriInCollectionDataBase{j}, adniCollectionData.Description(indicesMatchedName), numel(nameRsFmriInCollectionDataBase{j}))>0);
                    if ~isempty(indicesMatchedNameAndSequence)
                        break;
                    end
                end
                % Get the same scan as in the file:
                indexDate = find(adniCollectionData.AcqDate(indicesMatchedName(indicesMatchedNameAndSequence)) == dateTimeThisScan);
                if ~isempty(indexDate)
                    indexThisSubjectAndDate = indicesMatchedName(indicesMatchedNameAndSequence(indexDate));
                    if numel(indexThisSubjectAndDate) > 0 % If more than one with the same date and sequence (it can happen if downloaded more than once).
                        indexThisSubjectAndDate = indexThisSubjectAndDate(1);
                    end
                    subjectName{index} = adniCollectionData.Subject{indexThisSubjectAndDate};
                    group{index} = adniCollectionData.Group{indexThisSubjectAndDate};
                    age_years(index) = adniCollectionData.Age(indexThisSubjectAndDate);
                    sex(index) = adniCollectionData.Sex{indexThisSubjectAndDate};
                    visit{index} = adniCollectionData.Visit{indexThisSubjectAndDate};
                    date(index) = adniCollectionData.AcqDate(indexThisSubjectAndDate);
                    if date(index) ~= dateTimeThisScan
                        warning('Scan date of the fMRI does not match the database.')
                    end
                else
                    warning('Scan was not found in the collection database.')
                end
            end
            if sum(indiceRsFmriAdvanced) > 0
                casesToProcessAdvanced{end+1} = dirPath(i).name;
                index = numel(casesToProcessAdvanced);
                tempPath = [dataPath dirPath(i).name '/' dirSubpath(indiceRsFmriAdvanced).name '/'];
                % If more than one date, get the last one:
                tempDir = dir(tempPath);
                subfolderLatestScan = sort({tempDir(:).name});
                subfolderLatestScan = subfolderLatestScan{end};
                % Create datetime struct:
                dateTimeThisScan = datetime(subfolderLatestScan(1:10),'Format','y-MM-dd');
                % Continue with subfolders
                tempPath = [dataPath dirPath(i).name '/' dirSubpath(indiceRsFmriAdvanced).name '/' subfolderLatestScan '/'];
                tempDir = dir(tempPath);
                subfolderStudy = sort({tempDir(:).name});
                subfolderStudy = subfolderStudy{end}; % For resting state should be only one.
                pathRsFmriAdvanced{end+1} = [tempPath '/' subfolderStudy '/'];
                %% MATCH IMAGE DATA WITH IMAGE COLLECTION DATA
                indicesMatchedName = find((strcmpi(casesToProcessAdvanced{end}, adniCollectionData.Subject))); 
                % Search for fmri sequence:
                for j = 1 : numel(nameRsFmriAdvancedInCollectionDataBase)
                    indicesMatchedNameAndSequence = find(strncmpi(nameRsFmriAdvancedInCollectionDataBase{j}, adniCollectionData.Description(indicesMatchedName), numel(nameRsFmriAdvancedInCollectionDataBase{j}))>0);
                    if ~isempty(indicesMatchedNameAndSequence)
                        break;
                    end
                end
                 % Get the same scan as in the file:
                indexDate = find(adniCollectionData.AcqDate(indicesMatchedName(indicesMatchedNameAndSequence)) == dateTimeThisScan);
                if ~isempty(indexDate)
                    indexThisSubjectAndDate = indicesMatchedName(indicesMatchedNameAndSequence(indexDate));
                     if numel(indexThisSubjectAndDate) > 0 % If more than one with the same date and sequence (it can happen if downloaded more than once).
                        indexThisSubjectAndDate = indexThisSubjectAndDate(1);
                    end
                    subjectName{index} = adniCollectionData.Subject{indexThisSubjectAndDate};
                    group{index} = adniCollectionData.Group{indexThisSubjectAndDate};
                    age_years(index) = adniCollectionData.Age(indexThisSubjectAndDate);
                    sex(index) = adniCollectionData.Sex{indexThisSubjectAndDate};
                    visit{index} = adniCollectionData.Visit{indexThisSubjectAndDate};
                    date(index) = adniCollectionData.AcqDate(indexThisSubjectAndDate);
                    if date(index) ~= dateTimeThisScan
                        warning('Scan date of the MB fMRI does not match the database.')
                    end
                else
                    warning('Scan was not found in the collection database.')
                end
            end
        end
    end
end
%% FOLDERS FOR DPARSF
t1NameNifti = 'T1Img';
fmriNameNifti = 'FunImg';
fmriSliceCorrectedNameNifti = 'FunImgA';
fmriFullyPreprocessedNameNifti = 'FunImgARWSD';
fieldmapBaseNameNifti = 'FieldMap';
fieldmapPhaseNameNifti = 'PhaseDiffImg';
fieldmapMagNameNifti = 'Maginute1Img';

% t1DparsfPath = [dparsfPath '/' t1NameNifti '/'];
% if ~isdir(t1DparsfPath)
%     mkdir(t1DparsfPath);
% end

fmriDparsfPath = [dparsfPath '/' fmriNameNifti '/'];
if ~isdir(fmriDparsfPath)
    mkdir(fmriDparsfPath);
end

% fieldmapPhaseDparsfPath = [dparsfPath '/' fieldmapBaseNameNifti '/' fieldmapPhaseNameNifti '/'];
% if ~isdir(fieldmapPhaseDparsfPath)
%     mkdir(fieldmapPhaseDparsfPath);
% end
% 
% fieldmapMagDparsfPath = [dparsfPath '/' fieldmapBaseNameNifti '/' fieldmapMagNameNifti '/'];
% if ~isdir(fieldmapMagDparsfPath)
%     mkdir(fieldmapMagDparsfPath);
% end

%% PROCESS EACH CASE
for i = 1 : numel(casesToProcess)    
    %% CONVERT DICOM TO NIFTI
    pathThisSubject = [dataPath casesToProcess{i}];
    dicomPath = [pathThisSubject '/DICOM/'];
    niftiPathThisSubject = [niftiPath casesToProcess{i} '/'];
    if ~exist(niftiPathThisSubject) || overwriteNifti
        dcmTagsTemp = dicm2nii(pathRsFmri{i}, niftiPathThisSubject, format);
    end
    %% GET DATA FROM DICOM
    if exist([niftiPathThisSubject dcmHeadersFilename])
        dcmTags{i} = load([niftiPathThisSubject dcmHeadersFilename]);
        niftiFilesThisSubejct = fieldnames(dcmTags{i}.h); % Get all the nifti files converted, it should be only one.
        niftiRsFmriFilenames{i} = niftiFilesThisSubejct{1};
        dcmTagsRsFmri{i} = getfield(dcmTags{i}.h, niftiRsFmriFilenames{i});
        voxelSize_mm(i,:) = [dcmTagsRsFmri{i}.PixelSpacing' dcmTagsRsFmri{i}.SliceThickness];
        imageSize_mm(i,:) = [dcmTagsRsFmri{i}.Rows dcmTagsRsFmri{i}.Columns dcmTagsRsFmri{i}.LocationsInAcquisition]; 
        tRs(i) = dcmTagsRsFmri{i}.RepetitionTime;
        scannerManufacturer{i} = dcmTagsRsFmri{i}.Manufacturer;
        scannerModel{i} = dcmTagsRsFmri{i}.ManufacturerModelName;
        
        if isfield(dcmTagsRsFmri{i}, 'PatientAge')
            dcmAge_years(i) = str2num(dcmTagsRsFmri{i}.PatientAge(1:end-1));
        else
            dcmAge_years(i) = 0;
        end
        if isfield(dcmTagsRsFmri{i}, 'PatientSex')
            dcmSex(i) = dcmTagsRsFmri{i}.PatientSex;
        else
            dcmSex(i) = 'N';
        end
        image = niftiread([niftiPathThisSubject niftiRsFmriFilenames{i}]);
        imageRealSize_voxels(i,:) = size(image);
        timePoints(i) = size(image,4);
        if strncmpi(scannerManufacturer{i}, 'Philips', numel('Philips'))
            % Philips does not include slice order, by defaults:
            numSlices = imageRealSize_voxels(i,3);
            sliceOrder{i} = [1:2:numSlices 2:2:numSlices];
        else
            [sortedTimes, sortedIndices] = sort(dcmTagsRsFmri{i}.SliceTiming);
            sliceOrder{i} = sortedIndices(end:-1:1)';
        end
        %% CREATE FOLDERS FOR THIS SUBJECT DPARSF AND COPY FILES
        %t1DparsfPathThisSubject = [t1DparsfPath '/' casesToProcess{i} '/'];
        %% CHECK IF DATA ALREADY THERE
        fmriDparsfPathThisSubject = [fmriDparsfPath '/' casesToProcess{i} '/'];
        if ~exist([fmriDparsfPathThisSubject niftiRsFmriFilenames{i} format]) || overwriteNifti
            %fieldmapMagDparsfPathThisSubject = [fieldmapMagDparsfPath '/' casesToProcess{i} '/'];
            %fieldmapPhaseDparsfPathThisSubject = [fieldmapPhaseDparsfPath '/' casesToProcess{i} '/'];
            %mkdir(t1DparsfPathThisSubject);
            %copyfile([niftiPath nameT1 format], t1DparsfPathThisSubject);
            mkdir(fmriDparsfPathThisSubject);
            copyfile([niftiPathThisSubject niftiRsFmriFilenames{i} format], fmriDparsfPathThisSubject);
            %mkdir(fieldmapMagDparsfPathThisSubject);
            %copyfile([niftiPath nameFieldmappingMag1 format], fieldmapMagDparsfPathThisSubject);
            %mkdir(fieldmapPhaseDparsfPathThisSubject);
            %copyfile([niftiPath nameFieldmappingPhase format], fieldmapPhaseDparsfPathThisSubject);
        end
    else
        warning([niftiPathThisSubject ': Could not be converted to Nifti']);
    end
end
%% READ SOME ADDITIONAL DATA FROM THE HEADERS
for i = 1 : numel(casesToProcess)  
    %  spm_ms = (0.5 - s.SliceTiming) * s.RepetitionTime;
    %  [~, spm_order] = sort(-s.SliceTiming);
    if isfield(dcmTagsRsFmri{i}, 'SliceTiming')
        sliceTimingSPM{i} = (0.5 - dcmTagsRsFmri{i}.SliceTiming)*tRs(i);
        [~, spm_order{i}] = sort(sliceTimingSPM{i}); % Check if the same as sliceOrder
    end
end
%%
% for i = 1 : numel(casesToProcess)    
%     %% CONVERT DICOM TO NIFTI
%     pathThisSubject = [dataPath casesToProcess{i}];
%     dicomPath = [pathThisSubject '/DICOM/'];
%     niftiPathThisSubject = [niftiPath casesToProcess{i} '/'];
%     %% GET DATA FROM DICOM
%     dcmTags{i} = load([niftiPathThisSubject dcmHeadersFilename]);
%     niftiFilesThisSubejct = fieldnames(dcmTags{i}.h); % Get all the nifti files converted, it should be only one.
%     niftiRsFmriFilenames{i} = niftiFilesThisSubejct{1};
% end
%% FIELDMAPPING CORRECTION
% h.gre_field_mapping_2mm_e1.deltaTE
% fprintf('Parameters for fieldmap correction. TE1:%.1fms, TE2:%.1fms, dTE:%.1fms\n', h.gre_field_mapping_2mm_e1.EchoTime, ...
%     h.gre_field_mapping_2mm_e2.EchoTime, h.gre_field_mapping_2mm_e1.deltaTE);

%% SAVE DATA
save([outputPath 'workspaceProcessingScript']);
%% CREATE SPECFIC DPARSF FOLDERS FOR GROUPS
dparsfPerGroupPath = [outputPath '/DPARSFgroups/'];
% Divide them by scanner as they have different parameters.
% SIEMENS, same image size:
groups = {'AD', 'CN', 'MCI', 'EMCI', 'LMCI'};
scanners = {'SI', 'GE','Philips'};
for s = 1 : numel(scanners)
    for g = 1 : numel(groups)
        dparsfThisGroupPath = [dparsfPerGroupPath groups{g} '_' scanners{s} '/' fmriNameNifti '/'];
        mkdir(dparsfThisGroupPath);        
        indicesSubgroup = find( (strncmp(scanners{s}, scannerManufacturer,2) & strcmp(group, groups{g})) >0);
        spm_orderSubgroup{s,g} =  spm_order(indicesSubgroup);
        imageRealSizePerGroup_voxels{s,g} = imageRealSize_voxels(indicesSubgroup,:);
        timePointsPerGroup{s,g} = timePoints(indicesSubgroup);
        for i = indicesSubgroup
            niftiPathThisSubject = [niftiPath casesToProcess{i} '/'];
            fmriDparsfPathThisSubject = [dparsfThisGroupPath '/' casesToProcess{i} '/'];
            mkdir(fmriDparsfPathThisSubject);
            copyfile([niftiPathThisSubject niftiRsFmriFilenames{i} format], fmriDparsfPathThisSubject);
        end
    end
end

%% CHECK DATA FOR CONFIG
% Check slice order and number of slices for each manufacturer:
for s = 1 : numel(scanners)
    indicesDifferentPerManufacturer{s} = [];

    indicesManufacturer = find( strncmp(scanners{s}, scannerManufacturer,2) >0);
    spm_orderManufacturer{s} =  spm_order(indicesManufacturer);
    imageRealSizeManufacturer_voxels{s} = imageRealSize_voxels(indicesManufacturer,:);
    timePointsManufacturer{s} = timePoints(indicesManufacturer);
    sliceOrderManufacturer{s} =  sliceOrder(indicesManufacturer);
    % Check if all time points for this manufacturer are the same:
    if unique(timePointsManufacturer{s})
        timePointsPerManufacturer(s) = timePointsManufacturer{s}(1);
    else
        warning(sprintf('Not all of the fMRI scans have the same time points for %s manufacturer', scanners{s}))
    end
%     if unique(spm_orderManufacturer{s})
%         spmOrderPerManufacturer{s} = spm_orderManufacturer{s}{1};
%     else
%         warning(sprintf('Not all of the fMRI scans have the same time slice order for %s manufacturer', scanners{s}))
%     end
    j = 0;
    for i = 1 : numel(sliceOrderManufacturer{s})
        if(numel(sliceOrderManufacturer{s}{1}) == numel(sliceOrderManufacturer{s}{i}))
            if all(sliceOrderManufacturer{s}{1} == sliceOrderManufacturer{s}{i})
                j = j + 1;
            else
                indicesDifferentPerManufacturer{s} = [indicesDifferentPerManufacturer{s} indicesManufacturer(i)];
            end
        else
            indicesDifferentPerManufacturer{s} = [indicesDifferentPerManufacturer{s} indicesManufacturer(i)];
        end
    end
    if j == numel(sliceOrderManufacturer{s})
        sliceOrderPerManufacturer{s} = sliceOrderManufacturer{s}{1};
    else
        if numel(indicesDifferentPerManufacturer{s}) < j % If different are few
            sliceOrderPerManufacturer{s} = sliceOrderManufacturer{s}{1};
        else    % i == 1 is the different.
            sliceOrderPerManufacturer{s} = sliceOrderManufacturer{s}{indicesDifferentPerManufacturer{s}(1)};
        end
        warning(sprintf('%d of the fMRI scans have a different time slice order for %s manufacturer', numel(sliceOrderManufacturer{s})-j, scanners{s}))
    end
end
%% SAVE DATA
save([outputPath 'workspaceProcessingScript']);
%% RUN DPARSF
% Process each folder using the respective config file.
RunDparsf(dparsfConfigFilenamesPerManufacturer, [outputPath 'workspaceProcessingScript.mat'])

%% FUNCTION TO RUN DPARSF
function RunDparsf(dparsfConfigFilenamesPerManufacturer, workspaceFilename)
    if nargin == 2
        load(workspaceFilename)
    end
    for s = 1 : numel(scanners)
        configThisManufacturer = load(dparsfConfigFilenamesPerManufacturer{s});
        for g = 1 : numel(groups)
            dparsfThisGroupBasePath = [dparsfPerGroupPath groups{g} '_' scanners{s} '/'];
            dparsfThisGroupFunPath = [dparsfThisGroupBasePath fmriNameNifti '/'];
            listDir = dir(dparsfThisGroupFunPath);
            subjectNames = listDir(3:end); % Remove ., ..
            configThisManufacturer.Cfg.WorkingDir = dparsfThisGroupBasePath;
            configThisManufacturer.Cfg.DataProcessDir = dparsfThisGroupBasePath;
            configThisManufacturer.Cfg.SubjectID = {subjectNames(:).name};
            configThisManufacturer.Cfg.SubjectNum = numel(subjectNames);
            configThisManufacturer.Cfg.IsNeedReorientFunImgInteractively = 0; % Do not reorient manually.
            configThisManufacturer.Cfg.IsNeedReorientT1ImgInteractively = 0;
            %configThisManufacturer.Cfg.IsAllowGUI = 0;
            % File to the AAL atlas:
            configThisManufacturer.Cfg.CalFC.ROIDef = {[dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/Templates/aal.nii']};
            [Error]=DPARSF_run(configThisManufacturer.Cfg);
        end
    end
end
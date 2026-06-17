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
overwriteNifti = 0; % If 0 won't process existing Nfiti data in the DPARSF folder
format = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';
oasisPath = [adniPartitionPath '/OASIS/'];
dataPath = [oasisPath '/second_batch/OASIS_next_60CN_60AD_2026_05_26_RAW_DATA_120of120-001/OASIS_next_60CN_60AD_2026_05_26_RAW/'];

dparsfConfigFilenamesPerManufacturer = {['./DPARSF/configFiles/DPARSF_config_siemens.mat'], ...
    [dataPartitionPath './DPARSF/configFiles/DPARSF_config_ge.mat'],...
    [dataPartitionPath './DPARSF/configFiles/DPARSF_config_philips.mat']};

outputPath = [oasisPath '/second_batch/OASIS_next_60CN_60AD_2026_05_26_RAW_DATA_120of120-001/OASIS_next_60CN_60AD_2026_05_26_processed/'];
dparsfPath = [outputPath '/DPARSF/'];
niftiPath = [outputPath '/Nifti/'];
if ~isdir(dparsfPath)
    mkdir(dparsfPath)
end
if ~isdir(niftiPath)
    mkdir(niftiPath)
end
%% PARAMETERS THAT MUST BE THE SAME
repetitionTimeDefault_sec = 2.2;
timePointsDefault = 164;
%% NAME SEQUENCES
nameT1 = 't1_mprage_1x1x1';
namefMRIfolder = 'func';
% nameFieldmappingMag1 = 'gre_field_mapping_2mm_e1';
% nameFieldmappingMag2 = 'gre_field_mapping_2mm_e2';
% nameFieldmappingPhase = 'gre_field_mapping_2mm_phase';
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
%% CASES TO PROCESS
casesToProcess = [];
k = 1;
if isempty(casesToProcess)
    casesToProcess = {};
    casesToProcessAdvanced = {};
    pathRsFmri = {};
    pathRsFmriAdvanced = {};
    dirPath = dir(dataPath);
    for i = 3 : numel(dirPath)
        if(dirPath(i).isdir) 
            dirsFunc = dir([fullfile( dataPath, dirPath(i).name, namefMRIfolder), '*']);
            for j = 1 : numel(dirsFunc)
                % Get the nifti file:
                filenameNifti = dir(fullfile( dataPath, dirPath(i).name, dirsFunc(j).name, '*.nii.gz'));
                filenameJson = dir(fullfile( dataPath, dirPath(i).name, dirsFunc(j).name, '*.json'));
                subjectNames{k} = dirPath(i).name;
                filenamesfMRI{k} = fullfile( dataPath, dirPath(i).name, dirsFunc(j).name, filenameNifti.name);
                [path,namefMRIscan,ext] = fileparts(filenamesfMRI{k}); % double extension .nii.gz
                [path, namefMRIscan,ext] = fileparts(namefMRIscan);
                % Read the JSON file
                fid = fopen(fullfile( dataPath, dirPath(i).name, dirsFunc(j).name, filenameJson.name), 'r');
                raw = fread(fid, inf, 'uint8=>char')';
                fclose(fid);
                
                % Parse into a struct
                params{k} = jsondecode(raw);
                niftiParams = niftiinfo(filenamesfMRI{k});

                voxelSize_mm(k,:) = niftiParams.PixelDimensions(1:3);
                imageSize_voxels(k,:) = niftiParams.ImageSize(1:3); 
                tRs(k) = params{k}.RepetitionTime;
                scannerManufacturer{k} = params{k}.Manufacturer;
                scannerModel{k} = params{k}.ManufacturersModelName;
                numSlices(k) = imageSize_voxels(k,3);
                timePoints(k) = niftiParams.ImageSize(4);
                if strcmpi(scannerManufacturer{k}, 'Siemens')
                    % Philips does not include slice order, by defaults:
                    if rem(numSlices,2)
                        sliceOrder{k} = [2:2:numSlices 1:2:numSlices];
                    else
                        sliceOrder{k} = [1:2:numSlices 2:2:numSlices];
                    end
                else
                    error('Unexpected manufacturer.')
                end

                % Copy to DPARSF folder only if they match the required
                % parameters:
                if (tRs(k) == repetitionTimeDefault_sec) && (timePoints(k) == timePointsDefault)
                    fmriDparsfPathThisSubject = [fmriDparsfPath '/' namefMRIscan '/'];
                    mkdir(fmriDparsfPathThisSubject);
                    copyfile(filenamesfMRI{k}, fmriDparsfPathThisSubject);
                    copyfile(fullfile( dataPath, dirPath(i).name, dirsFunc(j).name, filenameJson.name), fmriDparsfPathThisSubject);        
                end
                k = k + 1;
            end
        end
    end
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

%% SAVE DATA
save([outputPath 'workspaceProcessingScript']);


%% CHECK DATA FOR CONFIG
% Check slice order and number of slices for each manufacturer:


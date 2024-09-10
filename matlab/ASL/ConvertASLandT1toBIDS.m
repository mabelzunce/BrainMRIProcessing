clear all
close all

%% NAME SEQUENCES AND VARIABLES
niftiExtension = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';

imagesDir = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Organized_DICOM/';
niftiDir = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/BIDS/rawdata/';

%% CREATE DATASET DESCRIPTION JSON 
datasetDescriptionJsonPath = [niftiDir 'dataset_description.json'];

% Fill with fields and values
datasetDescription = struct();
datasetDescription.Name = 'Long COVID Dataset';
datasetDescription.BIDSVersion = '1.9.0';

jsonString = jsonencode(datasetDescription);
prettyJsonString = prettyPrintJson(jsonString);

% Create and Write json
fileID = fopen(datasetDescriptionJsonPath, 'w');
fprintf(fileID, '%s', prettyJsonString);
fclose(fileID);


%% CASES TO PROCESS 
%subject = 'CP0061';
casesToProcess = []; 
if isempty(casesToProcess)
    dirPaths = dir(imagesDir);
    casesToProcess = {dirPaths(3:end).name};
end

%% ADD PATHS
addpath('/home/sol/BrainTools/dicm2nii/')   
addpath(genpath( '/home/sol/BrainTools/ExploreASL/'));


%% Convert DICOM to Nifti

%% Delete participants.tsv
filenameParticipantsTSV= [niftiDir 'participants.tsv'];
if exist(filenameParticipantsTSV, 'file') == 2
    % Delete the file
    delete(filenameParticipantsTSV);
end


% Define the source and destination directories
for i = 1 : numel(casesToProcess)
    subject = char(casesToProcess(i));
    subjectBIDS = ['sub-' subject]; %Subject Name using BIDS format 
    subjectBIDSDir = [niftiDir subjectBIDS '/']; %Path to subject in BIDS dir 
    subjectDirDICOM = [imagesDir subject '/DICOM/'];
    %subjectDirNifti = [niftiDir subject '/Nifti/'];
    
    converterOut = dicm2nii(subjectDirDICOM, niftiDir, 'bids'); 

    % Rename T2 with FLAIR
    subjectT2NameJson = [subjectBIDSDir 'anat/' subjectBIDS '_T2w.json'];
    subjectFLAIRNameJson = [subjectBIDSDir 'anat/' subjectBIDS '_FLAIR.json'];

    subjectT2NameImage = [subjectBIDSDir 'anat/' subjectBIDS '_T2w.nii.gz'];
    subjectFLAIRNameImage= [subjectBIDSDir 'anat/' subjectBIDS '_FLAIR.nii.gz'];
    
    movefile(subjectT2NameJson, subjectFLAIRNameJson);
    movefile(subjectT2NameImage, subjectFLAIRNameImage);

    % Modify Json
    subjectDICOMHeaders =  [subjectBIDSDir dcmHeadersFilename];

    % Read Tags 
    ASLTagName = ['sub_' subject '_asl' ];
    dcmTagsSubject = load(subjectDICOMHeaders);
    dcmTagsASL = getfield(dcmTagsSubject.h,ASLTagName);


    % Read and Modify Json
    subjectASLJson =  [subjectBIDSDir 'perf/' subjectBIDS '_asl.json' ];

    fid = fopen(subjectASLJson, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);

    ASLJsonData = jsondecode(str);
    %ASLJsonData.RepetitionTimePreparation = repetitionTime; 
    ASLJsonData.ArterialSpinLabelingType = 'PASL';
    
    %ASLJsonData.PostLabelingDelay	=
    ASLJsonData.BackgroundSuppression	= true; 
    ASLJsonData.M0Type = 'Absent'; 
    ASLJsonData.TotalAcquiredPairs	= 10; 
    ASLJsonData.BolusCutOffFlag = true; 
    ASLJsonData.BolusCutOffDelayTime = dcmTagsASL.CSAImageHeaderInfo.QCData(4)/1000; % cambiar esto
    ASLJsonData.PostLabelingDelay = dcmTagsASL.CSAImageHeaderInfo.QCData(5)/1000; % cambiar esto
    ASLJsonData.BolusCutOffTechnique  =  'QUIPSS-II';
    ASLJsonData.MagneticFieldStrength = dcmTagsASL.MagneticFieldStrength; 
    ASLJsonData.ProcedureStepDescription = dcmTagsASL.BodyPartExamined;
    
    ASLJsonData.BodyPartExamined = dcmTagsASL.BodyPartExamined;
    ASLJsonData.AcquisitionNumber = dcmTagsASL.AcquisitionNumber;
    ASLJsonData.StationName = dcmTagsASL.StationName; 
    ASLJsonData.PatientPosition = dcmTagsASL.PatientPosition;
    ASLJsonData.SoftwareVersions = dcmTagsASL.SoftwareVersion;
    ASLJsonData.ProtocolName = dcmTagsASL.ProtocolName;
    ASLJsonData.AcquisitionNumber = dcmTagsASL.AcquisitionNumber; 
    ASLJsonData.SpacingBetweenSlices = dcmTagsASL.SpacingBetweenSlices; 
    ASLJsonData.SAR = dcmTagsASL.SAR; 
    ASLJsonData.RepetitionTime = dcmTagsASL.RepetitionTime / 1000; 
    %ASLJsonData.BolusDuration = dcmTagsASL.CSAImageHeaderInfo.QCData(4)/1000;
    %ASLJsonData.InversionTime = dcmTagsASL.CSAImageHeaderInfo.QCData(5)/1000;
    ASLJsonData.DwellTime = dcmTagsASL.RealDwellTime * 1000;
    ASLJsonData.ProcedureStepDescription =  dcmTagsASL.PerformedProcedureStepDescription;
    ASLJsonData.ProcedureStepDescription =  dcmTagsASL.PerformedProcedureStepDescription;

    % Write JSON 
    modifiedJson = jsonencode(ASLJsonData);
    prettyJsonString = prettyPrintJson(modifiedJson);
    
    % Write the modified JSON back to the file
    fid = fopen(subjectASLJson, 'w');
    if fid == -1
        error('Cannot open file for writing.');
    end
    %fwrite(fid, modifiedJson, 'char');
    fwrite(fid, prettyJsonString, 'char');
    fclose(fid);

    % Create aslcontext_tsv 
    aslContextFile = [subjectBIDSDir 'perf/' subjectBIDS '_aslcontext.tsv'];
    
    volumeType = repmat({'label', 'control'}, 1, 10);
    volumeType = volumeType';

    fileID = fopen(aslContextFile, 'w');

    % Write header
    fprintf(fileID, 'volume_type\n');
    
    % Write the data
    for j = 1:length(volumeType)
        fprintf(fileID, '%s\n', volumeType{j});
    end
    
    % Close the file
    fclose(fileID);
end


%% MODIFIED PARTICIPANTS TSV

% Read the table

tbl = readtable(filenameParticipantsTSV, 'FileType', 'text', 'Delimiter', '\t');

% Check if the first element starts with 'sub-'
if ~ startsWith(tbl.participant_id{1}, 'sub-')
    % Identify IDs that do not start with 'sub-'
    mask = ~startsWith(tbl.participant_id, 'sub-');

    % Only modify and write back if there are changes
    if any(mask)
        tbl.participant_id(mask) = strcat('sub-', tbl.participant_id(mask));
        writetable(tbl, filenameParticipantsTSV, 'FileType', 'text', 'Delimiter', '\t');
    end
end
%% RUN ASL PIPELINE
%[x] = ExploreASL('/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/BIDS/', 0, 1);

%% PATH SUBJECT DICOM DIR 
%ConvertDicomFolderStructure_CarefulSlow(subjectDir)

%subjectDir = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/ASL_BIDS/Raw_DICOM/CP0060';
%% COPY TO NEW DIR 
% % Copy the directory
% copyfile(subjectDir, destDir);
% 
% % Verify that the copying process was successful
% if exist(destDir, 'dir')
%     disp('Directory copied successfully.');
% else
%     disp('There was an issue copying the directory.');
% end

clear all
close all

%% NAME SEQUENCES AND VARIABLES
niftiExtension = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';

imagesDir = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Organized_DICOM/';
niftiDir = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/BIDS/';

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
% Define the source and destination directories
for i = 1 : numel(casesToProcess)
    subject = char(casesToProcess(i));
    subjectBIDS = ['sub-' subject]; %Subject Name using BIDS format 
    subjectBIDSDir = [niftiDir subjectBIDS '/']; %Path to subject in BIDS dir 
    subjectDirDICOM = [imagesDir subject '/DICOM/'];
    %subjectDirNifti = [niftiDir subject '/Nifti/'];
    
    converterOut = dicm2nii(subjectDirDICOM, niftiDir, 'bids'); 
    
    subjectDICOMHeaders =  [subjectBIDSDir dcmHeadersFilename];

    % Read Tags 
    ASLTagName = ['sub_' subject '_asl' ];
    dcmTagsSubject = load(subjectDICOMHeaders);
    dcmTagsASL = getfield(dcmTagsSubject.h,ASLTagName);

    % Extract necessary fields 
    repetitionTime = dcmTagsASL.RepetitionTime/1000; 

    % Read and Modify Json
    subjectASLJson =  [subjectBIDSDir 'perf/' subjectBIDS '_asl.json' ];

    fid = fopen(subjectASLJson, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);

    ASLJsonData = jsondecode(str);
    ASLJsonData.RepetitionTimePreparation = repetitionTime; 
    ASLJsonData.ArterialSpinLabelingType = 'PASL';
    %ASLJsonData.PostLabelingDelay	=
    ASLJsonData.BackgroundSuppression	= true; 
    ASLJsonData.M0Type = 'Absent'; 
    ASLJsonData.TotalAcquiredPairs	= 10; 
    %ASLJsonData.MagneticFieldStrength =
    
    % Write JSON 
    modifiedJson = jsonencode(ASLJsonData);
    %prettyJsonString = prettyPrintJson(modifiedJson);
    
    % Write the modified JSON back to the file
    fid = fopen(subjectASLJson, 'w');
    if fid == -1
        error('Cannot open file for writing.');
    end
    fwrite(fid, modifiedJson, 'char');
    fclose(fid);


end


%% MODIFIED PARTICIPANTS TSV
filename = [niftiDir 'participants.tsv'];
tbl = readtable(filename , 'FileType', 'text', 'Delimiter', '\t');

mask = ~startsWith(tbl.participant_id, 'sub-');
tbl.participant_id(mask) = strcat('sub-', tbl.participant_id(mask));

writetable(tbl, filename, 'FileType', 'text', 'Delimiter', '\t');

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

%% ADD PATHS
addpath('D:\UNSAM\Brain\dicm2nii\')
addpath(genpath('D:\UNSAM\Brain\DPABI_V6.2_220915\'))
addpath('D:\UNSAM\Brain\spm12\spm12\')
addpath('D:\UNSAM\Brain\DPABI_V6.2_220915\DPARSF')
%% PATHS AND FILENAMES
format = '.nii.gz';
dataPath = 'D:/UNSAM/Brain/CovidProject/Estudio/MRIs/';
dparsfPath = [dataPath '/DPARSF/'];
if ~isdir(dparsfPath)
    mkdir(dparsfPath)
end
%% NAME SEQUENCES
nameT1 = 't1_mprage_1x1x1';
nameFmri = 'funcional_3s_30te_2_4_iso';
nameFieldmappingMag1 = 'gre_field_mapping_2mm_e1';
nameFieldmappingMag2 = 'gre_field_mapping_2mm_e2';
nameFieldmappingPhase = 'gre_field_mapping_2mm_phase';
%% CASES TO PROCESS
casesToProcess = {'CP0002', 'CP0006', 'CP0007', 'CP0008', 'CP0009', 'CP0010'};
%% FOLDERS FOR FREEFURSFER
freesurferSubjectsPath =  'D:/freesurfer_subjects/';
filenameFreesurferScript = [dataPath 'freesurferScript.txt'];
if ~isdir(freesurferSubjectsPath)
    mkdir(freesurferSubjectsPath);
end
freesurferLines = [];
%% FOLDERS FOR DPARSF
dataPath = 'D:/UNSAM/Brain/CovidProject/Estudio/MRIs/';
t1NameNifti = 'T1Img';
fmriNameNifti = 'FunImg';
fieldmapBaseNameNifti = 'FieldMap';
fieldmapPhaseNameNifti = 'PhaseDiffImg';
fieldmapMagNameNifti = 'Maginute1Img';

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

fieldmapMagDparsfPath = [dparsfPath '/' fieldmapBaseNameNifti '/' fieldmapMagNameNifti '/'];
if ~isdir(fieldmapMagDparsfPath)
    mkdir(fieldmapMagDparsfPath);
end

%% PROCESS EACH CASE
for i = 1 : numel(casesToProcess)
    %% CONVERT DICOM TO NIFTI
    pathThisSubject = [dataPath casesToProcess{i}];
    dicomPath = [pathThisSubject '/DICOM/'];
    niftiPath = [pathThisSubject '/Nifti/'];
    dcmTags = dicm2nii(dicomPath, niftiPath, format);
    %% CREATE FOLDERS FOR MELODIC
    melodicPath = [pathThisSubject '/Nifti/'];
    %% CREATE FOLDERS FOR THIS SUBJECT DPARSF AND COPY FILES
    t1DparsfPathThisSubject = [t1DparsfPath '/' casesToProcess{i} '/'];
    fmriDparsfPathThisSubject = [fmriDparsfPath '/' casesToProcess{i} '/'];
    fieldmapMagDparsfPathThisSubject = [fieldmapMagDparsfPath '/' casesToProcess{i} '/'];
    fieldmapPhaseDparsfPathThisSubject = [fieldmapPhaseDparsfPath '/' casesToProcess{i} '/'];
    mkdir(t1DparsfPathThisSubject);
    copyfile([niftiPath nameT1 format], t1DparsfPathThisSubject);
    mkdir(fmriDparsfPathThisSubject);
    copyfile([niftiPath nameFmri format], fmriDparsfPathThisSubject);
    mkdir(fieldmapMagDparsfPathThisSubject);
    copyfile([niftiPath nameFieldmappingMag1 format], fieldmapMagDparsfPathThisSubject);
    mkdir(fieldmapPhaseDparsfPathThisSubject);
    copyfile([niftiPath nameFieldmappingPhase format], fieldmapPhaseDparsfPathThisSubject);
    %% CREATE FREESURFER SCRIPT WITH ALL THE LINES
    freesurferLines = [freesurferLines sprintf('recon-all -subject %s -i %s\n', casesToProcess{i}, [niftiPath nameT1 format])];
end
%% FIELDMAPPING CORRECTION
% h.gre_field_mapping_2mm_e1.deltaTE
fprintf('Parameters for fieldmap correction. TE1:%.1fms, TE2:%.1fms, dTE:%.1fms\n', h.gre_field_mapping_2mm_e1.EchoTime, ...
    h.gre_field_mapping_2mm_e2.EchoTime, h.gre_field_mapping_2mm_e1.deltaTE);
%% RUN FREESURFER
freesurferLines = [freesurferLines 'recon-all -subject '];
for i = 1 : numel(casesToProcess)
    freesurferLines = [freesurferLines 'recon-all -subject ' casesToProcess{i} ' -all \n'];
end
freesurferLines = [freesurferLines ' -all\n'];
% If in windows, save the script and run it later in Ubuntu for windows:
fid = fopen(filenameFreesurferScript,'w');
fprintf(fid, freesurferLines);
fclose(fid);
%% RUN DPARSF




% getDemographicsData
clear all
close all

%% ADD PATHS
addpath('D:\UNSAM\Brain\dicm2nii\')
addpath(genpath('D:\UNSAM\Brain\DPABI_V6.2_220915\'))
addpath('D:\UNSAM\Brain\spm12\spm12\')
addpath('D:\UNSAM\Brain\DPABI_V6.2_220915\DPARSF')

%% PATHS AND FILENAMES
overwriteNifti = 1; % If 0 won't process existing Nfiti data in the DPARSF folder
format = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';
adniPath = 'F:/ADNIdata/';
dataPath = [adniPath '/ADNIdicom/'];

processedPath = 'K:/procesados/';

%% ADNI DATABASE COLLECTION
filenameADNICollection = 'ADNI3_fmri_12_19_2022.csv';
filenameADNICollection = 'ADNI_fMRI_TR=3000ms_2_16_2023';
nameRsFmriInCollectionDataBase = {'Axial rsfMRI (Eyes Open)', 'Axial fcMRI', 'Axial - Advanced fMRI'}; %{'Axial rsfMRI (Eyes Open)', 'Axial fcMRI'};
nameRsFmriAdvancedInCollectionDataBase = {'Axial MB rsfMRI'};%'Axial MB rsfMRI (Eyes Open)';
% Load adni dicom collection csv:
adniCollectionData = readtable([adniPath filenameADNICollection]);

%% NAME SEQUENCES
nameT1 = 't1_mprage_1x1x1';
nameRsFmri = {'Axial_rsfMRI', 'Axial_fcMRI'};
nameRsFmriAdvanced = 'Axial_MB_rsfMRI';
% nameFieldmappingMag1 = 'gre_field_mapping_2mm_e1';
% nameFieldmappingMag2 = 'gre_field_mapping_2mm_e2';
% nameFieldmappingPhase = 'gre_field_mapping_2mm_phase';

%% GET DEMOGRAPHICS
% Goes through the processed data, ordered in subfolders with batches and
% then by DPARSF subfolder
dirPath = dir(processedPath);
% Only directories:
dirPath = dirPath([dirPath(:).isdir]);
% Go through each subfolder
casesToProcess = {};
indexInDatabase = 1;
indexOutOfDatabase = 1;
for i = 3 : numel(dirPath)
    batchPath = [processedPath dirPath(i).name '/'];
    funImgPath = [batchPath 'FunImg/'];
    subjectPaths = dir(funImgPath);
    for j = 3 : numel(subjectPaths)
        indicesSubject = find(strcmp(subjectPaths(j).name, adniCollectionData.Subject) > 0);
        if ~isempty(indicesSubject)
            [sortedDates, indicesDates] = sort(adniCollectionData.AcqDate(indicesSubject));
            % Use the last date:
            indexToUse = indicesSubject(indicesDates(end));
            subjectName{indexInDatabase} = adniCollectionData.Subject{indexToUse};
            sex{indexInDatabase} = adniCollectionData.Sex{indexToUse};
            age_years(indexInDatabase) = adniCollectionData.Age(indexToUse);
            group{indexInDatabase} = adniCollectionData.Group{indexToUse};
            visit{indexInDatabase} = adniCollectionData.Visit{indexToUse};
            date(indexInDatabase) = adniCollectionData.AcqDate(indexToUse);
            indexInDatabase = indexInDatabase + 1;
        else
            casesNotFound{indexOutOfDatabase} = subjectPaths(j).name;
            indexOutOfDatabase = indexOutOfDatabase + 1;
            warning(sprintf('Subject %s was not found in the collection database.', subjectPaths(j).name));
        end
    end
end

%%
figure;
set(gcf, 'Position', [100 100 800 600]);
% Alternative plot:
ha = boxchart(age_years,'GroupByColor',group);
legend
ylabel('Age [years]')
saveas(gca, [processedPath 'demographics_age'], 'epsc');
saveas(gca, [processedPath 'demographics_age'], 'tif');

%% Check each sex
indicesMaleAD = find((strcmp(sex, 'M')>0 & strcmp(group, 'AD')) > 0);
indicesMaleCN = find((strcmp(sex, 'M')>0 & strcmp(group, 'CN')) > 0);
indicesFemaleAD = find((strcmp(sex, 'F')>0 & strcmp(group, 'AD')) > 0);
indicesFemaleCN = find((strcmp(sex, 'F')>0 & strcmp(group, 'CN')) > 0);
numMalesAD = numel(indicesMaleAD);
numMalesCN = numel(indicesMaleCN);
numFemalesAD = numel(indicesFemaleAD);
numFemalesCN = numel(indicesFemaleCN);
disp(sprintf('AD group. Males: %d. Females: %d.\n', numMalesAD, numFemalesAD))
disp(sprintf('CN group. Males: %d. Females: %d.\n', numMalesCN, numFemalesCN))
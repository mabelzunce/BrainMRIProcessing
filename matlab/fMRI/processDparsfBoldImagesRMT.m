clear all
%close all
dataPartitionPath = '/data/'; %'D:/'
adniPartitionPath = '/data_imaging/'; %'F:/'
%% ADD PATHS
dpabiPath = [dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/'];
addpath([dataPartitionPath 'UNSAM/Brain/dicm2nii/'])
addpath(genpath(dpabiPath))
addpath([dataPartitionPath 'UNSAM/Brain/spm12/spm12/'])
addpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/DPARSF/'])
%% DATA PATHS
adniPath = [adniPartitionPath '/ADNIdata/ADNI3_Advanced_fMRI/ADNI3_AdvancedFmri/'];
processedDataPath = [adniPath '/Processed/DPARSF/'];
dparsfPreprocessedFolder = 'FunImgRWSDCF';
dparsfPreprocessedDataPath = fullfile(processedDataPath, dparsfPreprocessedFolder);
outputPath = [adniPath '/RMT_Analysis/'];
%% FOLDERS FOR DPARSF
fmriNameNifti = 'FunImg';
fmriFullyPreprocessedNameNiftiFolder = 'FunImgRWSDCF'; %Filtered_4DVolume.nii
fmriFullyPreprocessedBoldSignalsFolder = '/Results/ROISignals_FunImgARWSDCF/';
fmriFullyPreprocessedBoldSignalsPath = fullfile(processedDataPath, dparsfPreprocessedFolder);
clear all
%close all
dataPartitionPath = '/data/'; %'D:/'
imagingPartitionPath = '/data_imaging/'; %'F:/'
%% ADD PATHS
dpabiPath = [dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/'];
addpath([dataPartitionPath 'UNSAM/Brain/dicm2nii/'])
addpath(genpath(dpabiPath))
addpath([dataPartitionPath 'UNSAM/Brain/spm12/spm12/'])
addpath([dataPartitionPath 'UNSAM/Brain/DPABI_V6.2_220915/DPARSF/'])
%% CONFIG
createPlots = 0;
%% DATA PATHS
dataPath = [dataPartitionPath '/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/DataBaseDrive-20231031T195909Z-001/DataBase/ADNI_fMRI_screening/AAL/'];
dataPath = [imagingPartitionPath '/CovidProject/Estudio/PreprocessedMRI/DPARSF/'];
dataPath = [imagingPartitionPath '/ADNIdata/Trini/'];
numBatches = 5;
dcmHeadersFilename = 'dcmHeaders.mat';
%% GO THROUGH BATCHES AND SUBFOLDERS
k = 1;
for i = 1 : numBatches
    batchSubdir = fullfile(dataPath, sprintf('Batch%d', i), '2mm', 'FunImg');
    % Get list of subjects:
    fmriSubjects = dir(batchSubdir);
    fmriSubjects = fmriSubjects([fmriSubjects(:).isdir]);
    fmriSubjects = fmriSubjects(3:end);
    numel(fmriSubjects)
    for j = 1 : numel(fmriSubjects)
        niftiFilename = dir(fullfile(batchSubdir, fmriSubjects(j).name, '*.nii'));
        niftiFilename = fullfile(batchSubdir, fmriSubjects(j).name, niftiFilename.name);
        if exist(niftiFilename)
            subjectsNames{k} = fmriSubjects(j).name;
            niftiFilenames{k} = niftiFilename;
            dcmTags{k} = load(fullfile(batchSubdir, subjectsNames{k}, dcmHeadersFilename));
            sequencesPerSubject{k} = fieldnames(dcmTags{k}.h);
            dcmTagsRsFmri{k} = getfield(dcmTags{k}.h,sequencesPerSubject{k}{1});
            fMRI_tR(k) = dcmTagsRsFmri{k}.RepetitionTime;
            fMRI_tE(k) = dcmTagsRsFmri{k}.EchoTime;
            if isfield(dcmTagsRsFmri{k}, 'MosaicRefAcqTimes')
                fMRI_sliceAcqTimes{k} = dcmTagsRsFmri{k}.MosaicRefAcqTimes;
                fMRI_fslSliceOrcer{k} = dcmTagsRsFmri{k}.SliceTiming;
                [times, fMRI_sliceOrder{k}] = sort(dcmTagsRsFmri{k}.MosaicRefAcqTimes);
            elseif isfield(dcmTagsRsFmri{k}, 'SliceTiming')
                fMRI_sliceAcqTimes{k} = (0.5 - dcmTagsRsFmri{k}.SliceTiming) * dcmTagsRsFmri{k}.RepetitionTime;
                fMRI_fslSliceOrcer{k} = dcmTagsRsFmri{k}.SliceTiming;
                [times, fMRI_sliceOrder{k}] = sort(fMRI_sliceAcqTimes{k});
            end
            %image = niftiread([niftifMriFilenames{k}]);
            info = niftiinfo(niftiFilename);
            fMRI_imageSize_voxels(k,:) = info.ImageSize;
            fMRI_voxelSize_mm(k,:) = info.PixelDimensions;
            fMRI_inPlanePhaseEncodingDirection{k} = dcmTagsRsFmri{k}.InPlanePhaseEncodingDirection;
            fMRI_unwarpDirection(k,:) = dcmTagsRsFmri{k}.UnwarpDirection;
            if isfield(dcmTagsRsFmri{k}, 'EffectiveEPIEchoSpacing')
                fMRI_effectiveEPIEchoSpacing(k) = dcmTagsRsFmri{k}.EffectiveEPIEchoSpacing;
            end
            manufacturer{k} = dcmTagsRsFmri{k}.Manufacturer;
            model{k} = dcmTagsRsFmri{k}.ManufacturerModelName;
            k = k + 1;
        else
            disp(niftiFilename)
        end
    end
end
%% CHECK SLICE TIMINGSS
indexPhilips = find(strncmpi(manufacturer, 'Philips', numel('Philips'))>0);
indexSiemens = find(strncmpi(manufacturer, 'siemens', numel('siemens'))>0);
indexGE = find(strncmpi(manufacturer, 'GE', numel('GE'))>0);
if (numel(indexGE) + numel(indexSiemens) + numel(indexPhilips)) ~= numel(manufacturer)
    error('Not all scanners classified')
end
%/home/martin/data_imaging/ADNIdata/Trini/Batch1/2mm/FunImg/006_S_0498/
%% REWRITE FILES FOR REPROCESSING
outputPath = [imagingPartitionPath '/ADNIdata/DPARSF_initialvisit/Nifti/'];
if ~exist(outputPath)
    mkdir(outputPath);
end
for i = 1 : numel(subjectsNames)
    subjectOutputDir = fullfile(outputPath, subjectsNames{i});
    image = niftiread(niftiFilenames{i});
    info = niftiinfo(niftiFilenames{i});
    outputFilename = fullfile(outputPath, subjectsNames{i}, subjectsNames{i});
    mkdir(subjectOutputDir);
    niftiwrite(image, outputFilename, info,  'Compressed', 1);
    [path, file, ext] = fileparts(niftiFilenames{i});
    copyfile(fullfile(path, dcmHeadersFilename), fullfile(subjectOutputDir, dcmHeadersFilename));
end
% Copy dicom converter file:
for i = 1 : numel(subjectsNames)
    subjectOutputDir = fullfile(outputPath, subjectsNames{i});
    [path, file, ext] = fileparts(niftiFilenames{i});
    copyfile(fullfile(path, dcmHeadersFilename), fullfile(subjectOutputDir, dcmHeadersFilename));
end
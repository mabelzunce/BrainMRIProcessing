% Test slicing for each manufcaturer.
clear all
%% PHILIPS
dcmPath = 'D:\ADNIdata\ADNIdicom\006_S_6277\Axial_fcMRI__EYES_OPEN_\2019-04-25_11_39_17.0\I1158788\';
dcmFiles = dir([dcmPath '*.dcm']);
disp(sprintf('Number of files: %d', numel(dcmFiles)))
for i = 1 : numel(dcmFiles)
    disp(i)
    dcmHeader{i} = dicominfo([dcmPath dcmFiles(i).name]);
    triggerTimes(i) = dcmHeader{i}.TriggerTime;
    tempId(i) = dcmHeader{i}.TemporalPositionIdentifier;
    sliceLocation(i) = dcmHeader{i}.SliceLocation;
    sliceId(i) = dcmHeader{i}.Private_2001_100a;
    image(:,:, sliceId(i),tempId(i)) = dicomread([dcmPath dcmFiles(i).name]);
end
%% Slice order for philips
TRsec = 3;
nSlices = 48;
TA = TRsec/nSlices; %assumes no temporal gap between volumes
bidsSliceTiming=[0:TA:TRsec-TA]; %ascending
if false %descending
    bidsSliceTiming = flip(bidsSliceTiming);
end
if true %interleaved
    order = [1:2:nSlices 2:2:nSlices]
bidsSliceTiming(order) = bidsSliceTiming;
end
%report results
fprintf('SliceTiming: [\n');
for i = 1 : nSlices
    fprintf(' %g', bidsSliceTiming(i));
    if (i < nSlices)
        fprintf(',\n');
    else
        fprintf(' ],\n');
    end
end
%% SIEMENS
clear all
dcmPath = 'D:\ADNIdata\ADNIdicom\002_S_0413\Axial_rsfMRI__Eyes_Open_\2017-06-21_13_23_38.0\I863058\';
dcmFiles = dir([dcmPath '*.dcm']);
disp(sprintf('Number of files: %d', numel(dcmFiles)));
for i = 1 : numel(dcmFiles)
    %disp(i)
    dcmHeader{i} = dicominfo([dcmPath dcmFiles(i).name]);
    timeId(i) = dcmHeader{i}.InstanceNumber;
    timeId2(i) = dcmHeader{i}.AcquisitionNumber;
    sliceLocation(i) = dcmHeader{i}.SliceLocation;
    typecast(uint8(dcmHeader{end}.Private_0019_1028), 'double')
    typecast(uint8(dcmHeader{end}.Unknown_0040_0000), 'uint32')
    %sliceId(i) = dcmHeader{i}.Private_2001_100a;
    image(:,:,:,timeId(i)) = reshape(dicomread([dcmPath dcmFiles(i).name]), dcmHeader{i}.AcquisitionMatrix(1), dcmHeader{i}.AcquisitionMatrix(4), []);
end
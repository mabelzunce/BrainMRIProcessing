clear all
close all
%% GET DEMOGRAPHICS FOR PAPER
dataPath = '/home/martin/data_imaging/ADNIdata/fMRI_ADNI2_ADNI3_initial_visit/';
filenameSummaryTable = fullfile(dataPath, 'SubjctsDataAndTestsSchaefer2018_100Parcels_17Networks.csv');
%% READ DATA
summary = readtable(filenameSummaryTable);
% Remove low quality or data with preprocessing errors:
subjectsToRemove = {'114_S_6039', '035_S_6953', '128_S_2002', '031_S_4021', ...
    '130_S_5231', '130_S_6647'};
% Remove data to exclude:
indicesToExclude = [];
for i = 1 : numel(subjectsToRemove)
    ind=find(strcmp(subjectsToRemove{i}, summary.SubjectID)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind(1)];
    end
end
summary(indicesToExclude,:) = [];
%% GET STAS
scannerTypes = {'Philips', 'Siemens', 'GE'};
[numPerStudy, channelsPhase] = hist(categorical(summary.Phase));
% Get stats for the scanners:
for i = 1 : numel(summary.ImagingProtocol)
    splits = strsplit(summary.ImagingProtocol{i}, ';');
    manufacturerLines{i} = splits(contains(splits, 'Manufacturer'));
    for j = 1 : 3
        if contains(manufacturerLines{i}, scannerTypes{j},'IgnoreCase',true)
            scanner{i} = scannerTypes{j};
            break;
        end
    end
    if numel(scanner)<i % Didn;t match with any, proably GE because different names used
        scanner{i} = 'GE';
    end
end
[numPerScanner, channelsScanner] = hist(categorical(scanner));
%%
groups = {'AD', 'CN', 'MCI'}
for i = 1 : numel(groups)
    indicesThisGroup = contains(summary.ResearchGroup, groups{i});
    subjectsPerGroup(i) = sum(indicesThisGroup);
    sexByGroup{i} = hist(categorical(summary(indicesThisGroup,:).Sex), {'M','F'});
    meanAge(i) = mean(summary(indicesThisGroup,:).Age);
    stdAge(i) = std(summary(indicesThisGroup,:).Age);
end

% summary.Sex
% summary.Phase
% summary.ResearchGroup
% summary.ImagingProtocol


clear all
pathData = '/home/martin/data/UNSAM/CovidProject/Estudio/DataAnalysis/BrainMorphometry/';
studyDataPath = '/home/martin/data/UNSAM/CovidProject/Estudio/';
resultsPath = '../DataAnalysis/';

subjectsToExclude = {'CP0011', 'CP0015','CP0035', 'CP0144', 'CP0106'};

summaryExcelFilename = fullfile(studyDataPath, 'RespuestasCuestionarioEvaluacionCognitivaResonancia.xlsx');

leftParcDktAncovaFilename = 'left_parcellation_dkt_ancova.csv';
parcelationDktFilename = 'parcellation_DKT.csv';
segmentationFilename = 'segmentation.csv';
brainVolumes = 'brain_volumes.csv';
leftParcDktAncova = readtable([pathData, leftParcDktAncovaFilename]);
%% SUMMARY
summaryTable = readtable(summaryExcelFilename);
% Remove data to exclude:
indicesToExclude = [];
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, summaryTable.ID)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind'];
    end
end
summaryTable(indicesToExclude,:) = [];
% Remove individuals without MRI
indicesNoMri = ~strcmp(summaryTable.SeHizoResonancia, 'Si');
summaryTable(indicesNoMri,:) = [];
%% VOLUMES
indicesToExclude = [];
volumes = readtable([pathData, brainVolumes]);
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, volumes.subject)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind'];
    end
end
volumes(indicesToExclude,:) = [];
% Sort by name
volumes = sortrows(volumes,1);
%% SEGMENTATION
indicesToExclude = [];
segmentation = readtable([pathData, segmentationFilename]);
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, segmentation.subject)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind'];
    end
end
segmentation(indicesToExclude,:) = [];
% Sort by name
segmenation = sortrows(segmentation,1);
% Check if the names match:


indicesRoisToUse = [4:9 13 14 16 22:30];
matRoisToUse = segmentation(:,indicesRoisToUse);
groups = table2cell(segmentation(:,end-1));
matSegmentationsPerSubject = table2array(matRoisToUse);

for i = 1 : size(matSegmentationsPerSubject, 2)
    [p(i), table{i}, stats{i}] = anova1(matSegmentationsPerSubject(:,i), groups, 'off');
end
%% FDR
alpha = 0.5;
[p_ordered, indices] = sort(p);
n = numel(p);
i = [1 : n];
critical = i.*alpha./n;
% Get the max p value for p < critical:
p_corr = max(p_ordered(p_ordered < critical));
indices_sign = p_ordered<p_corr;
disp(p_corr)
% p signficativos
indices(indices_sign)
p(indices(indices_sign))
%% HOLM
p_cric = alpha./(n+1-i);
indice_p_greater = find(p_ordered > p_cric);
p_thresh = p_ordered(indice_p_greater(1))
indices_sign2 = p_ordered < p_thresh;
% p signficativos
indices(indices_sign2)
p(indices(indices_sign2))
%% PERMUTATION
indicesCOVID = strcmp(groups, 'COVID Prolongado');
indicesControl = ~indicesCOVID;
diffMeanGroups = mean(matSegmentationsPerSubject(indicesCOVID,:)) - mean(matSegmentationsPerSubject(indicesControl,:));
numPermutations = 5000;
for j = 1 : numPermutations
    reshuffled_groups = groups(randperm(numel(groups)));
    indicesCOVIDShuffled = strcmp(reshuffled_groups, 'COVID Prolongado');
    indicesControlShuffled = ~indicesCOVIDShuffled;
    for i = 1 : size(matSegmentationsPerSubject, 2)
        [p_perm(j,i), table{i}, stats{i}] = anova1(matSegmentationsPerSubject(:,i), reshuffled_groups, 'off');
    end
    diffMeanGroupsPerm(j,:) = mean(matSegmentationsPerSubject(indicesCOVIDShuffled,:)) - mean(matSegmentationsPerSubject(indicesControlShuffled,:));
end
%% compute p values for the measured difference:

p_perm = sum(diffMeanGroups > diffMeanGroupsPerm,1)./numPermutations;
% Z scores
z_scores = (diffMeanGroups - mean(diffMeanGroupsPerm,1))./std(diffMeanGroupsPerm, [], 1)
% p values norm dist
p_perm_norm = normcdf(diffMeanGroups, mean(diffMeanGroupsPerm,1), std(diffMeanGroupsPerm, [], 1))
%% COMPARE ALL p
p
p_perm
% Significant:
indices_sign_perm = find(p_perm < 0.05)
indices(indices_sign)
indices(indices_sign2)

%% NOW FOR PARCELLATIONS
parcellation = readtable([pathData, parcelationDktFilename]);

indicesLeft = strcmp(parcellation.Hemisphere, 'lh');
indicesRight = strcmp(parcellation.Hemisphere, 'rh');
parcellationLeft = parcellation(indicesLeft,:);
parcellationRight = parcellation(indicesRight,:);

indicesRoisToUse = [2:35];
matRoisToUseLeft = parcellationLeft(:,indicesRoisToUse);
matParcellationLeftPerSubject = table2array(matRoisToUseLeft);
matRoisToUseRight = parcellationRight(:,indicesRoisToUse);
matParcellationRightPerSubject = table2array(matRoisToUseRight);
parcellationBoth = [matParcellationLeftPerSubject + matParcellationRightPerSubject]./2;

groups = table2cell(parcellationRight(:,end-1));
indicesCOVID = strcmp(groups, 'COVID Prolongado');
indicesControl = ~indicesCOVID;



for i = 1 : size(parcellationBoth, 2)
    [p(i), table{i}, stats{i}] = anova1(parcellationBoth(:,i), groups, 'off');
end
% FDR
alpha = 0.5;
[p_ordered, indices] = sort(p);
n = numel(p);
i = [1 : n];
critical = i.*alpha./n;
% Get the max p value for p < critical:
if ~isempty(p_ordered(p_ordered < critical))
    p_corr = max(p_ordered(p_ordered < critical));
    indices_sign = p_ordered<p_corr;
    disp(p_corr)
    % p signficativos
    if ~isempty(indices_sign)
        indices(indices_sign)
        p(indices(indices_sign))
    end
end
% HOLM
p_cric = alpha./(n+1-i);
indice_p_greater = find(p_ordered > p_cric);
if ~isempty(indice_p_greater)
    p_thresh = p_ordered(indice_p_greater(1))
    indices_sign2 = p_ordered < p_thresh;
    % p signficativos
    if ~isempty(indices_sign2)
        indices(indices_sign2)
        p(indices(indices_sign2))
    end
end
% PERMUTATION
diffMeanGroups = mean(parcellationBoth(indicesCOVID,:)) - mean(parcellationBoth(indicesControl,:));
numPermutations = 5000;
for j = 1 : numPermutations
    reshuffled_groups = groups(randperm(numel(groups)));
    indicesCOVIDShuffled = strcmp(reshuffled_groups, 'COVID Prolongado');
    indicesControlShuffled = ~indicesCOVIDShuffled;
    for i = 1 : size(parcellationBoth, 2)
        [p_perm(j,i), table{i}, stats{i}] = anova1(parcellationBoth(:,i), reshuffled_groups, 'off');
    end
    diffMeanGroupsParcPerm(j,:) = mean(parcellationBoth(indicesCOVIDShuffled,:)) - mean(parcellationBoth(indicesControlShuffled,:));
end
% compute p values for the measured difference:

p_perm = sum(diffMeanGroups > diffMeanGroupsParcPerm,1)./numPermutations;
% Z scores
z_scores = (diffMeanGroups - mean(diffMeanGroupsParcPerm,1))./std(diffMeanGroupsParcPerm, [], 1)
% p values norm dist
p_perm_norm = normcdf(diffMeanGroups, mean(diffMeanGroupsParcPerm,1), std(diffMeanGroupsParcPerm, [], 1))

% COMPARE ALL p
p
p_perm
% Significant:
indices_sign_perm = find(p_perm < 0.05)
indices(indices_sign)
indices(indices_sign2)
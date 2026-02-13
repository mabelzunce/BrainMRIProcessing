clear all
pathData = '/home/martin/data/UNSAM/CovidProject/Estudio/DataAnalysis/BrainMorphometry/';
studyDataPath = '/home/martin/data/UNSAM/CovidProject/Estudio/';
resultsPath = [studyDataPath '/DataAnalysis/'];

subjectsToExclude = {'CP0011', 'CP0015','CP0035', 'CP0144', 'CP0106'};

summaryExcelFilename = fullfile(studyDataPath, 'RespuestasCuestionarioEvaluacionCognitivaResonancia.xlsx');

parcelationVolFilename = 'parcellation_volumes_clean.csv';
parcelationThicknessFilename = 'parcellation_thickness_clean.csv';
segmentationFilename = 'segmentation_clean.csv';
brainVolumes = 'brain_volumes.csv';
%leftParcDktAncova = readtable([pathData, leftParcDktAncovaFilename]);
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
if sum(~strcmp(segmenation.subject, summaryTable.ID)) > 0
    error('Tables dont match');
end
if sum(~strcmp(segmenation.Grupo, summaryTable.Grupo)) > 0
    error('Tables dont match');
end
if sum(~strcmp(segmenation.Edad, summaryTable.Edad)) > 0
    warning(sprintf('Tables dont match for age of % subject', sum(~strcmp(segmenation.Edad, summaryTable.Edad))));
end
age = segmenation.Edad;
eTIV = segmenation.Estimated_Total_Intracranial_Volume;
indicesRoisToUse = [5 : size(segmentation, 2)];
matRoisToUse = segmentation(:,indicesRoisToUse);
group = categorical(segmentation.Grupo);
matSegmentationsPerSubject = table2array(matRoisToUse);
tbl2 = table(group, age, eTIV);
tbl = table(group, eTIV);
for i = 1 : size(matSegmentationsPerSubject, 2)
    metric = matSegmentationsPerSubject(:,i);    
    %[p(i), table{i}, stats{i}] = anova(tbl, "metric ~ group + age + eTIV");
    anovaTable{i} = anova(tbl, metric,categoricalFactors=["group"]);
    anovaTable2{i} = anova(tbl2, metric,categoricalFactors=["group"]);
    p(i) = anovaTable{i}.stats.pValue(1);
    p_age(i) = anovaTable{i}.stats.pValue(2);
    p_eTIV(i) = anovaTable{i}.stats.pValue(3);
end
% The same with fitlm
for i = 1 : size(matSegmentationsPerSubject, 2)
    metric = matSegmentationsPerSubject(:,i); 
    tblfitlm = table(metric, group, age, eTIV);
    % Model with all covariates
    mdl_all{i} = fitlm(tblfitlm, 'metric ~ 1 + group + age + eTIV');
    mdl_eTIV{i} = fitlm(tblfitlm, 'metric ~ 1 + group + eTIV');
end
% Get regressed out data:
age_demean = age - mean(age);
eTIV_demean = eTIV - mean(eTIV);
for i = 1 : size(matSegmentationsPerSubject, 2)
    metric = matSegmentationsPerSubject(:,i);
    tbl_reg = table(metric, group, age_demean, eTIV_demean);
    mdl = fitlm(tbl_reg,'metric ~ age_demean + eTIV_demean'); %age was not significant mdl = fitlm(tbl_reg,'metric ~ age_demean + eTIV_demean'); 
    coefficients = mdl.Coefficients.Estimate;
    % regress out variables:
    metric_reg_age_etiv(:,i) = metric - coefficients(2)*age_demean - coefficients(3)*eTIV_demean;
end

% FDR
alpha = 0.05;
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
    indices(indices_sign)
    p(indices(indices_sign))
else
    p_corr = [];
    indices_sign = [];
end
% HOLM
p_cric = alpha./(n+1-i);
indice_p_greater = find(p_ordered > p_cric);
p_thresh = p_ordered(indice_p_greater(1))
indices_sign2 = p_ordered < p_thresh;
% p signficativos
indices(indices_sign2)
p(indices(indices_sign2))

% PERMUTATION
group = segmentation.Grupo;
indicesCOVID = strcmp(group, 'COVID');
indicesControl = ~indicesCOVID;
diffMeanGroups = mean(matSegmentationsPerSubject(indicesCOVID,:)) - mean(matSegmentationsPerSubject(indicesControl,:));
diffMeanGroupsReg = mean(metric_reg_age_etiv(indicesCOVID,:)) - mean(metric_reg_age_etiv(indicesControl,:));
% original t_statistic:
[h_t_orig, p_t_orig, ci_t_orig,stats_t_orig] = ttest2(matSegmentationsPerSubject(indicesCOVID,:),matSegmentationsPerSubject(indicesControl,:));
tstats_orig = stats_t_orig.tstat;
[h_t_orig, p_t_reg_orig, ci_t_orig,stats_t_reg_orig] = ttest2(metric_reg_age_etiv(indicesCOVID,:),metric_reg_age_etiv(indicesControl,:));
tstatsReg_orig = stats_t_reg_orig.tstat;
numPermutations = 5000;
for j = 1 : numPermutations
    % Shuffle
    reshuffled_groups = group(randperm(numel(group)));
    indicesCOVIDShuffled = strcmp(reshuffled_groups, 'COVID');
    indicesControlShuffled = ~indicesCOVIDShuffled;
    % Another option, completely random
    indicesRandomGroups = randsample([1 -1], size(matSegmentationsPerSubject,1), true);
    indicesCOVIDShuffled = indicesRandomGroups == -1;
    indicesControlShuffled = indicesRandomGroups == 1;
    reshuffled_groups(indicesRandomGroups == -1) = {'COVID'};
    reshuffled_groups(indicesRandomGroups == 1) = {'CONTROL'};

    tbl_shuffled = table(reshuffled_groups, age, eTIV);
    covariates = [age, eTIV];
    % for i = 1 : size(matSegmentationsPerSubject, 2)
    %     anovaTable{i} = anova(tbl_shuffled, metric,categoricalFactors=["reshuffled_groups"]);
    %     mdl = fitlm(covariates,metric);
    %     p_anova_perm(i) = anovaTable{i}.stats.pValue(1);
    %     p_age_perm(i) = anovaTable{i}.stats.pValue(2);
    %     p_eTIV_perm(i) = anovaTable{i}.stats.pValue(3);
    %     %[p_perm(j,i), table{i}, stats{i}] = anova1(matSegmentationsPerSubject(:,i), reshuffled_groups, 'off');
    % end
    diffMeanGroupsPerm(j,:) = mean(matSegmentationsPerSubject(indicesCOVIDShuffled,:)) - mean(matSegmentationsPerSubject(indicesControlShuffled,:));
    [h_aux,p_aux,ci,stats] = ttest2(matSegmentationsPerSubject(indicesCOVIDShuffled,:),matSegmentationsPerSubject(indicesControlShuffled,:));
    tstats(j,:) = stats.tstat;
    pstats(j,:) = p_aux;

    diffMeanGroupsPermReg(j,:) = mean(metric_reg_age_etiv(indicesCOVIDShuffled,:)) - mean(metric_reg_age_etiv(indicesControlShuffled,:));
    [h_aux,p_aux,ci_aux,stats] = ttest2(metric_reg_age_etiv(indicesCOVIDShuffled,:),metric_reg_age_etiv(indicesControlShuffled,:));
    tstatsReg(j,:) = stats.tstat;
    pstatsReg(j,:) = p_aux;
end
tstats_max_perm = max(tstats,[],2); % max tstats for each permutation
tstats_min_perm = min(tstats,[],2); % min tstats for each permutation
tstats_max_perm_reg = max(tstatsReg,[],2); % max tstats for each permutation
tstats_min_perm_reg = min(tstatsReg,[],2); % min tstats for each permutation

% Get number of values lower than the tstats_max_perm for each of the
% original values:
num_greater = sum(tstats_max_perm>tstats_orig,1);
num_lower = sum(tstats_min_perm<tstats_orig,1);
num_greater_reg = sum(tstats_max_perm_reg >= tstats_orig,1);
num_lower_reg = sum(tstats_min_perm_reg <= tstats_orig,1);
p_greater = num_greater./ numPermutations;
p_lower = num_lower ./ numPermutations;
p_greater_reg = num_greater_reg./ numPermutations;
p_lower_reg = num_lower_reg ./ numPermutations;

% Westfall and young permutation
[p_orig_sorted, indices_p_orig_sorted] = sort(p_t_orig);
[p_reg_orig_sorted, indices_p_reg_orig_sorted] = sort(p_t_reg_orig);
pstats_sorted = sort(pstats,2);
pstatsReg_sorted = sort(pstatsReg,2);
tstats_sorted = sort(tstats,2);
tstatsReg_sorted = sort(tstatsReg,2);
[tstats_orig_sorted, indices_tstats_orig_sorted] = sort(tstats_orig);
[tstats_reg_orig_sorted, indices_tstats_reg_orig_sorted] = sort(tstatsReg_orig);
% For each ROI compare the ith highes permuted statistics:
for i = 1 : size(matSegmentationsPerSubject,2)
    perm_min = tstats_sorted(:, i);
    permReg_min = tstatsReg_sorted(:, i);
    pPerm_min = pstats_sorted(:, i);
    pRegPerm_min = pstatsReg_sorted(:, i);
    %adjusted_p_values(indices_tstats_orig_sorted(i)) = mean(perm_min <= tstats_orig_sorted(i));
    adjusted_p_values_sorted(i) = mean(perm_min <= tstats_orig_sorted(i));
    adjusted_p_values2_sorted(i) = mean(pPerm_min <= p_orig_sorted(i));
    adjusted_p_values_reg_sorted(i) = mean(permReg_min <= tstats_orig_sorted(i));
    adjusted_p_values2_reg_sorted(i) = mean(pRegPerm_min <= p_reg_orig_sorted(i));
end
csvwrite(fullfile(resultsPath, 'segmentation_pValues.csv'), p_t_orig);
csvwrite(fullfile(resultsPath, 'segmentationReg_pValues.csv'), p_t_reg_orig);
csvwrite(fullfile(resultsPath, 'segmentationPermutations_pValues.csv'), pstats);
csvwrite(fullfile(resultsPath, 'segmentationPermutationReg_pValues.csv'), pstatsReg);

csvwrite(fullfile(resultsPath, 'segmentation_tstats.csv'), tstats_orig);
csvwrite(fullfile(resultsPath, 'segmentationReg_tstats.csv'), tstatsReg_orig);
csvwrite(fullfile(resultsPath, 'segmentationPermutations_tstats.csv'), tstats);
csvwrite(fullfile(resultsPath, 'segmentationPermutationReg_tstats.csv'), tstatsReg);


% for i = 2 : size(matSegmentationsPerSubject,2)
%     if adjusted_p_values_sorted(i) < adjusted_p_values_sorted(i-1)
%         adjusted_p_values_sorted(i) = adjusted_p_values_sorted(i - 1);
%     end
%     if adjusted_p_values_reg_sorted(i) < adjusted_p_values_reg_sorted(i-1)
%         adjusted_p_values_reg_sorted(i) = adjusted_p_values_reg_sorted(i - 1);
%     end
%     if adjusted_p_values2_sorted(i) < adjusted_p_values2_sorted(i-1)
%         adjusted_p_values2_sorted(i) = adjusted_p_values2_sorted(i - 1);
%     end
%     if adjusted_p_values2_reg_sorted(i) < adjusted_p_values2_reg_sorted(i-1)
%         adjusted_p_values2_reg_sorted(i) = adjusted_p_values2_reg_sorted(i - 1);
%     end
% end
index_p_perm_critical = find(adjusted_p_values_sorted>0.05);
p_perm_critical = adjusted_p_values_sorted(index_p_perm_critical(1));
index_p_perm_critical2 = find(adjusted_p_values2_sorted>0.05);
p_perm_critical2 = adjusted_p_values2_sorted(index_p_perm_critical2(1));
index_p_reg_perm_critical = find(adjusted_p_values_reg_sorted>0.05);
p_reg_perm_critical = adjusted_p_values_reg_sorted(index_p_reg_perm_critical(1));
index_p_reg_perm_critical2 = find(adjusted_p_values2_reg_sorted>0.05);
p_reg_perm_critical2 = adjusted_p_values2_reg_sorted(index_p_reg_perm_critical2(1));

adjusted_p_values = adjusted_p_values_sorted(indices_tstats_orig_sorted);
adjusted_p_values2 = adjusted_p_values2_sorted(indices_p_orig_sorted);
adjusted_p_reg_values = adjusted_p_values_reg_sorted(indices_tstats_reg_orig_sorted);
adjusted_p_reg_values2 = adjusted_p_values2_reg_sorted(indices_tstats_reg_orig_sorted);

indices_sign_perm_WY = find(adjusted_p_values < p_perm_critical);
indices_sign_perm_2_WY = find(adjusted_p_values2 < p_perm_critical2);
indices_sign_perm_reg_WY = find(adjusted_p_reg_values < p_reg_perm_critical);
indices_sign_perm_reg_2_WY = find(adjusted_p_values2_reg_sorted < p_reg_perm_critical2);

segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_WY)
segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_2_WY)
segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg_WY)
segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg_2_WY)

% compute p values for the measured difference:
p_perm = sum(diffMeanGroups > diffMeanGroupsPerm,1)./numPermutations;
p_perm_reg = sum(diffMeanGroupsReg > diffMeanGroupsPermReg,1)./numPermutations;
% Z scores
z_scores = (diffMeanGroups - mean(diffMeanGroupsPerm,1))./std(diffMeanGroupsPerm, [], 1);
z_scores_reg = (diffMeanGroupsReg - mean(diffMeanGroupsPermReg,1))./std(diffMeanGroupsPermReg, [], 1);
% p values norm dist
p_perm_norm = normcdf(diffMeanGroups, mean(diffMeanGroupsPerm,1), std(diffMeanGroupsPerm, [], 1));
p_perm_norm_reg = normcdf(diffMeanGroupsReg, mean(diffMeanGroupsPermReg,1), std(diffMeanGroupsPermReg, [], 1));

% Significant:
disp('Regions significantly different without multiple comparisons:')
p
indices_sign_no_mc = find(p<0.05)
segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_no_mc)
disp('Regions significantly different with permutation:')
p_perm
indices_sign_perm = find(p_perm < 0.05)
indices_sign_perm_reg = find(p_perm_reg < 0.05)
segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm)
disp('Regions significantly different with FDR:')
indices(indices_sign)
segmentation(:,indicesRoisToUse).Properties.VariableNames(indices(indices_sign))
disp('Regions significantly different with Holm:')
indices(indices_sign2)
segmentation(:,indicesRoisToUse).Properties.VariableNames(indices(indices_sign2))

% Save results
resultsSegmentation.no_mc.p_values = p;
resultsSegmentation.no_mc.p_values_ord = p_ordered;
resultsSegmentation.no_mc.indices_sign = find(p<0.05);
resultsSegmentation.no_mc.p_sign = p(find(p<0.05));
resultsSegmentation.no_mc.names = segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_no_mc);
% fdr
resultsSegmentation.fdr.p_th= p_corr;
resultsSegmentation.fdr.indices_sign = indices(indices_sign);
resultsSegmentation.fdr.p_sign = p(indices(indices_sign));
resultsSegmentation.fdr.names = segmentation(:,indicesRoisToUse).Properties.VariableNames(indices(indices_sign));

% holm
resultsSegmentation.holm.p_th = p_thresh;
resultsSegmentation.holm.indices_sign = indices(indices_sign2);
resultsSegmentation.holm.p_sign = p(indices(indices_sign2));
resultsSegmentation.holm.names = segmentation(:,indicesRoisToUse).Properties.VariableNames(indices(indices_sign2));

% permutation
resultsSegmentation.permutation.p_values = p_perm;
resultsSegmentation.permutation.adjusted_p_values = adjusted_p_values;
resultsSegmentation.permutation.diff_means_perm = diffMeanGroupsPerm;
resultsSegmentation.permutation.indices_sign = indices_sign_perm;
resultsSegmentation.permutation.p_sign = p_perm(indices_sign_perm);
resultsSegmentation.permutation.z_scores = z_scores;
resultsSegmentation.permutation.p_norm = p_perm_norm;
resultsSegmentation.permutation.names = segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm);

% permutation
resultsSegmentation.permutationWestfallAndYoung.p_values = adjusted_p_values;
resultsSegmentation.permutationWestfallAndYoung.indices_sign = indices_sign_perm_WY;
resultsSegmentation.permutationWestfallAndYoung.p_sign = adjusted_p_values(indices_sign_perm_WY);
resultsSegmentation.permutationWestfallAndYoung.tstats = tstats;
resultsSegmentation.permutationWestfallAndYoung.names = segmentation(:,indicesRoisToUse).Properties.VariableNames(adjusted_p_values<0.05);

% permutation
resultsSegmentation.permutationWestfallAndYoung_reg.p_values = adjusted_p_reg_values;
resultsSegmentation.permutationWestfallAndYoung_reg.indices_sign =indices_sign_perm_reg_WY;
resultsSegmentation.permutationWestfallAndYoung_reg.p_sign = adjusted_p_reg_values(indices_sign_perm_reg_WY);
resultsSegmentation.permutationWestfallAndYoung_reg.tstats = tstatsReg;
resultsSegmentation.permutationWestfallAndYoung_reg.names = segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg_WY);


% permutation with variables regressed out
resultsSegmentation.permutation_reg.p_values = p_perm_reg;
resultsSegmentation.permutation_reg.adjusted_p_values = adjusted_p_reg_values;
resultsSegmentation.permutation_reg.diff_means_perm = diffMeanGroupsPermReg;
resultsSegmentation.permutation_reg.indices_sign = indices_sign_perm_reg;
resultsSegmentation.permutation_reg.p_sign = p_perm_reg(indices_sign_perm_reg);
resultsSegmentation.permutation_reg.z_scores = z_scores_reg;
resultsSegmentation.permutation_reg.p_norm = p_perm_norm_reg;
resultsSegmentation.permutation_reg.names = segmentation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg);

save([resultsPath 'segmentation_multiple_comparisons.mat'], 'resultsSegmentation')

% compare p values and regions:
[p' p_perm' adjusted_p_reg_values'] % p includes age an etvi as covariates

[p(indices_sign_perm)' p_perm(indices_sign_perm)' adjusted_p_values(indices_sign_perm)']
resultsSegmentation.permutation.names
%% NOW FOR PARCELLATIONS
clear p p_ordered anovaTable p_perm p_age_perm p_eTIV_perm diffMeanGroupsPerm diffMeanGroupsPermReg tstats tstatsReg pstats pstatsReg

parcellation = readtable([pathData, parcelationVolFilename]);
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, parcellation.subject)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind'];
    end
end
parcellation(indicesToExclude,:) = [];
% Sort by name
parcellation = sortrows(parcellation,1);

% Check if the names match:
if sum(~strcmp(parcellation.subject, summaryTable.ID)) > 0
    error('Tables dont match');
end
if sum(~strcmp(parcellation.Grupo, summaryTable.Grupo)) > 0
    error('Tables dont match');
end
if sum(~strcmp(parcellation.Edad, summaryTable.Edad)) > 0
    warning(sprintf('Tables dont match for age of %d subject', sum(~strcmp(segmenation.Edad, summaryTable.Edad))));
end
age = parcellation.Edad;
eTIV = parcellation.Estimated_Total_Intracranial_Volume;
indicesRoisToUse = [2 : size(parcellation, 2)-3];
matRoisToUse = parcellation(:,indicesRoisToUse);
group = categorical(parcellation.Grupo);
matSegmentationsPerSubject = table2array(matRoisToUse);
tbl = table(group, age, eTIV);
tbl2 = table(group, eTIV);
for i = 1 : size(matSegmentationsPerSubject, 2)
    metric = matSegmentationsPerSubject(:,i);    
    %[p(i), table{i}, stats{i}] = anova(tbl, "metric ~ group + age + eTIV");
    anovaTable{i} = anova(tbl, metric,categoricalFactors=["group"]);
    anovaTable2{i} = anova(tbl2, metric,categoricalFactors=["group"]);
    p_anova(i) = anovaTable{i}.stats.pValue(1);
    p_age(i) = anovaTable{i}.stats.pValue(2);
    p_eTIV(i) = anovaTable{i}.stats.pValue(3);
    p2(i) = anovaTable2{i}.stats.pValue(1);
    p_kw(i) = kruskalwallis(metric, group, 'off');
    [h(i), p(i)] = ttest2(matSegmentationsPerSubject(strcmp("CONTROL", string(group)),i), matSegmentationsPerSubject(strcmp("COVID", string(group)),i));
end

% Get regressed out data:
age_demean = age - mean(age);
eTIV_demean = eTIV - mean(eTIV);
for i = 1 : size(matSegmentationsPerSubject, 2)
    metric = matSegmentationsPerSubject(:,i);
    tbl_reg = table(metric, group, age_demean, eTIV_demean);
    mdl = fitlm(tbl_reg,'metric ~ age_demean + eTIV_demean');
    coefficients = mdl.Coefficients.Estimate;
    % regress out variables:
    metric_reg_age_etiv(:,i) = metric - coefficients(3)*eTIV_demean ;%- coefficients(2)*age_demean;
end

% FDR
alpha = 0.05;
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
diffMeanGroups = mean(matSegmentationsPerSubject(indicesCOVID,:)) - mean(matSegmentationsPerSubject(indicesControl,:));
diffMeanGroupsReg = mean(metric_reg_age_etiv(indicesCOVID,:)) - mean(metric_reg_age_etiv(indicesControl,:));
% original t_statistic:
[h_t_orig, p_t_orig, ci_t_orig,stats_t_orig] = ttest2(matSegmentationsPerSubject(indicesCOVID,:),matSegmentationsPerSubject(indicesControl,:));
tstats_orig = stats_t_orig.tstat;
[h_t_orig, p_t_reg_orig, ci_t_orig,stats_t_reg_orig] = ttest2(metric_reg_age_etiv(indicesCOVID,:),metric_reg_age_etiv(indicesControl,:));
tstatsReg_orig = stats_t_reg_orig.tstat;
numPermutations = 10000;
for j = 1 : numPermutations
    reshuffled_groups = group(randperm(numel(group)));
    indicesCOVIDShuffled = strcmp(string(reshuffled_groups), 'COVID');
    indicesControlShuffled = ~indicesCOVIDShuffled;
    tbl_shuffled = table(reshuffled_groups, age, eTIV);
    covariates = [age, eTIV];
    % Another option, completely random
    indicesRandomGroups = randsample([1 -1], size(matSegmentationsPerSubject,1), true);
    indicesCOVIDShuffled = indicesRandomGroups == -1;
    indicesControlShuffled = indicesRandomGroups == 1;
    reshuffled_groups(indicesRandomGroups == -1) = {'COVID'};
    reshuffled_groups(indicesRandomGroups == 1) = {'CONTROL'};
    % for i = 1 : size(matSegmentationsPerSubject, 2)
    %     anovaTable{i} = anova(tbl_shuffled, metric,categoricalFactors=["reshuffled_groups"]);
    %     mdl = fitlm(covariates,metric);
    %     p_anova_perm(i) = anovaTable{i}.stats.pValue(1);
    %     p_age_perm(i) = anovaTable{i}.stats.pValue(2);
    %     p_eTIV_perm(i) = anovaTable{i}.stats.pValue(3);
    %     %[p_perm(j,i), table{i}, stats{i}] = anova1(matSegmentationsPerSubject(:,i), reshuffled_groups, 'off');
    % end
    diffMeanGroupsPerm(j,:) = mean(matSegmentationsPerSubject(indicesCOVIDShuffled,:)) - mean(matSegmentationsPerSubject(indicesControlShuffled,:));
    [h_aux,p_aux,ci,stats] = ttest2(matSegmentationsPerSubject(indicesCOVIDShuffled,:),matSegmentationsPerSubject(indicesControlShuffled,:));
    tstats(j,:) = stats.tstat;
    pstats(j,:) = p_aux;
    diffMeanGroupsPermReg(j,:) = mean(metric_reg_age_etiv(indicesCOVIDShuffled,:)) - mean(metric_reg_age_etiv(indicesControlShuffled,:));
    [h_aux,p_aux,ci_aux,stats] = ttest2(metric_reg_age_etiv(indicesCOVIDShuffled,:),metric_reg_age_etiv(indicesControlShuffled,:));
    tstatsReg(j,:) = stats.tstat;
    pstatsReg(j,:) = p_aux;
end
tstats_max_perm = max(tstats,[],2); % max tstats for each permutation
tstats_min_perm = min(tstats,[],2); % min tstats for each permutation
tstats_max_perm_reg = max(tstatsReg,[],2); % max tstats for each permutation
tstats_min_perm_reg = min(tstatsReg,[],2); % min tstats for each permutation

% Get number of values lower than the tstats_max_perm for each of the
% original values:
num_greater = sum(tstats_max_perm>tstats_orig,1);
num_lower = sum(tstats_min_perm<tstats_orig,1);
num_greater_reg = sum(tstats_max_perm_reg >= tstats_orig,1);
num_lower_reg = sum(tstats_min_perm_reg <= tstats_orig,1);
p_greater = num_greater./ numPermutations;
p_lower = num_lower ./ numPermutations;
p_greater_reg = num_greater_reg./ numPermutations;
p_lower_reg = num_lower_reg ./ numPermutations;

% Westfall and young permutation
[p_orig_sorted, indices_p_orig_sorted] = sort(p_t_orig);
[p_reg_orig_sorted, indices_p_reg_orig_sorted] = sort(p_t_reg_orig);
pstats_sorted = sort(pstats,2);
pstatsReg_sorted = sort(pstatsReg,2);
tstats_sorted = sort(tstats,2);
tstatsReg_sorted = sort(tstatsReg,2);
[tstats_orig_sorted, indices_tstats_orig_sorted] = sort(tstats_orig);
[tstats_reg_orig_sorted, indices_tstats_reg_orig_sorted] = sort(tstatsReg_orig);
% For each ROI compare the ith highes permuted statistics:
for i = 1 : size(matSegmentationsPerSubject,2)
    perm_min = tstats_sorted(:, i);
    permReg_min = tstatsReg_sorted(:, i);
    pPerm_min = pstats_sorted(:, i);
    pRegPerm_min = pstatsReg_sorted(:, i);
    %adjusted_p_values(indices_tstats_orig_sorted(i)) = mean(perm_min <= tstats_orig_sorted(i));
    adjusted_p_values_sorted(i) = mean(perm_min <= tstats_orig_sorted(i));
    adjusted_p_values2_sorted(i) = mean(pPerm_min <= p_orig_sorted(i));
    adjusted_p_values_reg_sorted(i) = mean(permReg_min <= tstats_orig_sorted(i));
    adjusted_p_values2_reg_sorted(i) = mean(pRegPerm_min <= p_reg_orig_sorted(i));
end
csvwrite(fullfile(resultsPath, 'corticalThikness_pValues.csv'), p_t_orig);
csvwrite(fullfile(resultsPath, 'corticalThiknessReg_pValues.csv'), p_t_reg_orig);
csvwrite(fullfile(resultsPath, 'corticalThiknessPermutations_pValues.csv'), pstats);
csvwrite(fullfile(resultsPath, 'corticalThiknessPermutationReg_pValues.csv'), pstatsReg);

csvwrite(fullfile(resultsPath, 'corticalThikness_tstats.csv'), tstats_orig);
csvwrite(fullfile(resultsPath, 'corticalThiknessReg_tstats.csv'), tstatsReg_orig);
csvwrite(fullfile(resultsPath, 'corticalThiknessPermutations_tstats.csv'), tstats);
csvwrite(fullfile(resultsPath, 'corticalThiknessPermutationReg_tstats.csv'), tstatsReg);


% for i = 2 : size(matSegmentationsPerSubject,2)
%     if adjusted_p_values_sorted(i) < adjusted_p_values_sorted(i-1)
%         adjusted_p_values_sorted(i) = adjusted_p_values_sorted(i - 1);
%     end
%     if adjusted_p_values_reg_sorted(i) < adjusted_p_values_reg_sorted(i-1)
%         adjusted_p_values_reg_sorted(i) = adjusted_p_values_reg_sorted(i - 1);
%     end
%     if adjusted_p_values2_sorted(i) < adjusted_p_values2_sorted(i-1)
%         adjusted_p_values2_sorted(i) = adjusted_p_values2_sorted(i - 1);
%     end
%     if adjusted_p_values2_reg_sorted(i) < adjusted_p_values2_reg_sorted(i-1)
%         adjusted_p_values2_reg_sorted(i) = adjusted_p_values2_reg_sorted(i - 1);
%     end
% end
index_p_perm_critical = find(adjusted_p_values_sorted>0.05);
p_perm_critical = adjusted_p_values_sorted(index_p_perm_critical(1));
index_p_perm_critical2 = find(adjusted_p_values2_sorted>0.05);
p_perm_critical2 = adjusted_p_values2_sorted(index_p_perm_critical2(1));
index_p_reg_perm_critical = find(adjusted_p_values_reg_sorted>0.05);
p_reg_perm_critical = adjusted_p_values_reg_sorted(index_p_reg_perm_critical(1));
index_p_reg_perm_critical2 = find(adjusted_p_values2_reg_sorted>0.05);
p_reg_perm_critical2 = adjusted_p_values2_reg_sorted(index_p_reg_perm_critical2(1));

adjusted_p_values = adjusted_p_values_sorted(indices_tstats_orig_sorted);
adjusted_p_values2 = adjusted_p_values2_sorted(indices_p_orig_sorted);
adjusted_p_reg_values = adjusted_p_values_reg_sorted(indices_tstats_reg_orig_sorted);
adjusted_p_reg_values2 = adjusted_p_values2_reg_sorted(indices_tstats_reg_orig_sorted);

indices_sign_perm_WY = find(adjusted_p_values < p_perm_critical);
indices_sign_perm_2_WY = find(adjusted_p_values2 < p_perm_critical2);
indices_sign_perm_reg_WY = find(adjusted_p_reg_values < p_reg_perm_critical);
indices_sign_perm_reg_2_WY = find(adjusted_p_values2_reg_sorted < p_reg_perm_critical2);

parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_WY)
parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_2_WY)
parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg_WY)
parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg_2_WY)


% compute p values for the measured difference:
p_perm = sum(diffMeanGroups > diffMeanGroupsPerm,1)./numPermutations;
p_perm_reg = sum(diffMeanGroupsReg > diffMeanGroupsPermReg,1)./numPermutations;
% Z scores
z_scores = (diffMeanGroups - mean(diffMeanGroupsPerm,1))./std(diffMeanGroupsPerm, [], 1);
z_scores_reg = (diffMeanGroupsReg - mean(diffMeanGroupsPermReg,1))./std(diffMeanGroupsPermReg, [], 1);
% p values norm dist
p_perm_norm = normcdf(diffMeanGroups, mean(diffMeanGroupsPerm,1), std(diffMeanGroupsPerm, [], 1));
p_perm_norm_reg = normcdf(diffMeanGroupsReg, mean(diffMeanGroupsPermReg,1), std(diffMeanGroupsPermReg, [], 1));


% COMPARE ALL p
p
p_perm
% Significant:
indices_sign_perm = find(p_perm < 0.05)
indices_sign_perm_reg = find(p_perm_reg < 0.05)
indices(indices_sign)
indices(indices_sign2)

% Save results
resultsParcellationVol.no_mc.p_values = p;
resultsParcellationVol.no_mc.p_values_ord = p_ordered;
resultsParcellationVol.no_mc.indices_sign = find(p<0.05);
resultsParcellationVol.no_mc.p_sign = p(find(p<0.05));
resultsParcellationVol.no_mc.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(find(p<0.05));
% fdr
resultsParcellationVol.fdr.p_th= p_corr;
resultsParcellationVol.fdr.indices_sign = indices(indices_sign);
resultsParcellationVol.fdr.p_sign = p(indices(indices_sign));
resultsParcellationVol.fdr.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(indices(indices_sign));

% holm
resultsParcellationVol.holm.p_th = p_thresh;
resultsParcellationVol.holm.indices_sign = indices(indices_sign2);
resultsParcellationVol.holm.p_sign = p(indices(indices_sign2));
resultsParcellationVol.holm.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(indices(indices_sign2));

% permutation
resultsParcellationVol.permutation.p_values = p_perm;
resultsParcellationVol.permutation.diff_means_perm = diffMeanGroupsPerm;
resultsParcellationVol.permutation.indices_sign = indices_sign_perm;
resultsParcellationVol.permutation.p_sign = p_perm(indices_sign_perm);
resultsParcellationVol.permutation.z_scores = z_scores;
resultsParcellationVol.permutation.p_norm = p_perm_norm;
resultsParcellationVol.permutation.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm);

% permutation
resultsParcellationVol.permutationWestfallAndYoung.p_values = adjusted_p_values;
resultsParcellationVol.permutationWestfallAndYoung.indices_sign = find(adjusted_p_values<0.05);
resultsParcellationVol.permutationWestfallAndYoung.p_sign = p_perm(adjusted_p_values<0.05);
resultsParcellationVol.permutationWestfallAndYoung.tstats = tstats;
resultsParcellationVol.permutationWestfallAndYoung.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(adjusted_p_values<0.05);

% permutation with variables regressed out
resultsParcellationVol.permutation_reg.p_values = p_perm_reg;
resultsParcellationVol.permutation_reg.adjusted_p_values = adjusted_p_reg_values;
resultsParcellationVol.permutation_reg.diff_means_perm = diffMeanGroupsPermReg;
resultsParcellationVol.permutation_reg.indices_sign = indices_sign_perm_reg;
resultsParcellationVol.permutation_reg.p_sign = p_perm_reg(indices_sign_perm_reg);
resultsParcellationVol.permutation_reg.z_scores = z_scores_reg;
resultsParcellationVol.permutation_reg.p_norm = p_perm_norm_reg;
resultsParcellationVol.permutation_reg.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg);

% permutation
resultsParcellationVol.permutationWestfallAndYoung_reg.p_values = adjusted_p_reg_values;
resultsParcellationVol.permutationWestfallAndYoung_reg.indices_sign =indices_sign_perm_reg_WY;
resultsParcellationVol.permutationWestfallAndYoung_reg.p_sign = adjusted_p_reg_values(indices_sign_perm_reg_WY);
resultsParcellationVol.permutationWestfallAndYoung_reg.tstats = tstatsReg;
resultsParcellationVol.permutationWestfallAndYoung_reg.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg_WY);

save([resultsPath 'parcellation_vol_multiple_comparisons.mat'], 'resultsParcellationVol')

% compare p values and regions:
[p' p_perm' adjusted_p_values' adjusted_p_reg_values']

[p(p<0.05)' p_perm(p<0.05)' adjusted_p_values(p<0.05)' adjusted_p_reg_values(p<0.05)']
resultsParcellationVol.no_mc.names'
%% NOW FOR PARCELLATIONS thick
clear p p_ordered anovaTable p_perm p_age_perm p_eTIV_perm diffMeanGroupsPerm diffMeanGroupsPermReg tstats tstatsReg pstats pstatsReg
parcellation = readtable([pathData, parcelationThicknessFilename]);

% Check if the names match:
if sum(~strcmp(parcellation.subject, summaryTable.ID)) > 0
    error('Tables dont match');
end
if sum(~strcmp(parcellation.Grupo, summaryTable.Grupo)) > 0
    error('Tables dont match');
end
if sum(~strcmp(parcellation.Edad, summaryTable.Edad)) > 0
    warning(sprintf('Tables dont match for age of %d subject', sum(~strcmp(segmenation.Edad, summaryTable.Edad))));
end
age = parcellation.Edad;
indicesRoisToUse = [2 : size(parcellation, 2)-2];
matRoisToUse = parcellation(:,indicesRoisToUse);
groupStr = parcellation.Grupo;
group = categorical(parcellation.Grupo);
matSegmentationsPerSubject = table2array(matRoisToUse);
tbl = table(group, age); 
for i = 1 : size(matSegmentationsPerSubject, 2)
    metric = matSegmentationsPerSubject(:,i);    
    %[p(i), table{i}, stats{i}] = anova(tbl, "metric ~ group + age + eTIV");
    anovaTable{i} = anova(tbl, metric,categoricalFactors=["group"]);
    p_an(i) = anovaTable{i}.stats.pValue(1);
    p_age(i) = anovaTable{i}.stats.pValue(2);
    p_eTIV(i) = anovaTable{i}.stats.pValue(3);
    [h(i), p(i)] = ttest2(matSegmentationsPerSubject(strcmp("CONTROL", string(group)),i), matSegmentationsPerSubject(strcmp("COVID", string(group)),i));
end

% Get regressed out data:
age_demean = age - mean(age);
for i = 1 : size(matSegmentationsPerSubject, 2)
    metric = matSegmentationsPerSubject(:,i);
    tbl_reg = table(metric, group, age_demean);
    mdl = fitlm(tbl_reg,'metric ~ age_demean');
    coefficients = mdl.Coefficients.Estimate;
    % regress out variables:
    metric_reg_age(:,i) = metric - coefficients(2)*age_demean;
end


% FDR
alpha = 0.05;
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
diffMeanGroups = mean(matSegmentationsPerSubject(indicesCOVID,:)) - mean(matSegmentationsPerSubject(indicesControl,:));
diffMeanGroupsReg = mean(metric_reg_age(indicesCOVID,:)) - mean(metric_reg_age(indicesControl,:));
% original t_statistic:
[h_t_orig, p_t_orig, ci_t_orig,stats_t_orig] = ttest2(matSegmentationsPerSubject(indicesCOVID,:),matSegmentationsPerSubject(indicesControl,:));
tstats_orig = stats_t_orig.tstat;
[h_t_orig, p_t_reg_orig, ci_t_orig,stats_t_reg_orig] = ttest2(metric_reg_age(indicesCOVID,:),metric_reg_age(indicesControl,:));
tstatsReg_orig = stats_t_reg_orig.tstat;
numPermutations = 50000;
for j = 1 : numPermutations
    reshuffled_groups = group(randperm(numel(group)));
    indicesCOVIDShuffled = strcmp(string(reshuffled_groups), 'COVID');
    indicesControlShuffled = ~indicesCOVIDShuffled;
    % % Another option, completely random
    % indicesRandomGroups = randsample([1 -1], size(matSegmentationsPerSubject,1), true);
    % indicesCOVIDShuffled = indicesRandomGroups == -1;
    % indicesControlShuffled = indicesRandomGroups == 1;
    % reshuffled_groups(indicesRandomGroups == -1) = {'COVID'};
    % reshuffled_groups(indicesRandomGroups == 1) = {'CONTROL'};

    tbl_shuffled = table(reshuffled_groups, age, eTIV);
    covariates = [age, eTIV];
    % for i = 1 : size(matSegmentationsPerSubject, 2)
    %     anovaTable{i} = anova(tbl_shuffled, metric,categoricalFactors=["reshuffled_groups"]);
    %     mdl = fitlm(covariates,metric);
    %     p_anova_perm(i) = anovaTable{i}.stats.pValue(1);
    %     p_age_perm(i) = anovaTable{i}.stats.pValue(2);
    %     p_eTIV_perm(i) = anovaTable{i}.stats.pValue(3);
    %     %[p_perm(j,i), table{i}, stats{i}] = anova1(matSegmentationsPerSubject(:,i), reshuffled_groups, 'off');
    % end
    diffMeanGroupsPerm(j,:) = mean(matSegmentationsPerSubject(indicesCOVIDShuffled,:)) - mean(matSegmentationsPerSubject(indicesControlShuffled,:));
    [h_aux,p_aux,ci,stats] = ttest2(matSegmentationsPerSubject(indicesCOVIDShuffled,:),matSegmentationsPerSubject(indicesControlShuffled,:));
    tstats(j,:) = stats.tstat;
    pstats(j,:) = p_aux;
    diffMeanGroupsPermReg(j,:) = mean(metric_reg_age(indicesCOVIDShuffled,:)) - mean(metric_reg_age(indicesControlShuffled,:));
    [h_aux,p_aux,ci_aux,stats] = ttest2(metric_reg_age(indicesCOVIDShuffled,:),metric_reg_age(indicesControlShuffled,:));
    tstatsReg(j,:) = stats.tstat;
    pstatsReg(j,:) = p_aux;
end
% % Visually checking the results
% for i = 1 : size(diffMeanGroups,2)
%     figure;
%     set(gcf, 'Position', [100 100 1600 1200])
%     subplot(1,2,1);
%     hist(diffMeanGroupsPermReg(:,i));
%     hold on; plot(diffMeanGroups(i),0,'x');
%     subplot(1,2,2);
%     hist(diffMeanGroupsPermReg(:,i));
%     hold on; plot(diffMeanGroupsReg(i),0,'x');
%     pause
% end
tstats_max_perm = max(tstats,[],2); % max tstats for each permutation
tstats_min_perm = min(tstats,[],2); % min tstats for each permutation
tstats_max_perm_reg = max(tstatsReg,[],2); % max tstats for each permutation
tstats_min_perm_reg = min(tstatsReg,[],2); % min tstats for each permutation

% Get number of values lower than the tstats_max_perm for each of the
% original values:
num_greater = sum(tstats_max_perm>tstats_orig,1);
num_lower = sum(tstats_min_perm<tstats_orig,1);
num_greater_reg = sum(tstats_max_perm_reg >= tstats_orig,1);
num_lower_reg = sum(tstats_min_perm_reg <= tstats_orig,1);
p_greater = num_greater./ numPermutations;
p_lower = num_lower ./ numPermutations;
p_greater_reg = num_greater_reg./ numPermutations;
p_lower_reg = num_lower_reg ./ numPermutations;

% Westfall and young permutation
[p_orig_sorted, indices_p_orig_sorted] = sort(p_t_orig);
[p_reg_orig_sorted, indices_p_reg_orig_sorted] = sort(p_t_reg_orig);
pstats_sorted = sort(pstats,2);
pstatsReg_sorted = sort(pstatsReg,2);
tstats_sorted = sort(tstats,2);
tstatsReg_sorted = sort(tstatsReg,2);
[tstats_orig_sorted, indices_tstats_orig_sorted] = sort(tstats_orig);
[tstats_reg_orig_sorted, indices_tstats_reg_orig_sorted] = sort(tstatsReg_orig);
% For each ROI compare the ith highes permuted statistics:
for i = 1 : size(matSegmentationsPerSubject,2)
    perm_min = tstats_sorted(:, i);
    permReg_min = tstatsReg_sorted(:, i);
    pPerm_min = pstats_sorted(:, i);
    pRegPerm_min = pstatsReg_sorted(:, i);
    %adjusted_p_values(indices_tstats_orig_sorted(i)) = mean(perm_min <= tstats_orig_sorted(i));
    adjusted_p_values_sorted(i) = mean(perm_min <= tstats_orig_sorted(i));
    adjusted_p_values2_sorted(i) = mean(pPerm_min <= p_orig_sorted(i));
    adjusted_p_values_reg_sorted(i) = mean(permReg_min <= tstats_orig_sorted(i));
    adjusted_p_values2_reg_sorted(i) = mean(pRegPerm_min <= p_reg_orig_sorted(i));
end
csvwrite(fullfile(resultsPath, 'corticalThikness_pValues.csv'), p_t_orig);
csvwrite(fullfile(resultsPath, 'corticalThiknessReg_pValues.csv'), p_t_reg_orig);
csvwrite(fullfile(resultsPath, 'corticalThiknessPermutations_pValues.csv'), pstats);
csvwrite(fullfile(resultsPath, 'corticalThiknessPermutationReg_pValues.csv'), pstatsReg);

csvwrite(fullfile(resultsPath, 'corticalThikness_tstats.csv'), tstats_orig);
csvwrite(fullfile(resultsPath, 'corticalThiknessReg_tstats.csv'), tstatsReg_orig);
csvwrite(fullfile(resultsPath, 'corticalThiknessPermutations_tstats.csv'), tstats);
csvwrite(fullfile(resultsPath, 'corticalThiknessPermutationReg_tstats.csv'), tstatsReg);


for i = 2 : size(matSegmentationsPerSubject,2)
    if adjusted_p_values_sorted(i) < adjusted_p_values_sorted(i-1)
        adjusted_p_values_sorted(i) = adjusted_p_values_sorted(i - 1);
    end
    if adjusted_p_values_reg_sorted(i) < adjusted_p_values_reg_sorted(i-1)
        adjusted_p_values_reg_sorted(i) = adjusted_p_values_reg_sorted(i - 1);
    end
    if adjusted_p_values2_sorted(i) < adjusted_p_values2_sorted(i-1)
        adjusted_p_values2_sorted(i) = adjusted_p_values2_sorted(i - 1);
    end
    if adjusted_p_values2_reg_sorted(i) < adjusted_p_values2_reg_sorted(i-1)
        adjusted_p_values2_reg_sorted(i) = adjusted_p_values2_reg_sorted(i - 1);
    end
end
index_p_perm_critical = find(adjusted_p_values_sorted>0.05);
p_perm_critical = adjusted_p_values_sorted(index_p_perm_critical(1));
index_p_perm_critical2 = find(adjusted_p_values2_sorted>0.05);
p_perm_critical2 = adjusted_p_values2_sorted(index_p_perm_critical2(1));
index_p_reg_perm_critical = find(adjusted_p_values_reg_sorted>0.05);
p_reg_perm_critical = adjusted_p_values_reg_sorted(index_p_reg_perm_critical(1));
index_p_reg_perm_critical2 = find(adjusted_p_values2_reg_sorted>0.05);
p_reg_perm_critical2 = adjusted_p_values2_reg_sorted(index_p_reg_perm_critical2(1));

adjusted_p_values = adjusted_p_values_sorted(indices_tstats_orig_sorted);
adjusted_p_values2 = adjusted_p_values2_sorted(indices_p_orig_sorted);
adjusted_p_reg_values = adjusted_p_values_reg_sorted(indices_tstats_reg_orig_sorted);
adjusted_p_reg_values2 = adjusted_p_values2_reg_sorted(indices_tstats_reg_orig_sorted);

indices_sign_perm_WY = find(adjusted_p_values < p_perm_critical);
indices_sign_perm_2_WY = find(adjusted_p_values2 < p_perm_critical2);
indices_sign_perm_reg_WY = find(adjusted_p_reg_values < p_reg_perm_critical);
indices_sign_perm_reg_2_WY = find(adjusted_p_values2_reg_sorted < p_reg_perm_critical2);

parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_WY)
parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_2_WY)
parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg_WY)
parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg_2_WY)

% compute p values for the measured difference:
p_perm = sum(diffMeanGroups > diffMeanGroupsPerm,1)./numPermutations;
p_perm_reg = sum(diffMeanGroupsReg > diffMeanGroupsPermReg,1)./numPermutations;
% Z scores
z_scores = (diffMeanGroups - mean(diffMeanGroupsPerm,1))./std(diffMeanGroupsPerm, [], 1);
z_scores_reg = (diffMeanGroupsReg - mean(diffMeanGroupsPermReg,1))./std(diffMeanGroupsPermReg, [], 1);
% p values norm dist
p_perm_norm = normcdf(diffMeanGroups, mean(diffMeanGroupsPerm,1), std(diffMeanGroupsPerm, [], 1));
p_perm_norm_reg = normcdf(diffMeanGroupsReg, mean(diffMeanGroupsPermReg,1), std(diffMeanGroupsPermReg, [], 1));


% compute p values for the measured difference:
p_perm = sum(diffMeanGroups > diffMeanGroupsPerm,1)./numPermutations;
p_perm_reg = sum(diffMeanGroupsReg > diffMeanGroupsPermReg,1)./numPermutations;
% Z scores
z_scores = (diffMeanGroups - mean(diffMeanGroupsPerm,1))./std(diffMeanGroupsPerm, [], 1);
z_scores_reg = (diffMeanGroupsReg - mean(diffMeanGroupsPermReg,1))./std(diffMeanGroupsPermReg, [], 1);
% p values norm dist
p_perm_norm = normcdf(diffMeanGroups, mean(diffMeanGroupsPerm,1), std(diffMeanGroupsPerm, [], 1));
p_perm_norm_reg = normcdf(diffMeanGroupsReg, mean(diffMeanGroupsPermReg,1), std(diffMeanGroupsPermReg, [], 1));


% COMPARE ALL p
p
p_perm
% Significant:
indices_sign_perm = find(p_perm < 0.05)
indices_sign_perm_reg = find(p_perm_reg < 0.05)
indices(indices_sign)
indices(indices_sign2)

% Save results
resultsParcellationThick.no_mc.p_values = p;
resultsParcellationThick.no_mc.p_values_ord = p_ordered;
resultsParcellationThick.no_mc.indices_sign = find(p<0.05);
resultsParcellationThick.no_mc.p_sign = p(find(p<0.05));
resultsParcellationThick.no_mc.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(find(p<0.05));
% fdr
resultsParcellationThick.fdr.p_th= p_corr;
resultsParcellationThick.fdr.indices_sign = indices(indices_sign);
resultsParcellationThick.fdr.p_sign = p(indices(indices_sign));
resultsParcellationThick.fdr.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(indices(indices_sign));

% holm
resultsParcellationThick.holm.p_th = p_thresh;
resultsParcellationThick.holm.indices_sign = indices(indices_sign2);
resultsParcellationThick.holm.p_sign = p(indices(indices_sign2));
resultsParcellationThick.holm.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(indices(indices_sign2));

% permutation
resultsParcellationThick.permutation.p_values = p_perm;
resultsParcellationThick.permutation.diff_means_perm = diffMeanGroups;
resultsParcellationThick.permutation.indices_sign = indices_sign_perm;
resultsParcellationThick.permutation.p_sign = p_perm(indices_sign_perm);
resultsParcellationThick.permutation.z_scores = z_scores;
resultsParcellationThick.permutation.p_norm = p_perm_norm;
resultsParcellationThick.permutation.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm);

% permutation
resultsParcellationThick.permutationWestfallAndYoung.p_values = adjusted_p_values;
resultsParcellationThick.permutationWestfallAndYoung.indices_sign = find(adjusted_p_values<0.05);
resultsParcellationThick.permutationWestfallAndYoung.p_sign = p_perm(adjusted_p_values<0.05);
resultsParcellationThick.permutationWestfallAndYoung.tstats = tstats;
resultsParcellationThick.permutationWestfallAndYoung.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(adjusted_p_values<0.05);

% permutation with variables regressed out
resultsParcellationThick.permutation_reg.p_values = p_perm_reg;
resultsParcellationThick.permutation_reg.adjusted_p_values = adjusted_p_reg_values;
resultsParcellationThick.permutation_reg.diff_means_perm = diffMeanGroupsPermReg;
resultsParcellationThick.permutation_reg.indices_sign = indices_sign_perm_reg;
resultsParcellationThick.permutation_reg.p_sign = p_perm_reg(indices_sign_perm_reg);
resultsParcellationThick.permutation_reg.z_scores = z_scores_reg;
resultsParcellationThick.permutation_reg.p_norm = p_perm_norm_reg;
resultsParcellationThick.permutation_reg.names = parcellation(:,indicesRoisToUse).Properties.VariableNames(indices_sign_perm_reg);

save([resultsPath 'parcellation_thick_multiple_comparisons.mat'], 'resultsParcellationThick')

% compare p values and regions:
[p' p_perm' adjusted_p_values']

[p(indices_sign_perm)' p_perm(indices_sign_perm)' adjusted_p_values(indices_sign_perm)']
resultsParcellationThick.permutation.names'
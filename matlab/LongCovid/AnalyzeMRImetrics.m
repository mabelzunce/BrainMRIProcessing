clear all
close all

dataPartitionPath = '/data/'; %'D:/'
adniPartitionPath = '/data_imaging/'; %'F:/'
% Load data:
loadData = 0;

%% PATHS AND FILENAMES
studyDataPath = '../';
resultsPath = '../DataAnalysis/';
summaryExcelFilename = fullfile(studyDataPath, 'RespuestasCuestionarioEvaluacionCognitivaResonancia.xlsx');
volunteersExcelFilename = fullfile(studyDataPath, 'VoluntariosProyectoCovidProlongado.xlsx');
fslPreprocessedDataPath = '/home/martin/data/UNSAM/CovidProject/Estudio/PreprocessedMRI/FSL/';
freesurferBrainVolumesPath = '/home/martin/data/UNSAM/CovidProject/Estudio/DataAnalysis/AnatomicalMRI/freesurfer/';
outputPath = [resultsPath 'BrainMorphometry/'];
if ~exist(outputPath)
    mkdir(outputPath);
end
%% DATA DEFINITIONS 
groupNames = {'Control', 'Covid Prolongado'};
groupNamesForPlots = {'Controls', 'Long COVID'};
shortBroupNamesForPlots = {'C', 'LC'};
%% COLUMNS ORDER
% Groups:
colsDemographicData = 4:7;
colsCovidData = 8:36;
colsLongCovidSymptoms = 24:35;
colsHealthHistory = 37:54;
% Individual:
colId = 1; colGroup = 2; colAge = 4; colGender = 5; colHeight = 6; colWeight = 7; colVaccinated = 8; colNumVaccines = 9; 
colDateLastVaccine = 10; colCovidInfection = 11; colNumCovidInfections = 12; colDateFirstInfection = 13; colNumVaccinesCovidInfection = 14;
colDateLastSymptomsCovid = 15; colAdmitted = 16; colDateAdmitted = 17; colDateDischarged = 18; colAdmittedAfterCovid = 19; colNumAdmittedAfterCovid = 20;
colAdmittedICU = 21; colFeelRecovered = 22; colDiagnosedWithLongCovid = 23;
% Persisten symptoms post Covid:
colHeadache = 24; colFatigue = 25; colAnosmia = 26; colLossTaste = 27; colDisnea = 28; colMuscleWeakness = 29; colMusclePain = 30;
colFocusConfussion = 31; colTalkingComm = 32; colSleeping = 33; colMemory = 34; colAttention = 35; colOtherSymptoms = 36;
% Health background:
colHypertension = 37; colDiabetes = 38; colAsma = 39; colCholesterol = 40; colHeartAttack = 41;
colAngina = 42; colThrombosis = 43; colBloodPressureMedication = 44; colAspirine = 45; colBloodThinner = 46;
colOtherHealthIssues = 47; colOtherHealthIssuesDescr = 48; colSmoker = 49; colCigarettesPerDay = 50; colYearsSmoking = 51;
colYearsExSmoker = 52; colDrinksAlcohol = 53; colAlcoholUnitsPerWeek = 54;

% FAS QUEATIONNAIRE
colFAS = 75;
%% UBJECTS TO EXCLUDE:
% Apply some filters for unwanted volunteers:
subjectsToExclude = {'CP0011', 'CP0015'}; % 11 (foreign language) and 15 (elder and probaly with some form of dementia).
%% VOLUNTEERS TABLE
%volunteersTable = readtable(volunteersExcelFilename, 'NumHeaderLines', 1);
volunteersTable = readtable(volunteersExcelFilename, 'Range', 'A2:AI168');%, 'HeaderLines', 1);
% Remove data to exclude:
indicesToExclude = [];
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, volunteersTable.ID)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind];
    end
end
volunteersTable(indicesToExclude,:) = [];
subjectNamesAllVolunteers = volunteersTable.ID;
groupAllVolunteers = volunteersTable.Grupo;
indicesControlVolunteers = find(strcmp(volunteersTable.Grupo, groupNames{1}));
indicesCovidVolunteers = find(strcmp(volunteersTable.Grupo, groupNames{2}));
for i = 1 :  numel(groupNames)
    indicesGroupsVolunteers{i} = find(strcmp(volunteersTable.Grupo, groupNames{i}));
end

%% SUMMARY VOLUNTEERS
dateCuestionarios = datetime(volunteersTable.Cuestionarios);
dateResonancia = datetime(volunteersTable.Resonancia);
dateCuestionarios = datetime(volunteersTable.Cuestionarios);
indicesWithQuestionnaire = ~isnat(dateCuestionarios);
indicesWithMri = ~isnat(dateResonancia);
indicesWithQuestAndMri = indicesWithQuestionnaire & indicesWithMri;
numCovidWithQuest = sum(strcmp(volunteersTable.Grupo, 'Covid Prolongado') & indicesWithQuestionnaire);
numCovidWithMri = sum(strcmp(volunteersTable.Grupo, 'Covid Prolongado') & indicesWithMri);
numCovidWithQuestAndMri = sum(strcmp(volunteersTable.Grupo, 'Covid Prolongado') & indicesWithQuestAndMri);
numControlWithQuest = sum(strcmp(volunteersTable.Grupo, 'Control') & indicesWithQuestionnaire);
numControlWithMri = sum(strcmp(volunteersTable.Grupo, 'Control') & indicesWithMri);
numControlWithQuestAndMri = sum(strcmp(volunteersTable.Grupo, 'Control') & indicesWithQuestAndMri);
% Estadistica basada en estados:
statusLabels = {'Baja','Completo','Solo Cognitivo','Solo Resonancia','Turnos Pendientes'};
[countsStatus, labels] = hist(categorical(volunteersTable.Estado), statusLabels);
% Stats:
numRecruited = sum(countsStatus);
numCompleted = countsStatus(2);
numDropOut = countsStatus(1);
numOnlyCogn = countsStatus(3);
numOnlyMRI = countsStatus(4);
numPendingAppointments = countsStatus(5);
% See more what is going on with the pending appointments:
subjectsWithPendingAppointments = strcmp(volunteersTable.Estado, 'Turnos Pendientes');
subjectsWithPendingMRI = subjectsWithPendingAppointments & indicesWithQuestionnaire & ~indicesWithMri;
subjectsWithPendingCog = subjectsWithPendingAppointments & ~indicesWithQuestionnaire & indicesWithMri;
subjectsWithPendingAll = subjectsWithPendingAppointments & ~indicesWithQuestionnaire & ~indicesWithMri;
numSubjectsWithPendingMRI = sum(subjectsWithPendingMRI);
numSubjectsWithPendingCog = sum(subjectsWithPendingCog);
numSubjectsWithPendingAll = sum(subjectsWithPendingAll);
% Subjects by group:
subjectsControlCompleted = strcmp(volunteersTable.Grupo, 'Control') & strcmp(volunteersTable.Estado, 'Completo');
subjectsLongCovidCompleted = strcmp(volunteersTable.Grupo, 'Covid Prolongado') & strcmp(volunteersTable.Estado, 'Completo');
% Show the summary:
disp('Based on the volunteers Excel')
fprintf('Number of complete subjects: %d.\n', numCompleted);
fprintf('Number of drop outs: %d.\n', numDropOut);
fprintf('Number of subjects with only cognitive data: %d.\n', numOnlyCogn);
fprintf('Number of subjects with only MRI: %d.\n', numOnlyMRI);
%% LOAD MRI REPORTS EXCEL
%summaryTable = readtable(summaryExcelFilename, 'NumHeaderLines', 3);
summaryTable = readtable(summaryExcelFilename,'Sheet', 'respuestasCuestionario', 'Range', 'A4:EC144', 'ReadVariableNames', false, 'ReadRowNames', false);
% Remove data to exclude:
indicesToExclude = [];
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, summaryTable.Var1)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind];
    end
end
summaryTable(indicesToExclude,:) = [];

% Get fields:
subjectNames = summaryTable.Var1;
group = summaryTable.Var2;
age_years = summaryTable.Var4;
gender = summaryTable.Var5;
height_cm = summaryTable.Var6;
weight_kg = summaryTable.Var7;
vaxCovid_si_no = summaryTable.Var8;
numVaxDoses = summaryTable.Var9;
covid_si_no = summaryTable.Var11;
numCovidInfections = summaryTable.Var12;
dateFirstInfection = summaryTable.Var13;
admittedHospForCovid = summaryTable.Var16;
dateAdmitHospForCovid = summaryTable.Var17;
dateDischargeHospForCovid = summaryTable.Var18;
readmittedHospAfterCovid = summaryTable.Var19;

% Symptoms:
symptom = table2array(summaryTable(:,colsLongCovidSymptoms));
symptomsNames = {'Headaches', 'Fatigue', 'Anosmia', 'Loss of Taste', ...
    'Dyspnoea', 'Muscle weakness', 'Muscle ache', 'Brain fog', 'Communication problems',...
    'Sleeping problems', 'Memory problems', 'Problems with attention'};
labelsSymptoms = {'Yes, and they persist', 'Yes, but no anymore', 'No'};
%% LOAD MRI REPORTS AND COGNTIVE DATA EXCEL
%cognitiveTestsTable = readtable(summaryExcelFilename,'Sheet','EvaluacionCognitiva', 'NumHeaderLines', 2);
%mriReportsTable = readtable(summaryExcelFilename,'Sheet','InformeResonancia', 'NumHeaderLines', 0);
cognitiveTestsTable = readtable(summaryExcelFilename,'Sheet','EvaluacionCognitiva', 'Range', 'A3:AG141', 'ReadVariableNames', false, 'ReadRowNames', false);
mriReportsTable = readtable(summaryExcelFilename,'Sheet','InformeResonancia');
% Remove data to exclude:
indicesToExclude = [];
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, cognitiveTestsTable.Var1)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind];
    end
end
cognitiveTestsTable(indicesToExclude,:) = [];
% Remove data to exclude:
indicesToExclude = [];
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, mriReportsTable.ID)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind];
    end
end
mriReportsTable(indicesToExclude,:) = [];

groupTests = cognitiveTestsTable.Var2;
% Trail making test A
tmtA_score = cognitiveTestsTable.Var3;
tmtA_ref = cognitiveTestsTable.Var4;
tmtA_performance = cognitiveTestsTable.Var5;
% Trail making test B
tmtB_score = cognitiveTestsTable.Var6;
tmtB_ref = cognitiveTestsTable.Var7;
tmtB_performance = cognitiveTestsTable.Var8;
% WMS-R digit span forward
digitSpanForward_score = cognitiveTestsTable.Var9;
digitSpanForward_perc = cognitiveTestsTable.Var10;
digitSpanForward_performance = cognitiveTestsTable.Var11;
% WMS-R digit span backwards
digitSpanBackward_score = cognitiveTestsTable.Var12;
digitSpanBackward_perc = cognitiveTestsTable.Var13;
digitSpanBackward_performance = cognitiveTestsTable.Var14;
% Stroop P
stroopP_score = cognitiveTestsTable.Var15;
stroopP_ref = cognitiveTestsTable.Var16;
stroopP_performance = cognitiveTestsTable.Var17;
% Stroop C
stroopC_score = cognitiveTestsTable.Var18;
stroopC_ref = cognitiveTestsTable.Var19;
stroopC_performance = cognitiveTestsTable.Var20;
% Stroop P/C
stroopPC_score = cognitiveTestsTable.Var21;
stroopPC_ref = cognitiveTestsTable.Var22;
stroopPC_performance = cognitiveTestsTable.Var23;
% Stroop Interference
stroopInterference_score = cognitiveTestsTable.Var24;
stroopInterference_ref = cognitiveTestsTable.Var25;
stroopInterference_performance = cognitiveTestsTable.Var26;
% MOCA
moca_score = cognitiveTestsTable.Var27;
moca_ref = cognitiveTestsTable.Var28;
moca_performance = cognitiveTestsTable.Var29;

% SUMMARY IN A TABLE:
allCognitiveTests = cognitiveTestsTable(:,3:29);
allCognitiveTests_score = cognitiveTestsTable(:,3:3:29);
allCognitiveTests_ref = cognitiveTestsTable(:,4:3:29);
allCognitiveTests_perf = cognitiveTestsTable(:,5:3:29);
subjectNamesTests = cognitiveTestsTable.Var1;

subjectNamesMriReports = mriReportsTable.ID;
subjectNamesMriReportsWithMri = mriReportsTable.ID(strcmp(mriReportsTable.SeHizoResonancia,'Si'));
subjectNamesMriReportsWithReport = mriReportsTable.ID(strcmp(mriReportsTable.ReporteResonancia,'Si'));
%% CHECK IF DATA FOR QUESTIONNARIES IS THE SAME THAN FOR COGNTIIVE TESTS
% Get indices Questionnaires with data (using one field):
indicesSubjectsWithQ = find(~isnan(age_years)>0);
indicesSubjectsWithCT = find(~isnan(tmtA_score)>0);
subjectNamesWithQ = subjectNames(indicesSubjectsWithQ);
subjectNamesWithCT = subjectNamesTests(indicesSubjectsWithCT);
matchedLists = 1;
if numel(subjectNamesWithQ) ~= numel(subjectNamesWithCT)
    matchedLists = 0;
elseif sum(strcmp(subjectNamesWithQ, subjectNamesWithCT)) ~= numel(subjectNamesWithQ)
    matchedLists = 0;
end
if matchedLists == 0
    warning('The list of subjects with questionnaire data and cognitive tests does not match. We will try to match the data');
    indicesBothQ = [];
    indicesBothCT = [];
    for i = 1 : numel(subjectNamesWithQ)     
        indexCT = find(strcmp(subjectNamesWithQ{i}, subjectNamesWithCT));
        if ~isempty(indexCT)
            % Add it:
            indicesBothQ = [indicesBothQ i];
            indicesBothCT = [indicesBothCT indexCT];
        end
    end
    indicesSubjectsWithQ = indicesSubjectsWithQ(indicesBothQ);
    indicesSubjectsWithCT = indicesSubjectsWithCT(indicesBothCT);
    subjectNamesWithQ = subjectNames(indicesSubjectsWithQ);
    subjectNamesWithCT = subjectNamesTests(indicesSubjectsWithCT);
    if sum(strcmp(subjectNamesWithQ, subjectNamesWithCT)) ~= numel(subjectNamesWithQ)
        error('Lists could not be matched')
    else
        disp('Data matched.');
    end
end
% Get groups with tests:
indicesCovidQ = find(strcmp(group, 'COVID')>0);
indicesControlQ = find(strcmp(group, 'CONTROL')>0);
indicesCovidCT = find(strcmp(groupTests, 'COVID')>0);
indicesControlCT = find(strcmp(groupTests, 'CONTROL')>0);
indicesGroupsQ = {indicesControlQ, indicesCovidQ};
indicesGroupsCT = {indicesControlCT, indicesCovidCT};

%% MRI REPORTS
indicesSubjectsWithMriReport = strcmp(mriReportsTable.ReporteResonancia,'Si');
for i = 1 : numel(subjectNamesMriReportsWithReport)
    indexThisSubjectInVolunteers = find(strcmp(subjectNamesMriReports{i}, volunteersTable.ID));
    %indicesSubjectsWithMriReport(i) = indexThisSubjectInVolunteers;
    groupsSubjectsWithMriReport(i) = volunteersTable.Grupo(indexThisSubjectInVolunteers);
end
reports = categorical(mriReportsTable.Impresi_nDiagn_stica_(indicesSubjectsWithMriReport));

legendsForPlots = {'Incidental findings that do not requiere foloow-up', 'Incidental findings that requiere follow-up', ...
    'No findings'};
figure;
set(gcf, 'Position', [0 0 1600, 700]);
tiledlayout(1,2,'TileSpacing','compact');
for i = 1 : numel(groupNames)
    indicesThisGroup = strcmp(groupsSubjectsWithMriReport, groupNames{i});
    [counts, cat] = hist(reports(indicesThisGroup));
    nexttile
    pie(counts,'%.1f%%');
    %hl = piechart(counts);%,'%.1f%%');
    title(groupNamesForPlots{i})
    %hl.Labels(counts == 0) = '';
    % % Create legend
    % for i = 1 : numel(cat)
    %     hl.Names(i) = cat{i};
    % end
    % hl.LegendVisible = 1;
    if i == numel(groupNames)
        hl = legend(legendsForPlots, 'Location', 'South');
    end
end
hl.Position(2) = hl.Position(2) - 0.12;

%% SIENAX
indicesSubjectsWithSienax = [];
sienaxResults = readtable([fslPreprocessedDataPath '/SienaxResults2.csv']);
% Remove data to exclude:
indicesToExclude = [];
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, sienaxResults.SubjectNames)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind];
    end
end
sienaxResults(indicesToExclude,:) = [];
% Get the indices from the total population:
for i = 1 : numel(sienaxResults.SubjectNames)    
    % Using volunteers table:
    indexThisSubjectInVolunteers = find(strcmp(sienaxResults.SubjectNames{i}, volunteersTable.ID));
    groupsSubjectsWithSienax(i) = volunteersTable.Grupo(indexThisSubjectInVolunteers);
    indicesSubjectsWithSienax = [indicesSubjectsWithSienax; indexThisSubjectInVolunteers];

%     % Using questionnaires:
%     indexThisSubjectInVolunteers = find(strcmp(sienaxResults.SubjectNames{i}, subjectNames));
%     if ~isempty(indexThisSubjectInVolunteers)
%         groupsSubjectsWithSienax(i) = group(indexThisSubjectInVolunteers);
%         %indicesSubjectsWithMriReport(i) = indexThisSubjectInVolunteers;
%         indicesSubjectsWithSienax = [indicesSubjectsWithSienax; indexThisSubjectInVolunteers];
%     else
%         warning(sprintf('Data from %s is missing.\n', sienaxResults.SubjectNames{i}));
%     end
end

% Check for outliers to check the measurements:
indicesOutlierGMV = isoutlier(sienaxResults.GreyMatterVolume);
indicesOutlierWMV = isoutlier(sienaxResults.WhiteMatterVolume);
indicesOutlierTBV = isoutlier(sienaxResults.BrainVolume);
sienaxResults.SubjectNames(indicesOutlierGMV)
sienaxResults.SubjectNames(indicesOutlierTBV)
% Normalized statistical tests:
[pKW_normGMV, annovaTab_normGMV, stats_normGMV] = kruskalwallis(sienaxResults.GreyMatterVolume, groupsSubjectsWithSienax);
[pKW_normWMV, annovaTab_normWMV, stats_normWMV] = kruskalwallis(sienaxResults.WhiteMatterVolume, groupsSubjectsWithSienax);
[pKW_normBV, annovaTab_normBV, stats_normBV] = kruskalwallis(sienaxResults.BrainVolume, groupsSubjectsWithSienax);
[pKW_unnormGMV, annovaTab_unnormGMV, stats_unnormGMV] = kruskalwallis(sienaxResults.GreyMatterUnNormVolume, groupsSubjectsWithSienax);
[pKW_unnormWMV, annovaTab_unnormWMV, stats_unnormWMV] = kruskalwallis(sienaxResults.WhiteMatterUnNNormVolume, groupsSubjectsWithSienax);
[pKW_unnormUBV, annovaTab_unnormBV, stats_unnormBV] = kruskalwallis(sienaxResults.BrainUnNNormVolume, groupsSubjectsWithSienax);
tableSienax = table();
tableSienax.Metrics = cell2table({'Grey Matter Volume', 'White Matter Volume', 'Total Volume'}');
tableSienax.p_values = [pKW_normGMV; pKW_normWMV; pKW_normBV];
close all
%% MULTIVARIATE LINER REGRESSION FOR GMV
j = 0;
% I need the questionnaire data for this:
indicesSubjectsWithSienaxAndQ = [];
for i = 1 : numel(sienaxResults.SubjectNames)  
    % Using questionnaires:
    indexThisSubjectInVolunteers = find(strcmp(sienaxResults.SubjectNames{i}, subjectNames));
    if ~isempty(indexThisSubjectInVolunteers)
        j = j + 1;
        indicesSienaxWithW(j) = i;
        groupsSubjectsWithSienaxAndQ(j) = group(indexThisSubjectInVolunteers);
        %indicesSubjectsWithMriReport(i) = indexThisSubjectInVolunteers;
        indicesSubjectsWithSienaxAndQ = [indicesSubjectsWithSienaxAndQ; indexThisSubjectInVolunteers];
    else
        warning(sprintf('Data from %s is missing.\n', sienaxResults.SubjectNames{i}));
    end
end
% Remove subjects with missing data:
indicesToRemove = find(isnan(age_years(indicesSubjectsWithSienaxAndQ)) | ...
    isnan(tmtA_score(indicesSienaxWithW))>0);
indicesSubjectsWithSienaxAndQ(indicesToRemove) = [];
indicesSienaxWithW(indicesToRemove) = [];
ageSienax_years = age_years(indicesSubjectsWithSienaxAndQ); 
genderSienax_years = gender(indicesSubjectsWithSienaxAndQ); 
groupCategorical = categorical( group(indicesSubjectsWithSienaxAndQ));

tmtASienax_score = tmtA_score(indicesSienaxWithW);
tmtASienax_ref = tmtA_ref(indicesSienaxWithW);
tmtASienax_performance = categorical(tmtA_performance(indicesSienaxWithW));
mocaSienax_score = moca_score(indicesSienaxWithW);
mocaSienax_ref = moca_ref(indicesSienaxWithW);

gmVolume = sienaxResults.GreyMatterVolume(indicesSienaxWithW);
wmVolume = sienaxResults.WhiteMatterVolume(indicesSienaxWithW);
brainVolume = sienaxResults.BrainVolume(indicesSienaxWithW);


tbl = table(ageSienax_years, genderSienax_years, groupCategorical, tmtASienax_score, tmtASienax_ref, ...
    tmtASienax_performance, mocaSienax_score, mocaSienax_ref, gmVolume); % Last column, response variable.
mdlStep = stepwiselm(tbl, 'linear', 'Verbose',2);%'MPG ~ Weight'

mdl0 = fitlm(tbl, 'gmVolume~ageSienax_years');%+genderSienax_years');
variables = mdl0.Coefficients.Variables;
gmVolumeAgeCorr = gmVolume - variables(1,1) - variables(2,1)*ageSienax_years;

mdl1 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical');
variables = mdl1.Coefficients.Variables;

mdl2 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical+mocaSienax_score');
variables = mdl2.Coefficients.Variables;

mdl3 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical+mocaSienax_ref');
variables = mdl3.Coefficients.Variables;

mdl4 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical+tmtASienax_score');
variables = mdl4.Coefficients.Variables;

mdl5 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical+tmtASienax_ref');
variables = mdl4.Coefficients.Variables;

mdl6 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical+tmtASienax_performance');
variables = mdl4.Coefficients.Variables;

figure;
plot(ageSienax_years(groupCategorical == 'CONTROL'), gmVolume(groupCategorical == 'CONTROL'), 's', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
plot(ageSienax_years(groupCategorical == 'COVID'), gmVolume(groupCategorical == 'COVID'), 'o', 'LineWidth', 2, 'MarkerSize', 10);


legend(groupNamesForPlots);
saveas(gca, [resultsPath 'SienaxGreyMatterVolume.png']);
%% FREESURFER
% Process table with a summary of all the results (created by Sol script):
freesurferBrainVolumesCsv = [freesurferBrainVolumesPath 'brain_volumes.csv'];
freesurferSegmentationCsv = [freesurferBrainVolumesPath 'segmentation.csv'];
freesurferParcCsv = [freesurferBrainVolumesPath 'parcellation.csv'];
freesurferParcThickCsv = [freesurferBrainVolumesPath 'parcellation_thick.csv'];
freesurferParcDKTCsv = [freesurferBrainVolumesPath 'parcellation_DKT.csv'];
freesurferParcThickDKTCsv = [freesurferBrainVolumesPath 'parcellation_thick_DKT.csv'];

freesurferBrainVolumesTable = readtable(freesurferBrainVolumesCsv);
freesurferSegmentationTable = readtable(freesurferSegmentationCsv);
freesurferParcTable = readtable(freesurferParcCsv);
freesurferParcThickTable = readtable(freesurferParcThickCsv);
freesurferParcDKTTable = readtable(freesurferParcDKTCsv);
freesurferParcThickDKTTable = readtable(freesurferParcThickDKTCsv);

% Check if all the names agree
freesurferSubjects = freesurferBrainVolumesTable.subject;
if (strcmp(freesurferSegmentationTable.subject,freesurferSubjects) == 0) > 0
    error('Missmatch between tables');
end
% Parcellations have two entreies for subejct
if numel(freesurferParcTable.subject(1:2:end)) ~= numel(freesurferSubjects)    
    for i = 1 : 2: numel(freesurferParcTable.subject)
        if sum(strcmp(freesurferParcTable.subject(i),freesurferSubjects)) == 0
            warning(sprintf('Removing subject %s', freesurferParcTable.subject{i}))
        end
    end
elseif (strcmp(freesurferParcTable.subject(1:2:end), freesurferSubjects) == 0) > 0
    error('Missmatch between tables');
end
if (strcmp(freesurferParcThickTable.subject(1:2:end),freesurferSubjects) == 0) > 0
    error('Missmatch between tables');
end
if (strcmp(freesurferParcDKTTable.subject(1:2:end),freesurferSubjects) == 0) > 0
    error('Missmatch between tables');
end
if (strcmp(freesurferParcThickDKTTable.subject(1:2:end),freesurferSubjects) == 0) > 0
    error('Missmatch between tables');
end

% Remove data to exclude:
indicesToExclude = [];
for i = 1 : numel(subjectsToExclude)
    ind=find(strcmp(subjectsToExclude{i}, freesurferSubjects)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind];
    end
end
indicesToExcludeTwoSides = [indicesToExclude*2-1 indicesToExclude*2];
freesurferSubjects(indicesToExclude) = [];
freesurferBrainVolumesTable(indicesToExclude,:) = [];
freesurferSegmentationTable(indicesToExclude,:) = [];
freesurferParcTable(indicesToExcludeTwoSides,:) = [];
freesurferParcThickTable(indicesToExcludeTwoSides,:) = [];
freesurferParcDKTTable(indicesToExcludeTwoSides,:) = [];
freesurferParcThickDKTTable(indicesToExcludeTwoSides,:) = [];

% Parcellations have entries for left and right:
freesurferParcLeftTable = freesurferParcTable(strcmp(freesurferParcTable.Hemisphere, 'lh'),:);
freesurferParcRightTable = freesurferParcTable(strcmp(freesurferParcTable.Hemisphere, 'rh'),:);
freesurferParcThickLeftTable = freesurferParcThickTable(strcmp(freesurferParcThickTable.Hemisphere, 'lh'),:);
freesurferParcThickRightTable = freesurferParcThickTable(strcmp(freesurferParcThickTable.Hemisphere, 'rh'),:);
freesurferParcDKTLeftTable = freesurferParcDKTTable(strcmp(freesurferParcDKTTable.Hemisphere, 'lh'),:);
freesurferParcDKTRightTable = freesurferParcDKTTable(strcmp(freesurferParcDKTTable.Hemisphere, 'rh'),:);
freesurferParcThickDKTLeftTable = freesurferParcThickDKTTable(strcmp(freesurferParcThickDKTTable.Hemisphere, 'lh'),:);
freesurferParcThickDKTRightTable = freesurferParcThickDKTTable(strcmp(freesurferParcThickDKTTable.Hemisphere, 'rh'),:);
% Reduce tables to numeric data:
freesurferBrainVolumes = table2array(freesurferBrainVolumesTable(:,2:end-2));
freesurferSegmentation = table2array(freesurferSegmentationTable(:,2:end-2));
freesurferParcLeft = table2array(freesurferParcLeftTable(:,2:end-3));
freesurferParcRight = table2array(freesurferParcRightTable(:,2:end-3));
freesurferParcThickLeft = table2array(freesurferParcThickLeftTable(:,2:end-3));
freesurferParcThickRight = table2array(freesurferParcThickRightTable(:,2:end-3));
freesurferParcDKTLeft = table2array(freesurferParcDKTLeftTable(:,2:end-3));
freesurferParcDKTRight = table2array(freesurferParcDKTRightTable(:,2:end-3));
freesurferParcThickDKTLeft = table2array(freesurferParcThickDKTLeftTable(:,2:end-3));
freesurferParcThickDKTRight = table2array(freesurferParcThickDKTRightTable(:,2:end-3));

freesurferParcMean = (freesurferParcLeft + freesurferParcRight)./2;
freesurferParcThickMean = (freesurferParcThickLeft + freesurferParcThickRight)./2;
freesurferParcDKTMean = (freesurferParcDKTLeft + freesurferParcDKTRight)./2;
freesurferParcThickDKTMean = (freesurferParcThickDKTLeft + freesurferParcThickDKTRight)./2;

indicesSubjectsWithFreesurfer = [];
for i = 1 : numel(freesurferSubjects)
    % Using volunteers table:
    indexThisSubjectInVolunteers = find(strcmp(freesurferSubjects{i}, volunteersTable.ID));
    groupsSubjectsWithFreesurfer(i) = volunteersTable.Grupo(indexThisSubjectInVolunteers);
    indicesSubjectsWithFreesurfer = [indicesSubjectsWithFreesurfer; indexThisSubjectInVolunteers];
    
end
% Check for outliers:
freesurferSubjects(isoutlier(freesurferBrainVolumesTable.BrainSegmentationVolume))
% Get variables and normalize:
totalBrainVolume = freesurferBrainVolumesTable.BrainSegmentationVolume;
volWithoutVent = freesurferBrainVolumesTable.BrainSegmentationVolumeWithoutVentricles;
subcorticalGMV = freesurferBrainVolumesTable.SubcorticalGrayMatterVolume;
corticalGMV = freesurferBrainVolumesTable.TotalCorticalGrayMatterVolume;
totalGMV = freesurferBrainVolumesTable.TotalGrayMatterVolume;
totalWMV = freesurferBrainVolumesTable.TotalCerebralWhiteMatterVolume;
totalVentVol = freesurferBrainVolumesTable.VolumeOfVentriclesAndChoroidPlexus;
eTIV_mm3 = freesurferBrainVolumesTable.EstimatedTotalIntracranialVolume;
meanGMThickness = freesurferBrainVolumesTable.MeanThickness;


% Statistical tests without normalization:
[pKW_surfVolWithoutVent, annovaTab_surfVolWithoutVent, stats_surfVolWithoutVent] = kruskalwallis(volWithoutVent, groupsSubjectsWithFreesurfer);
[pKW_surfSubcorticalGMV, annovaTab_surfSubcorticalGMV, stats_surfSubcorticalGMV] = kruskalwallis(subcorticalGMV, groupsSubjectsWithFreesurfer);
[pKW_surfCorticalGMV, annovaTab_surfCorticalGMV, stats_surfCorticalGMV] = kruskalwallis(corticalGMV, groupsSubjectsWithFreesurfer);
[pKW_surfTotalGMV, annovaTab_surfTotalGMV, stats_surfTotalGMV] = kruskalwallis(totalGMV, groupsSubjectsWithFreesurfer);
[pKW_surfTotalWMV, annovaTab_surfTotalWMV, stats_surfTotalWMV] = kruskalwallis(totalWMV, groupsSubjectsWithFreesurfer);
[pKW_surfTotalVentVol, annovaTab_surfTotalVentVol, stats_surfTotalVentVol] = kruskalwallis(totalVentVol, groupsSubjectsWithFreesurfer);
[pKW_surfeTIV, annovaTab_surfeTIV, stats_surfeTIV] = kruskalwallis(eTIV_mm3, groupsSubjectsWithFreesurfer);
[pKW_surfMeanGMThickness, annovaTab_surfMeanGMThickness, stats_surfMeanGMThickness] = kruskalwallis(meanGMThickness, groupsSubjectsWithFreesurfer);
close all

% Normalized by eTIV_mm3
normTBV = totalBrainVolume./eTIV_mm3;
normVolWithoutVent = volWithoutVent./eTIV_mm3;
normSubcorticalGMV = subcorticalGMV./eTIV_mm3;
normCorticalGMV = corticalGMV./eTIV_mm3;
normTotalGMV = totalGMV./eTIV_mm3;
normTotalWMV = totalWMV./eTIV_mm3;

[pKW_surfNormVolWithoutVent, annovaTab_surfNormVolWithoutVent, stats_surfNormVolWithoutVent] = kruskalwallis(normVolWithoutVent, groupsSubjectsWithFreesurfer);
[pKW_surfNormSubcorticalGMV, annovaTab_surfNormSubcorticalGMV, stats_surfNormSubcorticalGMV] = kruskalwallis(normSubcorticalGMV, groupsSubjectsWithFreesurfer);
[pKW_surfNormCorticalGMV, annovaTab_surfNormCorticalGMV, stats_surfNormCorticalGMV] = kruskalwallis(normCorticalGMV, groupsSubjectsWithFreesurfer);
[pKW_surfNormTotalGMV, annovaTab_surfNormTotalGMV, stats_surfNormTotalGMV] = kruskalwallis(normTotalGMV, groupsSubjectsWithFreesurfer);
[pKW_surfNormTotalWMV, annovaTab_surfNormTotalWMV, stats_surfNormTotalWMV] = kruskalwallis(normTotalWMV, groupsSubjectsWithFreesurfer);
[pKW_surfNormTBV, annovaTab_surfNormTBV, stats_surfNormTBV] = kruskalwallis(normTBV, groupsSubjectsWithFreesurfer);
close all
%% Statistical tests for all structures
indicesFreesurferCovid = strcmp(groupsSubjectsWithFreesurfer, 'Covid Prolongado');
indicesFreesurferControls = strcmp(groupsSubjectsWithFreesurfer, 'Control');
indicesFreesurferGroups = {indicesFreesurferCovid, indicesFreesurferControls};

% Brain volumes:
for i = 1 : size(freesurferBrainVolumes, 2)
    pKW_freesurferBrainVol(i) = kruskalwallis(freesurferBrainVolumes(:,i), groupsSubjectsWithFreesurfer, 'off');
    for j = 1 : numel(indicesFreesurferGroups)
        medians_freesurferBrainVol(i,j) = median(freesurferBrainVolumes(indicesFreesurferGroups{j},i));
        iqrInf_freesurferBrainVol(i,j) = prctile(freesurferBrainVolumes(indicesFreesurferGroups{j},i),25);
        iqrSup_freesurferBrainVol(i,j) = prctile(freesurferBrainVolumes(indicesFreesurferGroups{j},i),75);
    end
end
% Check the ones that are significative: 
roisDifferentBrainVol = find(pKW_freesurferBrainVol < 0.05);
% Show them in the table (i =2, first column are the names and the last one
% the groups):
disp('Regions with size differences from freesurfer summary results:')
freesurferBrainVolumesTable(1,roisDifferentBrainVol+1)
for i = 1 : numel(roisDifferentBrainVol)
    figure;
    set(gcf, 'Position', [100 100 800 600])
    boxchart(categorical(groupsSubjectsWithFreesurfer), freesurferBrainVolumesTable(:,roisDifferentBrainVol(i)));
    title(freesurferBrainVolumesTable.Properties.VariableNames(roisDifferentBrainVol(i)+1));
    saveas(gca, [outputPath 'DiffRegionsBrainVol_' freesurferBrainVolumesTable.Properties.VariableNames{roisDifferentBrainVol(i)+1}], 'png');
end

% Segmentation:
for i = 1 : size(freesurferSegmentation, 2) % 
    pKW_freesurferSegm(i) = kruskalwallis(freesurferSegmentation(:,i), groupsSubjectsWithFreesurfer, 'off');
    for j = 1 : numel(indicesFreesurferGroups)
        medians_freesurferSegm(i,j) = median(freesurferSegmentation(indicesFreesurferGroups{j},i));
        iqrInf_freesurferSegm(i,j) = prctile(freesurferSegmentation(indicesFreesurferGroups{j},i),25);
        iqrSup_freesurferSegm(i,j) = prctile(freesurferSegmentation(indicesFreesurferGroups{j},i),75);
    end
end
% Check the ones that are significative: 
roisDifferentSegm = find(pKW_freesurferSegm < 0.05);
% Show them in the table (i =2, first column are the names and the last one
% the groups):
disp('Regions with size differences from freesurfer segmentation:')
freesurferSegmentationTable(1,roisDifferentSegm+1)
for i = 1 : numel(roisDifferentSegm)
    figure;
    set(gcf, 'Position', [100 100 800 600])
    %boxchart(categorical(groupsSubjectsWithFreesurfer), table2array(freesurferSegmentationTable(:,roisDifferentSegm(i)+1)),...
    %    'GroupByColor', categorical(groupsSubjectsWithFreesurfer));
    boxchart(categorical(groupsSubjectsWithFreesurfer), freesurferSegmentation(:,roisDifferentSegm(i)));
    title(freesurferSegmentationTable.Properties.VariableNames(roisDifferentSegm(i)+1));
    saveas(gca, [outputPath 'DiffRegionsSegm_' freesurferSegmentationTable.Properties.VariableNames{roisDifferentSegm(i)+1}], 'png');
end

% Parcellation structures:
for i = 1 : size(freesurferParcMean, 2)
    pKW_freesurferParc(i) = kruskalwallis(freesurferParcMean(:,i), groupsSubjectsWithFreesurfer, 'off');
    for j = 1 : numel(indicesFreesurferGroups)
        medians_freesurferParc(i,j) = median(freesurferParcMean(indicesFreesurferGroups{j},i));
        iqrInf_freesurferParc(i,j) = prctile(freesurferParcMean(indicesFreesurferGroups{j},i),25);
        iqrSup_freesurferParc(i,j) = prctile(freesurferParcMean(indicesFreesurferGroups{j},i),75);
    end
end
% Check the ones that are significative: 
roisDifferentParc = find(pKW_freesurferParc < 0.05);
% Show them in the table (i =2, first column are the names and the last one
% the groups):
disp('Regions with size differences from freesurfer parcelation:')
freesurferParcLeftTable.Properties.VariableNames(roisDifferentParc+1)
for i = 1 : numel(roisDifferentParc)
    figure;
    set(gcf, 'Position', [100 100 800 600])
    boxchart(categorical(groupsSubjectsWithFreesurfer), freesurferParcMean(:,roisDifferentParc(i)));
    title(freesurferParcLeftTable.Properties.VariableNames(roisDifferentParc(i)+1));
    saveas(gca, [outputPath 'DiffRegionsParc_' freesurferParcLeftTable.Properties.VariableNames{roisDifferentParc(i)+1}], 'png');
end


% Parcellation thickness structures:
for i = 1 : size(freesurferParcThickMean, 2)
    pKW_freesurferParcThick(i) = kruskalwallis(freesurferParcThickMean(:,i), groupsSubjectsWithFreesurfer, 'off');
    for j = 1 : numel(indicesFreesurferGroups)
        medians_freesurferParcThick(i,j) = median(freesurferParcThickMean(indicesFreesurferGroups{j},i));
        iqrInf_freesurferParcThick(i,j) = prctile(freesurferParcThickMean(indicesFreesurferGroups{j},i),25);
        iqrSup_freesurferParcThick(i,j) = prctile(freesurferParcThickMean(indicesFreesurferGroups{j},i),75);
    end
end
% Check the ones that are significative: 
roisDifferentParcThick = find(pKW_freesurferParcThick < 0.05);
% Show them in the table (i =2, first column are the names and the last one
% the groups):
disp('Regions with size differences from freesurfer parcelation:')
freesurferParcLeftTable.Properties.VariableNames(roisDifferentParcThick+1)
for i = 1 : numel(roisDifferentParcThick)
    figure;
    set(gcf, 'Position', [100 100 800 600])
    boxchart(categorical(groupsSubjectsWithFreesurfer), freesurferParcThickMean(:,roisDifferentParcThick(i)));
    title(freesurferParcThickLeftTable.Properties.VariableNames(roisDifferentParcThick(i)+1));
    saveas(gca, [outputPath 'DiffRegionsParcThick_' freesurferParcThickLeftTable.Properties.VariableNames{roisDifferentParcThick(i)+1}], 'png');
end

% freesurferBrainVolumes(indicesToExclude,:) = [];
% freesurferSegmentation(indicesToExclude,:) = [];
% freesurferParc(indicesToExclude,:) = [];
% freesurferParcThick(indicesToExclude,:) = [];
% freesurferParcDKT(indicesToExclude,:) = [];
% freesurferParcThickDKT(indicesToExclude,:) = [];

%% PLOT BTV VS eTIV
figure;
subplot(1,3,1);
plot(volWithoutVent, eTIV_mm3, 'o');
xlabel('eTIV [mm^3]');
ylabel('Volume without Ventricles [mm^3]')
subplot(1,3,2);
plot(totalBrainVolume, eTIV_mm3, 'o');
xlabel('eTIV [mm^3]');
ylabel('Total Brain Volume [mm^3]')
subplot(1,3,3);
plot(totalGMV, eTIV_mm3, 'o');
xlabel('eTIV [mm^3]');
ylabel('Total Gray Matter Volume [mm^3]')
saveas(gca, [outputPath 'FreesurferCorrelationWitheTIV'], 'png');
%% COMPARE FREESURFER VS SIENAX
indicesSubjectsSienaxInFreesurfer = [];
% Match Subject Names:
for i = 1 : numel(sienaxResults.SubjectNames)
   indexThisSubjectInVolunteers = find(strcmp(sienaxResults.SubjectNames{i}, freesurferSubjects));
    if ~isempty(indexThisSubjectInVolunteers)
        indicesSubjectsSienaxInFreesurfer = [indicesSubjectsSienaxInFreesurfer; indexThisSubjectInVolunteers];
    else
        warning(sprintf('Freesurfer missing for %s.\n', sienaxResults.SubjectNames{i}));
    end

end
% Sienax in cm3
sienaxNormGMV = (sienaxResults.GreyMatterUnNormVolume./1000) ./ eTIV_mm3(indicesSubjectsSienaxInFreesurfer);
sienaxNormWMV = (sienaxResults.WhiteMatterUnNNormVolume./1000) ./ eTIV_mm3(indicesSubjectsSienaxInFreesurfer);
sienaxNormTBV = (sienaxResults.BrainUnNNormVolume./1000) ./ eTIV_mm3(indicesSubjectsSienaxInFreesurfer);

% Plot one against each other:
figure;
plot(sienaxNormTBV, normTBV(indicesSubjectsSienaxInFreesurfer),'x');
mdl = fitlm(sienaxNormTBV, normTBV(indicesSubjectsSienaxInFreesurfer), 'y ~ x1', 'Intercept', 0);
hold on;
x = [min([sienaxNormTBV; normTBV(indicesSubjectsSienaxInFreesurfer)]) : 0.05 : 1];
plot(x, x.*mdl.Coefficients.Estimate, 'LineWidth',2)
xlabel('Sienax');
ylabel('Freesurfer');
title('TBV/eTIV')

labelsMethod = {'Sienax', 'Freesurfer'};
lengthData = [numel(sienaxNormGMV); numel(normTotalGMV)];
allVolumeNormGMV = [sienaxNormGMV; normTotalGMV];
allVolumeNormWMV = [sienaxNormWMV; normTotalWMV];
allVolumeNormTBV = [sienaxNormTBV; normTBV];
labelsGroups = cell(numel(allVolumeNormGMV),1);
index = 1;
for i = 1 : numel(labelsMethod)
    labelsGroups(index:index+lengthData(i)-1) = labelsMethod(i);
    index = index + lengthData(i);
end
%boxchart(, normTotalGMV,'GroupByColor',tbl.Year)
figure;
set(gcf, 'Position', [100 100 2400 1000])
subplot(1,3,1);
boxchart(categorical(labelsGroups), allVolumeNormGMV);
title('Grey Matter Volume');
subplot(1,3,2);
boxchart(categorical(labelsGroups), allVolumeNormWMV);
title('White Matter Volume');
subplot(1,3,3);
boxchart(categorical(labelsGroups), allVolumeNormTBV);
title('Total Brain Volume');
saveas(gca, [outputPath 'SienaxVsFreesurferVolumeComparison'], 'png');
%% FREESURFER AND MULTIVARIATE ANALYSIS
j = 0;
k = 0;
% I need the questionnaire data for this:
indicesSubjectsWithFreesurferAndQ = [];
indicesCognitiveWithFreesurferAndQ = [];
for i = 1 : numel(freesurferSubjects)  
    % Using questionnaires:
    indexThisSubjectInVolunteers = find(strcmp(freesurferSubjects{i}, subjectNames));
    indexThisSubjectInCognitive = find(strcmp(freesurferSubjects{i}, subjectNamesTests));
    % Only leave subjects with all the data (Q and C)
    if ~isempty(indexThisSubjectInVolunteers) && ~isempty(indexThisSubjectInCognitive)
        j = j + 1;
        indicesVolunteersForFreesurferWithQ(j) = i;
        groupsSubjectsWithFreesurferAndQ(j) = group(indexThisSubjectInVolunteers);
        %indicesSubjectsWithMriReport(i) = indexThisSubjectInVolunteers;
        indicesSubjectsWithFreesurferAndQ = [indicesSubjectsWithFreesurferAndQ; indexThisSubjectInVolunteers];
        indicesCognitiveForFreesurfer(j) = i;
        indicesCognitiveWithFreesurferAndQ = [indicesCognitiveWithFreesurferAndQ; indexThisSubjectInCognitive];
    else
        warning(sprintf('Data from %s is missing in the volunteers or cognitive table.\n', freesurferSubjects{i}));
    end
end
% Remove subjects with missing data:
indicesToRemoveQ = find(isnan(age_years(indicesSubjectsWithFreesurferAndQ))>0); ...
indicesToRemoveCognitive = find(isnan(tmtA_score(indicesSienaxWithW))>0);


if ~isempty(indicesToRemoveQ)
    disp('Subjects to remove because of missing data:')
    subjectNames(indicesSubjectsWithFreesurferAndQ)
end
if ~isempty(indicesToRemoveCognitive)
    disp('Subjects to remove because of missing cognitive data:')
    subjectNamesTests(indicesToRemoveCognitive)
end
subjectNamesFreesurfer = subjectNames(indicesSubjectsWithFreesurferAndQ);
ageFreesurfer_years = age_years(indicesSubjectsWithFreesurferAndQ); 
groupCategorical = categorical( group(indicesSubjectsWithFreesurferAndQ));
groupBinary = zeros(size(groupCategorical)); groupBinary(groupCategorical == 'COVID') = 1;

tmtAFreesurfer_score = tmtA_score(indicesCognitiveWithFreesurferAndQ);
mocaFreesurfer_score = moca_score(indicesCognitiveWithFreesurferAndQ);
mocaFreesurfer_ref = moca_ref(indicesCognitiveWithFreesurferAndQ);

freesurferAllCognitiveTests = cognitiveTestsTable(indicesCognitiveWithFreesurferAndQ,3:29);
freesurferAllCognitiveTests_score = cognitiveTestsTable(indicesCognitiveWithFreesurferAndQ,3:3:29);
freesurferAllCognitiveTests_ref = cognitiveTestsTable(indicesCognitiveWithFreesurferAndQ,4:3:29);
freesurferAllCognitiveTests_perf = cognitiveTestsTable(indicesCognitiveWithFreesurferAndQ,5:3:29);

gmVolume = normTotalGMV(indicesCognitiveForFreesurfer);
wmVolume = normTotalWMV(indicesCognitiveForFreesurfer);
brainVolume = normTBV(indicesCognitiveForFreesurfer);

tbl = table(ageFreesurfer_years, groupCategorical, tmtAFreesurfer_score, ...
    mocaFreesurfer_score, mocaFreesurfer_ref, gmVolume); % Last column, response variable.
%% EXPLORE THE DATA
X = [table2array(freesurferAllCognitiveTests_score) gmVolume wmVolume brainVolume ageFreesurfer_years];
varNames = {'TMT_A'; 'TMT_B'; 'DSF'; 'DSB'; 'Stroop-W'; 'Stroop-C';...
    'Stroop-WC'; 'Stroop-interf'; 'MOCA'; 'GMV [cm^3]'; 'WMV [cm^3]'; 'TBV [cm^3]'; 'Age [years]'};
figure;
set(gcf, 'Position', [100 100 2600 2400]);
[H,AX,BigAx] = gplotmatrix(X,[],groupCategorical,['b' 'r' ],[],[20],true, 'grpbars', varNames);%'m' 'g' 'r'
%% CORRELATION MATRIX
X = [table2array(freesurferAllCognitiveTests_score) gmVolume wmVolume brainVolume ageFreesurfer_years groupBinary];
varNames = {'TMT_A'; 'TMT_B'; 'DSF'; 'DSB'; 'Stroop-W'; 'Stroop-C';...
    'Stroop-WC'; 'Stroop-interf'; 'MOCA'; 'GMV'; 'WMV'; 'TBV'; 'Age'; 'Long COVID'};
figure;
set(gcf, 'Position', [100 100 1200 1200]);
corrMatrix = corrcoef(X,'Rows','complete'); % To avoid NaNs in missing data.
% scatter plot
n = size(corrMatrix, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
% matrix with the correlation coeficients, but setting zeros in the lower
% diag, to have a full range of the colormap between [-1 1]
colorsMatrix = corrMatrix;
colorsMatrix(~triu(ones(size(corrMatrix)))) = -1;
scatter(x(:), y(:), 800.*abs(corrMatrix(:)), colorsMatrix(:), 'filled', 'MarkerFaceAlpha', 0.6)
% enclose markers in a grid
xl = [1:n+1;repmat(n+1, 1, n+1)];
xl = [xl(:, 1), xl(:, 1:end-1)];
yl = repmat(n+1:-1:1, 2, 1);
line(xl, yl, 'color', 'k') % horizontal lines
line(yl, xl, 'color', 'k') % vertical lines
% show labels
text([1:n] - 0.1, (n:-1:1) + 0.5, varNames  , 'HorizontalAlignment', 'right')
text((1:n) + 0.5, repmat(n + 1.2, n, 1) , varNames, ...
    'HorizontalAlignment', 'left', 'Rotation', 45)
h = gca;
hc = colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.9;
axis(h, 'equal')
colormap('jet')
set(hc, 'FontSize',16);
fullFilename = [outputPath 'CorrelationMatrixTestScoresAndImaging'];
saveas(gca, [fullFilename], 'tif');
%% CORRELATION MATRIX
X = [table2array(freesurferAllCognitiveTests_ref) gmVolume wmVolume brainVolume ageFreesurfer_years groupBinary];
varNames = {'TMT_A'; 'TMT_B'; 'DSF'; 'DSB'; 'Stroop-W'; 'Stroop-C';...
    'Stroop-WC'; 'Stroop-interf'; 'MOCA'; 'GMV'; 'WMV'; 'TBV'; 'Age'; 'Long COVID'};
figure;
set(gcf, 'Position', [100 100 1200 1200]);
corrMatrix = corrcoef(X,'Rows','complete'); % To avoid NaNs in missing data.
% scatter plot
n = size(corrMatrix, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
colorsMatrix = corrMatrix;
colorsMatrix(~triu(ones(size(corrMatrix)))) = -1;
scatter(x(:), y(:), 800.*abs(corrMatrix(:)), colorsMatrix(:), 'filled', 'MarkerFaceAlpha', 0.6)
% enclose markers in a grid
xl = [1:n+1;repmat(n+1, 1, n+1)];
xl = [xl(:, 1), xl(:, 1:end-1)];
yl = repmat(n+1:-1:1, 2, 1);
line(xl, yl, 'color', 'k') % horizontal lines
line(yl, xl, 'color', 'k') % vertical lines
% show labels
text([1:n] - 0.1, (n:-1:1) + 0.5, varNames  , 'HorizontalAlignment', 'right')
text((1:n) + 0.5, repmat(n + 1.2, n, 1) , varNames, ...
    'HorizontalAlignment', 'left', 'Rotation', 45)
h = gca;
hc = colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.9;
axis(h, 'equal')
colormap('jet')
set(hc, 'FontSize',16);
fullFilename = [outputPath 'CorrelationMatrixTestRefsAndImaging'];
saveas(gca, [fullFilename], 'tif');
close all
%% CORRELATION MORPHOMETRY VS SYMPTOMS
symptoms = table2array(summaryTable(:,colsLongCovidSymptoms));
symptomsFreesurfer = symptoms(indicesSubjectsWithFreesurferAndQ,:); 
symptomsNames = {'Headaches', 'Fatigue', 'Anosmia', 'Loss of Taste', ...
    'Dyspnoea', 'Muscle weakness', 'Muscle ache', 'Brain fog', 'Communication problems',...
    'Sleeping problems', 'Memory problems', 'Attention problems'};
labelsSymptoms = {'Yes, and they persist', 'Yes, but no anymore', 'No'};
% The symptoms are only for covid:
indicesFreesurferCovid = (groupCategorical == 'COVID');
subjectNamesFreesurferCovid = subjectNamesFreesurfer(indicesFreesurferCovid);
symptomsFreesurfer = symptomsFreesurfer(indicesFreesurferCovid,:);
%% CLUSTERING OF SYMPTOMS
maxClusters = 8;
for k = 1 : maxClusters
    [idc, C, sumDist] = kmeans(symptomsFreesurfer,k);
    meanDistPerCluster(k) = mean(sumDist);
end
figure; plot(meanDistPerCluster);
% between 3 and 4 is the optimal
%%
numClusters = 3;
[idc, C, sum] = kmeans(symptomsFreesurfer_years,numClusters);
range = abs(max(C)-min(C));
symptomsOrderedByRange = sort(range);
[symptomsOrderedByRange, indexSymptomsOrderedByRange] = sort(range, 'descend');
[histSymptoms, channels] = hist3([symptomsFreesurfer(:,indexSymptomsOrderedByRange(1)), symptomsFreesurfer(:,indexSymptomsOrderedByRange(2))]);
%figure; plot(symptomsFreesurfer(:,indexSymptomsOrderedByRange(1)), symptomsFreesurfer(:,indexSymptomsOrderedByRange(2)), 'o');
[X,Y] = meshgrid(channels{1}, channels{2});
indicesNonZeros = histSymptoms(:) ~=0;
scatter(X(indicesNonZeros), Y(indicesNonZeros), 800.*abs(histSymptoms(indicesNonZeros)), histSymptoms(indicesNonZeros), 'filled', 'MarkerFaceAlpha', 0.6);
%
markerStyles = {'o','s', 'd', 'x'};
figure; 
numSymptomsToShow = 4;
for s = 1 : numSymptomsToShow
    subplot(2,2,s);
    for i = 1 : numClusters
        indicesThisCluster = find((idc == i) > 0);
        plot(symptomsFreesurfer(indicesThisCluster,indexSymptomsOrderedByRange(2*s-1)) + 0.2*(rand(size(indicesThisCluster))-0.5), ...
            symptomsFreesurfer(indicesThisCluster,indexSymptomsOrderedByRange(2*s))  + 0.2*(rand(size(indicesThisCluster))-0.5), ...
            markerStyles{i}, 'LineWidth', 3, 'MarkerSize', 5);
        xlabel(symptomsNames(indexSymptomsOrderedByRange(2*s-1)));
        ylabel(symptomsNames(indexSymptomsOrderedByRange(2*s)));
        hold on;
    end
end
%% USING PCA TO ANALYZE SYMPTOMS
[components, scores, eigenevalues, tsquared, varExplained, mu] = pca(symptomsFreesurfer, 'Rows','pairwise');
figure;
subplot(1,3,1);
plot(components, 'LineWidth', 2);
subplot(1,3,2);
plot(eigenevalues, 'LineWidth', 2);
subplot(1,3,3);
plot(cumsum(varExplained), 'LineWidth', 2);

weightedComponents = components(:,1:4).*repmat(eigenevalues(1:4)',size(components, 1),1);
figure;
plot(weightedComponents, '-','LineWidth', 2);
legend(num2str([1:4]'));
%% CORRELATION SYMPTOMS VS COGNITIVE TESTS
freesurferCovidCognitiveTests_ref = freesurferAllCognitiveTests_ref(indicesFreesurferCovid,:);

X = [table2array(freesurferCovidCognitiveTests_ref) symptomsFreesurfer];
varNames = {'TMT_A'; 'TMT_B'; 'DSF'; 'DSB'; 'Stroop-W'; 'Stroop-C';...
    'Stroop-WC'; 'Stroop-interf'; 'MOCA'};
varNames = cat(1, varNames, symptomsNames');
figure;
set(gcf, 'Position', [100 100 1200 1200]);
corrMatrix = corrcoef(X,'Rows','complete'); % To avoid NaNs in missing data.
% scatter plot
n = size(corrMatrix, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
colorsMatrix = corrMatrix;
colorsMatrix(~triu(ones(size(corrMatrix)))) = -1;
scatter(x(:), y(:), 800.*abs(corrMatrix(:)), colorsMatrix(:), 'filled', 'MarkerFaceAlpha', 0.6)
% enclose markers in a grid
xl = [1:n+1;repmat(n+1, 1, n+1)];
xl = [xl(:, 1), xl(:, 1:end-1)];
yl = repmat(n+1:-1:1, 2, 1);
line(xl, yl, 'color', 'k') % horizontal lines
line(yl, xl, 'color', 'k') % vertical lines
% show labels
text([1:n] - 0.1, (n:-1:1) + 0.5, varNames  , 'HorizontalAlignment', 'right')
text((1:n) + 0.5, repmat(n + 1.2, n, 1) , varNames, ...
    'HorizontalAlignment', 'left', 'Rotation', 45)
h = gca;
hc = colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.9;
axis(h, 'equal')
colormap('jet')
set(hc, 'FontSize',16);
fullFilename = [outputPath 'CorrelationMatrixTestRefsAndSymptoms'];
saveas(gca, [fullFilename], 'tif');
%%
X = [table2array(freesurferAllCognitiveTests_ref) gmVolume wmVolume brainVolume ageFreesurfer_years groupBinary];
varNames = {'TMT_A'; 'TMT_B'; 'DSF'; 'DSB'; 'Stroop-W'; 'Stroop-C';...
    'Stroop-WC'; 'Stroop-interf'; 'MOCA'; 'GMV'; 'WMV'; 'TBV'; 'Age'; 'Long COVID'};
figure;
set(gcf, 'Position', [100 100 1200 1200]);
corrMatrix = corrcoef(X,'Rows','complete'); % To avoid NaNs in missing data.
% scatter plot
n = size(corrMatrix, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
colorsMatrix = corrMatrix;
colorsMatrix(~triu(ones(size(corrMatrix)))) = -1;
scatter(x(:), y(:), 800.*abs(corrMatrix(:)), colorsMatrix(:), 'filled', 'MarkerFaceAlpha', 0.6)
% enclose markers in a grid
xl = [1:n+1;repmat(n+1, 1, n+1)];
xl = [xl(:, 1), xl(:, 1:end-1)];
yl = repmat(n+1:-1:1, 2, 1);
line(xl, yl, 'color', 'k') % horizontal lines
line(yl, xl, 'color', 'k') % vertical lines
% show labels
text([1:n] - 0.1, (n:-1:1) + 0.5, varNames  , 'HorizontalAlignment', 'right')
text((1:n) + 0.5, repmat(n + 1.2, n, 1) , varNames, ...
    'HorizontalAlignment', 'left', 'Rotation', 45)
h = gca;
hc = colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.9;
axis(h, 'equal')
colormap('jet')
set(hc, 'FontSize',16);
fullFilename = [outputPath 'CorrelationMatrixTestRefsAndImaging'];
saveas(gca, [fullFilename], 'tif');
%%
mdlStep = stepwiselm(tbl, 'linear', 'Verbose',2);%'MPG ~ Weight'

mdl1 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical');
variables = mdl1.Coefficients.Variables;

mdl2 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical+mocaFreesurfer_score');
variables = mdl2.Coefficients.Variables;

mdl3 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical+mocaFreesurfer_ref');
variables = mdl3.Coefficients.Variables;
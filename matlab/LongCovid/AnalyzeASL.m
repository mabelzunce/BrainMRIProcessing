clear all
close all

dataPartitionPath = '/data/'; %'D:/'
imagingPartitionPath = '/data_imaging/'; %'F:/'
% Load data:
loadData = 0;

%% PATHS AND FILENAMES
studyDataPath = [dataPartitionPath '/UNSAM/CovidProject2/'];
resultsPath = [studyDataPath '/DataAnalysis/'];
responseExcelFilename = fullfile(studyDataPath, 'Respuestas.xlsx');
summaryExcelFilename = fullfile(studyDataPath, 'ResumenRespuestas.xlsx');
volunteersExcelFilename = fullfile(studyDataPath, 'VoluntariosProyectoCovidProlongado.xlsx');
% ASL
aslResultsPath = fullfile(resultsPath, 'ASL', 'Abstract');
cbfDKFilename = fullfile(aslResultsPath, 'mean_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_PVC2.csv');
cbfHammersFilename = fullfile(aslResultsPath, 'mean_qCBF_StandardSpace_Hammers_PVC2.csv');
cbfHOcortFilename = fullfile(aslResultsPath, 'mean_qCBF_StandardSpace_HOcort_CONN_2024_PVC2.csv');
cbfH0subFilename = fullfile(aslResultsPath, 'mean_qCBF_StandardSpace_HOsub_CONN_PVC2.csv');
expAslStructFilename = fullfile(aslResultsPath, 'mean_qCBF_StandardSpace_MNI_Structural_PVC2.csv');
expAslThalamusFilename = fullfile(aslResultsPath, 'mean_qCBF_StandardSpace_Thalamus_2024_PVC2.csv');
expAslTotalGMFilename = fullfile(aslResultsPath, 'mean_qCBF_StandardSpace_TotalGM_2024_PVC2.csv');
expAslTotalWMFilename = fullfile(aslResultsPath, 'mean_qCBF_StandardSpace_TotalWM_PVC2.csv');
expAslDeepWMFilename = fullfile(aslResultsPath, 'mean_qCBF_StandardSpace_DeepWM_2024_PVC2.csv');
% ASL Images
aslImagesPath = [imagingPartitionPath 'CovidProject/Estudio2/PreprocessedMRI/ASL/'];
aslAllImagesFilename = fullfile(aslImagesPath, 'asl_total_covid_4d.nii.gz');
aslAllImagesNamesFilename = fullfile(aslImagesPath, 'subjects_asl_total_covid_4d.csv');
% Structural
fslPreprocessedDataPath = [imagingPartitionPath '/CovidProject2/PreprocessedMRI/FSL/'];
freesurferBrainVolumesPath = [resultsPath '/AnatomicalMRI/freesurfer/'];
outputPath = [resultsPath 'BrainPerfusion/'];
parcelationVolFilename = 'parcellation_volumes_clean_all.csv';
parcelationThicknessFilename = 'parcellation_thickness_clean_all.csv';
segmentationFilename = 'segmentation_all.csv';
brainVolumes = 'brain_volumes_all.csv';
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
colsLongCovidSymptoms = 25:36;
colsHealthHistory = 38:55;
% % Individual:
% colId = 1; colGroup = 2; colAge = 4; colGender = 5; colHeight = 6; colWeight = 7; colVaccinated = 8; colNumVaccines = 9; 
% colDateLastVaccine = 10; colCovidInfection = 11; colNumCovidInfections = 12; colDateFirstInfection = 13; colNumVaccinesCovidInfection = 14;
% colDateLastSymptomsCovid = 15; colAdmitted = 16; colDateAdmitted = 17; colDateDischarged = 18; colAdmittedAfterCovid = 19; colNumAdmittedAfterCovid = 20;
% colAdmittedICU = 21; colFeelRecovered = 22; colDiagnosedWithLongCovid = 23;
% % Persisten symptoms post Covid:
% colHeadache = 24; colFatigue = 25; colAnosmia = 26; colLossTaste = 27; colDisnea = 28; colMuscleWeakness = 29; colMusclePain = 30;
% colFocusConfussion = 31; colTalkingComm = 32; colSleeping = 33; colMemory = 34; colAttention = 35; colOtherSymptoms = 36;
% % Health background:
% colHypertension = 37; colDiabetes = 38; colAsma = 39; colCholesterol = 40; colHeartAttack = 41;
% colAngina = 42; colThrombosis = 43; colBloodPressureMedication = 44; colAspirine = 45; colBloodThinner = 46;
% colOtherHealthIssues = 47; colOtherHealthIssuesDescr = 48; colSmoker = 49; colCigarettesPerDay = 50; colYearsSmoking = 51;
% colYearsExSmoker = 52; colDrinksAlcohol = 53; colAlcoholUnitsPerWeek = 54;

% FAS QUEATIONNAIRE
colFAS = 75;
%% SUBJECTS TO EXCLUDE:
excludeVascular = 1;
% Apply some filters for unwanted volunteers:
subjectsToExclude = {'CP0011', 'CP0015',''}; % 11 (foreign language) and 15 (elder and probaly with some form of dementia). Remove empty names.
% subjects to exclude because of artefacts or problems with the images:
subjectsToExcludeAsl = {'CP0015', 'CP0075', 'CP0076','CP0078', 'CP0079', 'CP0062', 'CP0105', 'CP0216'};
subjectsToExcludeAsl = {'CP0015', 'CP0075', 'CP0076','CP0078', 'CP0079', 'CP0062', 'CP0105', 'CP0216',...
    'CP0114', 'CP0140', 'CP0141', 'CP0144', 'CP0196', 'CP0219', 'CP0220', 'CP0221', 'CP0222'};% With missing data

subjectsWithVascularSginal = { 'CP0044', 'CP0053', 'CP0054', 'CP0061', 'CP0107', 'CP0108', 'CP0123', 'CP0142', ...
    'CP0147', 'CP0154', 'CP0176', 'CP0178', 'CP0183', 'CP0188', 'CP0193', 'CP0205',};% With missing data
subjectsWithVascularSginal = {'CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061', 'CP0096','CP0107','CP0108','CP0118','CP0123', 'CP0142', ...
    'CP0147', 'CP0154','CP0167', 'CP0176','CP0178', 'CP0183', 'CP0188','CP0193','CP0199', 'CP0205'};
if excludeVascular
    subjectsToExcludeAsl = [subjectsToExcludeAsl subjectsWithVascularSginal];
end
%% READ ASL DATA
cbfDK = readtable(cbfDKFilename);
cbfHammers = readtable(cbfHammersFilename);
cbfGM = readtable(expAslTotalGMFilename);
expAslStruct = readtable(expAslStructFilename);
expAslDeepWM = readtable(expAslDeepWMFilename);
expAslTotalWM = readtable(expAslTotalWMFilename);
% Sort by name the tables:
cbfDK = sortrows(cbfDK,1);
cbfHammers = sortrows(cbfHammers,1);
cbfGM = sortrows(cbfGM,1);
expAslStruct = sortrows(expAslStruct,1);
expAslDeepWM = sortrows(expAslDeepWM,1);
expAslTotalWM = sortrows(expAslTotalWM,1);

% Ger gourps of subjects with vascular:
for i = 1 : numel(subjectsWithVascularSginal)
    ind=find(strcmp(subjectsWithVascularSginal{i}, cbfDK.ID)>0);
    groupVascular{i} = cbfDK.Grupo{ind};
end

% Remove data to exclude:
indicesToExclude = [];
for i = 1 : numel(subjectsToExcludeAsl)
    ind=find(strcmp(subjectsToExcludeAsl{i}, cbfDK.ID)>0);
    if ~isempty(ind)
       indicesToExclude = [indicesToExclude ind];
    end
end
cbfDK(indicesToExclude,:) = [];
cbfHammers(indicesToExclude,:) = [];
expAslStruct(indicesToExclude,:) = [];
expAslDeepWM(indicesToExclude,:) = [];
expAslTotalWM(indicesToExclude,:) = [];
% FILTER COLUMNS WITH NANS
indicesNaNs = find(sum(isnan(table2array(cbfDK(:,2:end-1))),1) > 10); % Remove columns with more than 10 NaN.
cbfDK(:,indicesNaNs+1) = []; % +1 because the indices does not count the first column
indicesNaNs = find(sum(isnan(table2array(cbfHammers(:,2:end-1))),1) > 10); % Remove columns with more than 10 NaN.
cbfHammers(:,indicesNaNs+1) = []; % +1 because the indices does not count the first column

% Matrix with only the regions:
indexFirstRegion = 10;
hammersRegionNames = cbfHammers.Properties.VariableNames(indexFirstRegion:end-1);
cbfDK_regions = table2array(cbfDK(:,indexFirstRegion:end-1));
cbfHammers_regions = table2array(cbfHammers(:,indexFirstRegion:end-1));

% indicesCovid and controls:
group = cbfDK.Grupo;
%group{strcmp(group, 'DUDA')} = 'COVID';
indicesControls = find(strcmp(cbfDK.Grupo, groupNames{1}));
indicesCovid = find(strcmp(cbfDK.Grupo, groupNames{2}));
indicesControlsHammers = find(strcmp(cbfHammers.Grupo, groupNames{1}));
indicesCovidHammers = find(strcmp(cbfHammers.Grupo, groupNames{2}));

if sum(indicesControls ~= indicesControlsHammers) ~= 0
    error('Indices for controls dont match for different tables');
end
if sum(indicesCovid ~= indicesCovidHammers) ~= 0
    error('Indices for controls dont match for different tables');
end

%% LOAD QUESTIONNAIRE
%responsesTable = readtable(summaryExcelFilename, 'NumHeaderLines', 3);
responsesTable = readtable(responseExcelFilename,'Sheet', 'respuestasCuestionario', 'Range', 'A4:EJ193', 'ReadVariableNames', false, 'ReadRowNames', false);
summaryTable = readtable(summaryExcelFilename,'Sheet', 'resumenTotal','ReadVariableNames', true, 'ReadRowNames', false);

% Leave only the data of the subjects with ASL
indicesWithAsl = [];
indicesAslWithData = []; % Both lists need to be matched, there can be missing data in both.
indicesSummaryWithData = [];
for i = 1 : numel(responsesTable.Var1)
    ind=find(strcmp(responsesTable.Var1{i}, cbfDK.ID)>0);
    if ~isempty(ind)
       indicesWithAsl = [indicesWithAsl i];
       indicesAslWithData = [indicesAslWithData ind'];
    else
        warning(sprintf('Subject %s without ASL or demographics data.\n', responsesTable.Var1{i}));
        % responsesMatchedTable(i,1) = cbfDK.ID(i);
        % %responsesMatchedTable(i,2:end) = NaN;
        % for j = 2 : size(responsesMatchedTable,2)
        %     if isnumeric(table2array(responsesMatchedTable(i,j)))
        %         responsesMatchedTable(i,j) = array2table(NaN);
        %     end
        % end
    end

end

cbfDK = cbfDK(indicesAslWithData,:);
cbfHammers = cbfHammers(indicesAslWithData,:);
cbfHammers_regions = cbfHammers_regions(indicesAslWithData,:); %TODO: CHECK
cbfDK_regions = cbfDK_regions(indicesAslWithData,:);
responsesMatchedTable = responsesTable(indicesWithAsl,:);
for i = 1 : numel(responsesTable.Var1)
    ind=find(strcmp(responsesTable.Var1{i}, summaryTable.ID)>0);
    if ~isempty(ind)
       indicesSummaryWithData = [indicesSummaryWithData ind'];
    else
        warning(sprintf('Subject %s not in summary.\n', responsesTable.Var1{i}));
    end
end
summaryTable = summaryTable(indicesSummaryWithData);
%Re estimate indices for groups:
%group{strcmp(group, 'DUDA')} = 'COVID';
indicesControls = find(strcmp(cbfDK.Grupo, groupNames{1}));
indicesCovid = find(strcmp(cbfDK.Grupo, groupNames{2}));

% Get fields:
subjectNames = responsesMatchedTable.Var1;
group = responsesMatchedTable.Var2;
age_years = responsesMatchedTable.Var4;
gender = responsesMatchedTable.Var5;
height_cm = responsesMatchedTable.Var6;
weight_kg = responsesMatchedTable.Var7;
vaxCovid_si_no = responsesMatchedTable.Var8;
numVaxDoses = responsesMatchedTable.Var9;
covid_si_no = responsesMatchedTable.Var11;
numCovidInfections = responsesMatchedTable.Var12;
dateFirstInfection = responsesMatchedTable.Var13;
admittedHospForCovid = responsesMatchedTable.Var16;
dateAdmitHospForCovid = responsesMatchedTable.Var17;
dateDischargeHospForCovid = responsesMatchedTable.Var18;
readmittedHospAfterCovid = responsesMatchedTable.Var19;
icuAdmittedHospForCovid = responsesMatchedTable.Var21;
%daysSinceFirstEnfection = responsesMatchedTable.Var137; % Not reliables:
daysSinceFirstEnfection = days(responsesMatchedTable.Var3-responsesMatchedTable.Var13);
daysSinceLastEnfection = days(responsesMatchedTable.Var3-responsesMatchedTable.Var15);
numVaccinesWhenInfection = responsesMatchedTable.Var14;
unvaccinatedWhenInfection = responsesMatchedTable.Var14 == 0;
severityNames = {'No infectado', 'Paciente ambulatorio', 'Internación moderada', 'Internación severa (UCI)'};
for i = 1 : numel(responsesMatchedTable.Var22)
    if isempty(responsesMatchedTable.Var22{i})
        if strcmp(group{i}, 'CONTROL')
            responsesMatchedTable.Var22{i} = 'Paciente ambulatorio';
        else
            responsesMatchedTable.Var22{i} = 'No infectado';
        end
    end
end
severityCovid = categorical(responsesMatchedTable.Var22, severityNames);
% Symptoms:
symptom = table2array(responsesMatchedTable(:,colsLongCovidSymptoms));
symptomsNames = {'Headaches', 'Fatigue', 'Anosmia', 'Loss of Taste', ...
    'Dyspnoea', 'Muscle weakness', 'Muscle ache', 'Brain fog', 'Communication problems',...
    'Sleeping problems', 'Memory problems', 'Problems with attention'};

labelsSymptoms = {'Yes, and they persist', 'Yes, but no anymore', 'No'};
educationLevels = {'Primario incompleto', 'Primario completo','Secundario incompleto', ...
    'Secundario completo', 'Terciario incompleto', 'Universitario incompleto', 'Terciario completo', ...
    'Universitario completo', 'Posgrado', 'NaN'};
education = categorical(responsesMatchedTable.Var131, educationLevels);
% Health:
smoker = responsesMatchedTable.Var50;
smokerCategories = {'SI', 'EX FUMADOR', 'NO'};
smoker = categorical(smoker, smokerCategories);
alcoholConsumption = responsesMatchedTable.Var54;
alcoholConsumptionCategories = {'todos los días', 'los fines de semana', 'No'};
alcoholConsumption = categorical(alcoholConsumption, alcoholConsumptionCategories);
healthIssuesLabels = {'Presion alta', 'Diabetes', 'Asma', 'Colesterol alto', ...
    'Infarto cardiaco', 'Angina de pecho', 'Embolia o trombosis?', 'Medicacion por presion arterial', 'Aspirinas', 'Antiagregantes'};
healthIssues = table2array(responsesMatchedTable(:,38:47));

% Fatigue
fas = table2array(responsesMatchedTable(:,75));
% Eqol
eqolDimensiones = table2array(responsesMatchedTable(:,56:60));
eqolVAS = table2array(responsesMatchedTable(:,61));

%% SYMPTOMS
symptomCovid = symptom(strcmp(group, 'COVID'),:);
symptomCovid(isnan(symptomCovid)) = 3;
% Stats:
countsSymptoms = hist(symptomCovid, [1:3]);
percSymptoms = countsSymptoms./ size(symptomCovid,1)*100;
array2table(percSymptoms, 'VariableNames', symptomsNames, 'RowNames', {'Persistent', 'Not anymore', 'No'})
disp('Symptoms');
for i = 1 : size(percSymptoms,2)
    fprintf("\t %s: %d(%.1f%%) / %d(%.1f%%) / %d(%.1f%%)\n", symptomsNames{i},[countsSymptoms(:,i), percSymptoms(:,i)]');
end
%% DEMOGRAPHIC STATS
numTotalSubjects = numel(subjectNames);
groupNamesThisTable = unique(group);
for i = 1 : numel(groupNamesThisTable)
    indicesThisGroup = strcmp(group, groupNamesThisTable{i});
    numSubjectsGroup(i) = sum(indicesThisGroup)';
    sexFemale(i) = sum(strcmp(gender(indicesThisGroup),'F'));
    sexFemalePerc(i) = sexFemale(i)./numSubjectsGroup(i);
    age_mean(i) = mean(age_years(indicesThisGroup));
    age_std(i) = std(age_years(indicesThisGroup));
    meanWeight(i) = mean(weight_kg(indicesThisGroup));
    meanHeight(i) = mean(height_cm(indicesThisGroup));
    meanBmi(i) = mean(weight_kg(indicesThisGroup)./((height_cm(indicesThisGroup)/100).^2));
    stdBmi(i) = std(weight_kg(indicesThisGroup)./((height_cm(indicesThisGroup)/100).^2));
    % education:
    histEducation(:,i) = hist(education(indicesThisGroup), educationLevels)';
    histEducationPerc(:,i) = histEducation(:,i)./numSubjectsGroup(i);
    % Covid:
    covidInfection(i) = sum(strcmp(covid_si_no, 'SI'));
    covidInfectionPerc(i) = covidInfection(i)/numSubjectsGroup(i);
    meanNumCovidInfections(i) = mean(numCovidInfections(indicesThisGroup));
    stdNumCovidInfections(i) = std(numCovidInfections(indicesThisGroup));
    histSeverityCovid(:,i) = hist(severityCovid(indicesThisGroup))';
    histSeverityCovidPerc(:,i) = 100*histSeverityCovid(:,i)./numSubjectsGroup(i);
    meanDaysSinceFirstEnfection(i) = mean(daysSinceFirstEnfection(indicesThisGroup), "omitnan");
    stdDaysSinceFirstEnfection(i) = std(daysSinceFirstEnfection(indicesThisGroup), "omitnan");
    meanDaysSinceLastEnfection(i) = mean(daysSinceLastEnfection(indicesThisGroup), "omitnan");
    stdDaysSinceLastEnfection(i) = std(daysSinceLastEnfection(indicesThisGroup), "omitnan");
    numUnvaccinatedWhenInfection(i) = sum(unvaccinatedWhenInfection(indicesThisGroup));
    numUnvaccinatedFieldWithData(i) = sum(~isnan(numVaccinesWhenInfection(indicesThisGroup)));
    numUnvaccinatedWhenInfectionPerc(i) = 100*numUnvaccinatedWhenInfection(i)/numUnvaccinatedFieldWithData(i);
    % Health 
    histAlcohol(:,i) = hist(alcoholConsumption(indicesThisGroup), alcoholConsumptionCategories)';
    histAlcoholPerc(:,i) = 100*histAlcohol(:,i)./numSubjectsGroup(i);
    histSmoker(:,i) = hist(smoker(indicesThisGroup), smokerCategories)';
    histSmokerPerc(:,i) = 100*histSmoker(:,i)./numSubjectsGroup(i);
    numHealthIssues(:,i) = sum(strcmp(healthIssues(indicesThisGroup,:), 'SI'),1);
    numHealthIssuesPerc(:,i) = 100*numHealthIssues(:,i)./numSubjectsGroup(i);
end
% Education levels in only three categories: incomplete secundary or below,
% secundary or incomplete terciary, univeristary
for i = 1 : numel(groupNamesThisTable)
    histEducationThreeLevels(1,i) = sum(histEducation(1:3,i));
    histEducationThreeLevels(2,i) = sum(histEducation(4:6,i));
    histEducationThreeLevels(3,i) = sum(histEducation(7:end,i));
    histEducationThreeLevelsPerc(:,i) = 100*histEducationThreeLevels(:,i)./numSubjectsGroup(i);
end
for i = 1 : numel(groupNamesThisTable)
    disp(groupNamesThisTable{i});
    fprintf("Sex: %d (%.1f)\n", sexFemale(i), sexFemalePerc(i));
    fprintf("age: %.1f (%.1f)\n", age_mean(i), age_std(i));
    fprintf("bmi: %d (%.1f)\n", meanBmi(i), stdBmi(i));
    fprintf("age: %d (%.1f)\n", sexFemale(i), sexFemalePerc(i));
    fprintf("education: (%d,%d,%d) (%.1f,%.1f,%.1f)\n", histEducationThreeLevels(:,i), histEducationThreeLevelsPerc(:,i));
    fprintf("Severity covid: (%d,%d,%d,%d) (%.1f,%.1f,%.1f,%.1f)\n", histSeverityCovid(:,i), histSeverityCovidPerc(:,i));
    fprintf("days since first infection: %.1f (%.1f)\n", meanDaysSinceFirstEnfection(i), stdDaysSinceFirstEnfection(i));
    fprintf("days since first infection: %.1f (%.1f)\n", meanDaysSinceLastEnfection(i), stdDaysSinceLastEnfection(i));
    fprintf("unvaccinated when infected: %d/%d (%.1f)\n", numUnvaccinatedWhenInfection(i), numUnvaccinatedFieldWithData(i), numUnvaccinatedWhenInfectionPerc(i));
    % Health:
    fprintf("smoker: %d(%.1f%%) / %d(%.1f%%) / %d(%.1f%%)\n", [histSmoker(:,i), histSmokerPerc(:,i)]');
    fprintf("alcohol: %d(%.1f%%) / %d(%.1f%%) / %d(%.1f%%)\n", [histAlcohol(:,i), histAlcoholPerc(:,i)]');
    for j = 1 : size(healthIssues,2)
        fprintf("\t%s: %d(%.1f%%) / %d(%.1f%%) / %d(%.1f%%)\n", healthIssuesLabels{j}, [numHealthIssues(j,i), numHealthIssuesPerc(j,i)]');
    end
end


%% LOAD COGNTIVE DATA EXCEL
%cognitiveTestsTable = readtable(summaryExcelFilename,'Sheet','EvaluacionCognitiva', 'NumHeaderLines', 2);
%mriReportsTable = readtable(summaryExcelFilename,'Sheet','InformeResonancia', 'NumHeaderLines', 0);
cognitiveTable = readtable(summaryExcelFilename,'Sheet','EvaluacionCognitiva', 'Range', 'A3:AM191', 'ReadVariableNames', false, 'ReadRowNames', false);

% Leave only the data of the subjects with ASL
indicesAsl = [];
for i = 1 : numel(cbfDK.ID)
    ind=find(strcmp(cbfDK.ID{i}, cognitiveTable.Var1)>0);
    if ~isempty(ind)
       indicesAsl = [indicesAsl ind'];
       cognitiveMatchedTable(i,:) = cognitiveTable(ind, :);
    else
        warning(sprintf('Subject %s without data.\n', cbfDK.ID{i}));
        cognitiveMatchedTable(i,1) = cbfDK.ID(i);
        for j = 2 : size(cognitiveMatchedTable,2)
            if isnumeric(table2array(cognitiveMatchedTable(i,j)))
                cognitiveMatchedTable(i,j) = array2table(NaN);
            end
        end
    end
end

cognitiveTable = cognitiveTable(indicesAsl,:);

% SUMMARY IN A TABLE:
allCognitiveTests = cognitiveMatchedTable(:,3:29);
allCognitiveTests_score = table2array(cognitiveMatchedTable(:,3:3:29));
allCognitiveTests_perc = table2array(cognitiveMatchedTable(:,4:3:29));
allCognitiveTests_perf = table2array(cognitiveMatchedTable(:,5:3:29));

%% ANCOVA
moca_perc = allCognitiveTests_perc(:,end);
tbl = table(group, age_years, gender, moca_perc, fas);
for i = 1 : size(cbfHammers_regions, 2)
    metric = cbfHammers_regions(:,i);    
    %[p(i), table{i}, stats{i}] = anova(tbl, "metric ~ group + age + eTIV");
    anovaTable{i} = anova(tbl, metric,categoricalFactors=["group", "gender"]);
    p(i) = anovaTable{i}.stats.pValue(1);
    p_age(i) = anovaTable{i}.stats.pValue(2);
    p_sex(i) = anovaTable{i}.stats.pValue(3);
    p_moca(i) = anovaTable{i}.stats.pValue(4);
    p_fas(i) = anovaTable{i}.stats.pValue(5);
end

% Identify regions with signficant p values
indicesSign_cbfHammers = find(p < 0.05);
% Get the names:
disp('Regions with signifcant differences:')


for i = 1 : numel(indicesSign_cbfHammers)
    fprintf("Hammers. Region %s, p=%.3f, p_age=%.3f, p_sex=%.3f, p_moca=%.3f, p_fas=%.3f. Median cbf: controls=%.1f, covid=%.1f. Num NaN: %d\n", ...
        cbfHammers.Properties.VariableNames{indexFirstRegion+indicesSign_cbfHammers(i)-1}, p(indicesSign_cbfHammers(i)), p_age(indicesSign_cbfHammers(i)), p_sex(indicesSign_cbfHammers(i)), p_moca(indicesSign_cbfHammers(i)), p_fas(indicesSign_cbfHammers(i)), ...
        median(cbfHammers_regions(indicesControls,indicesSign_cbfHammers(i)), "omitnan"), median(cbfHammers_regions(indicesCovid, indicesSign_cbfHammers(i)), "omitnan"),...
        sum(isnan(cbfHammers_regions(:, indicesSign_cbfHammers(i)))));
end
%% GET A TABLE WITH THESE RESULTS
tbl = table(cbfDK.ID, group, age_years, gender, moca_perc, fas);

for i = 1 : numel(indicesSign_cbfHammers)
    tbl(:,end+1) = array2table(cbfHammers_regions(:,indicesSign_cbfHammers(i)));
    tbl.Properties.VariableNames{end}= cbfHammers.Properties.VariableNames{indexFirstRegion+indicesSign_cbfHammers(i)-1};
end
writetable(tbl, fullfile(outputPath, 'AncovaSignificant.csv'))

%% TABLE TO EXPLORE MACHINE LEARNING METHODS
tbl = table(group, age_years, gender, moca_perc, fas);
%allCognitiveTests_perc);
%% READ 4D IMAGE
imageAllAsl = niftiread(aslAllImagesFilename);
imageInfo = niftiinfo(aslAllImagesFilename);
% CREATE SMOOTHED IMAGES
fwhm_mm = [4 4 4]; sigma_mm = fwhm_mm./2.35; sigma_voxels = sigma_mm./imageInfo.PixelDimensions(1:3);
fwhm8_mm = [8 8 8]; sigma8_mm = fwhm8_mm./2.35; sigma8_voxels = sigma8_mm./imageInfo.PixelDimensions(1:3);
for i = 1 : size(imageAllAsl,5)
    imageAllAslSmoothed_4mm(:,:,:,i) = imgaussfilt3(imageAllAsl(:,:,:,1,i), sigma_voxels);
    imageAllAslSmoothed_8mm(:,:,:,i) = imgaussfilt3(imageAllAsl(:,:,:,1,i), sigma8_voxels);
end
imageInfo4D = imageInfo;
imageInfo4D.ImageSize = size(imageAllAslSmoothed_4mm);
imageInfo4D.PixelDimensions = [imageInfo.PixelDimensions(1:3) 1];
imageInfo4D.raw.dim = [4 imageInfo4D.ImageSize];
imageInfo4D.raw.pixdim = [-1 imageInfo4D.PixelDimensions];
niftiwrite(imageAllAslSmoothed_4mm, fullfile(aslImagesPath, 'asl_all_subjects_smoothed_4mm.nii'), imageInfo4D, 'Compressed', 1);
niftiwrite(imageAllAslSmoothed_8mm, fullfile(aslImagesPath, 'asl_all_subjects_smoothed_8mm.nii'), imageInfo4D, 'Compressed', 1);

indicesImages = readtable(aslAllImagesNamesFilename);
indicesImagesToInclude = [];
for i = 1 : numel(cbfDK.ID)
    ind=find(strcmp(cbfDK.ID{i}, indicesImages.Subject)>0);
    if ~isempty(ind)
       indicesImagesToInclude = [indicesImagesToInclude ind'];
    else
        warning(sprintf('Subject %s without image.\n', cbfDK.ID{i}));
    end
end
imageAllAsl = imageAllAsl(:,:,:,indicesImagesToInclude);
imageAllAslSmoothed_4mm = imageAllAslSmoothed_4mm(:,:,:,indicesImagesToInclude);
imageAllAslSmoothed_8mm = imageAllAslSmoothed_8mm(:,:,:,indicesImagesToInclude);
indicesImages = indicesImages(indicesImagesToInclude,:);

% Get an atlas for each group:
subjectsControls = cbfDK.ID(indicesControls);
subjectsCovid = cbfDK.ID(indicesCovid);
imagesAslControls = imageAllAslSmoothed_4mm(:,:,:,indicesControls);
imagesAslCovid = imageAllAslSmoothed_4mm(:,:,:,indicesCovid);
meanCbfImageControls = mean(imagesAslControls,4);
meanCbfImageCovid = mean(imagesAslCovid,4);
stdCbfImageControls = std(imagesAslControls,0,4);
stdCbfImageCovid = std(imagesAslCovid,0,4);

medianCbfImageControls = median(imagesAslControls,4);
medianCbfImageCovid = median(imagesAslCovid,4);

imageInfo.ImageSize = size(meanCbfImageControls);
imageInfo.PixelDimensions = imageInfo.PixelDimensions(1:3);
imageInfo.raw.dim = [3 imageInfo.ImageSize];
imageInfo.raw.pixdim = [-1 imageInfo.PixelDimensions];
niftiwrite(meanCbfImageControls, fullfile(aslImagesPath, 'meanCbfImageControls.nii'), imageInfo, 'Compressed', 1);
niftiwrite(meanCbfImageCovid, fullfile(aslImagesPath, 'meanCbfImageCovid.nii'), imageInfo, 'Compressed', 1);
niftiwrite(stdCbfImageControls, fullfile(aslImagesPath, 'stdCbfImageControls.nii'), imageInfo, 'Compressed', 1);
niftiwrite(stdCbfImageCovid, fullfile(aslImagesPath, 'stdCbfImageCovid.nii'), imageInfo, 'Compressed', 1);
niftiwrite(medianCbfImageControls, fullfile(aslImagesPath, 'medianCbfImageControls.nii'), imageInfo, 'Compressed', 1);
niftiwrite(medianCbfImageCovid, fullfile(aslImagesPath, 'medianCbfImageCovid.nii'), imageInfo, 'Compressed', 1);
%% BOXPLOTS FOR HAMMERS
% Concatenate data from both groups
cbf_all = [meanRoisHammersControls; meanRoisHammersCovid];  % Matrix ((N1 + N2) x M)

% Create group labels
group_labels = [repmat({'CONTROLS'}, size(meanRoisHammersControls, 1), 1); repmat({'COVID'}, size(meanRoisHammersCovid, 1), 1)];

% Create ROI labels
num_rois = size(cbf_all, 2);
roi_labels = repmat(1:num_rois, size(cbf_all, 1), 1);
roi_labels = roi_labels(:);  % Vectorize for boxplot

% Reshape data and group labels to match the boxplot input format
%cbf_all = cbf_all(:);  % Convert matrix to column vector
%group_labels = repmat(group_labels, num_rois, 1);  % Repeat group labels for each ROI

figure; subplot(1,3,1);boxchart(categorical(group_labels), cbf_all(:,85))
subplot(1,3,2);boxchart(categorical(group_labels), cbf_all(:,12))
subplot(1,3,3);boxchart(categorical(group_labels), cbf_all(:,33))
%% LOAD ATLASES
mni152 = niftiread('/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz');
mni152_info = niftiread('/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz');
mni152_brain_mask = niftiread('/usr/local/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii.gz');
vba_mask = niftiread(fullfile(aslImagesPath,'VBA_mask_final.nii.gz'));
mni152_gm = analyze75read('/usr/local/fsl/data/standard/tissuepriors/avg152T1_gray.hdr')>100; % Create binary mask
mni152_gm = uint8(permute(mni152_gm, [2 1 3]));
% Hammers atlas:
hammers = niftiread(fullfile(aslImagesPath,'Hammers_mith-n30r95-MaxProbMap-gm-MNI152-SPM12_resampled.nii.gz'));
mni152_gm = imresize3(mni152_gm, size(meanCbfImageControls), 'nearest');
hammers_gm = uint8(hammers).* mni152_gm;
mni152_brain_mask_resampled = single(imresize3(mni152_brain_mask, size(meanCbfImageControls)));
niftiwrite(mni152_brain_mask_resampled, fullfile(aslImagesPath, 'mni152_resampled.nii'), imageInfo, 'Compressed', 1);
imageInfoUint8 = imageInfo;
imageInfoUint8.Datatype = 'uint8';
niftiwrite(uint8(mni152_gm), fullfile(aslImagesPath, 'mni152_gm_resampled.nii'), imageInfoUint8, 'Compressed', 1);

%% GET IMAGES FOR HAMMERS ATLASES
meanHammersCovid = zeros(size(hammers_gm));
meanHammersControls = zeros(size(hammers_gm));
labels = unique(hammers_gm);
for i = 1 : numel(labels)
    label = labels(i);
    voxelsRoi = find(hammers_gm==label);
    for j = 1 :size(imagesAslControls, 4)
        thisImage = imagesAslControls(:,:,:,j);
        meanRoisHammersControls(j,i) = mean(thisImage(voxelsRoi));
    end
    meanHammersControls(voxelsRoi) = mean(meanRoisHammersControls(:,i));
    for j = 1 :size(imagesAslCovid, 4)
        thisImage = imagesAslCovid(:,:,:,j);
        meanRoisHammersCovid(j,i) = mean(thisImage(voxelsRoi));
    end
    meanHammersCovid(voxelsRoi) = mean(meanRoisHammersCovid(:,i));
end
tableControls = table();
tableControls.ID = subjectsControls;
for i = 1 : numel(labels)
    tableControls(:,i+1) = array2table(meanRoisHammersControls(:,i));
    tableControls.Properties.VariableNames(i+1) = sprintf("Labels %d", labels(i));
end
tableCovid = table();
tableCovid.ID = subjectsCovid;
for i = 1 : numel(labels)
    tableCovid(:,i+1) = array2table(meanRoisHammersCovid(:,i));
    tableCovid.Properties.VariableNames(i+1) = sprintf("Labels %d", labels(i));
end

niftiwrite(single(meanHammersControls), fullfile(aslImagesPath, 'meanHammersControls.nii'), imageInfo, 'Compressed', 1);
niftiwrite(single(meanHammersCovid), fullfile(aslImagesPath, 'meanHammersCovid.nii'), imageInfo, 'Compressed', 1);
writetable(tableControls, fullfile(aslImagesPath, 'meanRoisHammersControls.csv'))
writetable(tableCovid, fullfile(aslImagesPath, 'meanRoisHammersCovid.csv'))

%% PROCESS INDIVIDUAL IMAGES FROM EXPLORE ASL USING MASK
pathPopulationAnalysis = '/media/martin/PortableSSD/ASLDirSol/BIDS/derivatives/ExploreASL/Population/';

for j = 1 : numel(cbfHammers.ID)
    cbfMasked = niftiread(fullfile(pathPopulationAnalysis, sprintf('qCBF_masked_sub-%s_1_ASL_1.nii.gz', cbfHammers.ID{j})));
    gmMask = niftiread(fullfile(pathPopulationAnalysis, sprintf('PV_pGM_sub-%s_1.nii.gz', cbfHammers.ID{j})))>0.5;
    for i = 1 : numel(labels)
        label = labels(i);
        voxelsRoi = find(hammers_gm==label & gmMask);
        meanRoisHammersMaskedImages(j,i) = mean(cbfMasked(voxelsRoi),  "omitnan");
    end
    j
end

%% ANCOVA HAMMERS IMAGES
moca_perc = allCognitiveTests_perc(:,end);
tbl = table(group, age_years, gender, moca_perc, fas);
for i = 1 : size(meanRoisHammersMaskedImages, 2)
    metric = meanRoisHammersMaskedImages(:,i);    
    %[p(i), table{i}, stats{i}] = anova(tbl, "metric ~ group + age + eTIV");
    anovaTable{i} = anova(tbl, metric,categoricalFactors=["group", "gender"]);
    p(i) = anovaTable{i}.stats.pValue(1);
    p_age(i) = anovaTable{i}.stats.pValue(2);
    p_sex(i) = anovaTable{i}.stats.pValue(3);
    p_moca(i) = anovaTable{i}.stats.pValue(4);
    p_fas(i) = anovaTable{i}.stats.pValue(5);
end

% Identify regions with signficant p values
indicesSign_cbfHammers = find(p < 0.05);
% Get the names:
disp('Regions with signifcant differences:')


for i = 1 : numel(indicesSign_cbfHammers)
    fprintf("Hammers. Label %d, p=%.3f, p_age=%.3f, p_sex=%.3f, p_moca=%.3f, p_fas=%.3f. Median cbf: controls=%.1f, covid=%.1f. Num NaN: %d\n", ...
        indicesSign_cbfHammers(i), p(indicesSign_cbfHammers(i)), p_age(indicesSign_cbfHammers(i)), p_sex(indicesSign_cbfHammers(i)), p_moca(indicesSign_cbfHammers(i)), p_fas(indicesSign_cbfHammers(i)), ...
        median(cbfHammers_regions(indicesControls,indicesSign_cbfHammers(i)), "omitnan"), median(cbfHammers_regions(indicesCovid, indicesSign_cbfHammers(i)), "omitnan"),...
        sum(isnan(cbfHammers_regions(:, indicesSign_cbfHammers(i)))));
end
%% COMPARE HAMMERS FROM IMAGES TO ROIS
% Labels in hammers
labelMiddleInferiotTemporalLobeR = 13;
labelAngularGyrusR = 33;
labelSuperiorParietalGyrusR = 63;
indicesHammersImages = [labelMiddleInferiotTemporalLobeR labelAngularGyrusR labelSuperiorParietalGyrusR];
% From ExploreASL table:
indicesHammersExploreAsl = [20 46 76];
figure;
for i = 1 : numel(indicesHammersExploreAsl)
    subplot(1, numel(indicesHammersExploreAsl), i);
    plot([cbfHammers_regions(indicesControls,indicesHammersExploreAsl(i)); cbfHammers_regions(indicesCovid,indicesHammersExploreAsl(i))],...
        [meanRoisHammersControls(:, indicesHammersImages(i));meanRoisHammersCovid(:, indicesHammersImages(i))],'o')
    xlabel('CBF Rois Explore ASL')
    ylabel('CBF Rois from Images')
end
figure;
for i = 1 : numel(indicesHammersExploreAsl)
    subplot(1, numel(indicesHammersExploreAsl), i);
    plot([cbfHammers_regions(indicesControls,indicesHammersExploreAsl(i)); cbfHammers_regions(indicesCovid,indicesHammersExploreAsl(i))],...
        [meanRoisHammersMaskedImages(indicesControls,indicesHammersImages(i)); meanRoisHammersMaskedImages(indicesCovid,indicesHammersImages(i))],'o')
    xlabel('CBF Rois Explore ASL')
    ylabel('CBF Rois from Masked Images')
end

figure;
for i = 1 : numel(indicesHammersExploreAsl)
    subplot(1, numel(indicesHammersExploreAsl), i);
    plot([meanRoisHammersControls(:, indicesHammersImages(i));meanRoisHammersCovid(:, indicesHammersImages(i))],...
        [meanRoisHammersMaskedImages(indicesControls,indicesHammersImages(i)); meanRoisHammersMaskedImages(indicesCovid,indicesHammersImages(i))],'o')
    xlabel('CBF Rois Explore ASL')
    ylabel('CBF Rois from Masked Images')
end
%% 
% figure;
% boxplot(cbf_all, {roi_labels, group_labels}, 'FactorSeparator', 1, ...
%     'Labels', repmat(string(1:num_rois), 1, 2), 'LabelOrientation', 'inline');
% title('CBF per ROI Comparison between Groups');
% xlabel('ROI');
% ylabel('CBF');

% %% PERFORM T-TEST AT VOXEL LEVEL
% % [h_voxels, p_voxels, ci_voxels, stats_voxels] = ttest2(permute(imagesAslControls, [4 1 2 3]), ...
% %     permute(imagesAslCovid, [4 1 2 3]), 'tail', 'right');
% % Hipoperfusion
% [h_voxels, p_voxels, ci_voxels, stats_voxels] = ttest2(permute(imagesAslControls.*repmat(mni152_brain_mask_resampled,[1 1 1 size(imagesAslControls,4)]), [4 1 2 3]), ...
%     permute(imagesAslCovid.*repmat(mni152_brain_mask_resampled,[1 1 1 size(imagesAslCovid,4)]), [4 1 2 3]), 'tail', 'right', 'alpha', 0.005);
% h_voxels = permute(h_voxels, [2 3 4 1]);
% p_voxels = permute(p_voxels, [2 3 4 1]);
% tstats_orig = permute(stats_voxels.tstat, [2 3 4 1]);
% tstats_orig(isnan(tstats_orig)) = 0;
% niftiwrite(h_voxels, fullfile(aslImagesPath, 'h_voxels.nii'), imageInfo, 'Compressed', 1);
% niftiwrite(p_voxels, fullfile(aslImagesPath, 'p_voxels.nii'), imageInfo, 'Compressed', 1);
% 
% % Hyperperfusion
% [h_hyper_voxels, p_hyper_voxels, ci_hyper_voxels, stats_hyper_voxels] = ttest2(permute(imagesAslControls.*repmat(mni152_brain_mask_resampled,[1 1 1 size(imagesAslControls,4)]), [4 1 2 3]), ...
%     permute(imagesAslCovid.*repmat(mni152_brain_mask_resampled,[1 1 1 size(imagesAslCovid,4)]), [4 1 2 3]), 'tail', 'left', 'alpha', 0.005);
% h_hyper_voxels = permute(h_hyper_voxels, [2 3 4 1]);
% p_hyper_voxels = permute(p_hyper_voxels, [2 3 4 1]);
% niftiwrite(h_hyper_voxels, fullfile(aslImagesPath, 'h_hyper_voxels.nii'), imageInfo, 'Compressed', 1);
% niftiwrite(p_hyper_voxels, fullfile(aslImagesPath, 'p_hyper_voxels.nii'), imageInfo, 'Compressed', 1);

%% PERFORM ANCOVA AT VOXEL LEVEL
sizeImage = size(imageAllAslSmoothed_8mm);
sizeImage = sizeImage(1:3);
tbl = table(group, age_years, gender);
% Anova doesn't work at voxel level:
indicesBrain = find(mni152_brain_mask_resampled > 0);
indicesBrain = find(vba_mask > 0); % Use vba mask
p_anova = zeros(sizeImage, 'single');
p_age = zeros(sizeImage, 'single');
p_gender = zeros(sizeImage, 'single');
diff_means = zeros(sizeImage, 'single');
meanControls = zeros(sizeImage, 'single');
meanCovid = zeros(sizeImage, 'single');
h_anova = zeros(sizeImage, 'logical');
h_hypo_anova = zeros(sizeImage, 'logical');
h_hyper_anova = zeros(sizeImage, 'logical');
vecImage = reshape(imageAllAslSmoothed_8mm, [], size(imageAllAslSmoothed_8mm,4));
for i = 1 :numel(indicesBrain)
    anovaTable{i} = anova(tbl, vecImage(indicesBrain(i),:));
    p_anova(indicesBrain(i)) = anovaTable{i}.stats.pValue(1);
    p_age(indicesBrain(i)) = anovaTable{i}.stats.pValue(2);
    p_gender(indicesBrain(i)) = anovaTable{i}.stats.pValue(3);
    % Hipoperfusion
    h_anova(indicesBrain(i)) = anovaTable{i}.stats.pValue(1) < 0.005;
    if anovaTable{i}.stats.pValue(1) < 0.005
        h_hypo_anova(indicesBrain(i)) = (anovaTable{i}.groupmeans('group').Mean(1)>anovaTable{i}.groupmeans('group').Mean(2));
        h_hyper_anova(indicesBrain(i)) = (anovaTable{i}.groupmeans('group').Mean(1)<anovaTable{i}.groupmeans('group').Mean(2));
        meanControls(indicesBrain(i)) = anovaTable{i}.groupmeans('group').Mean(1);
        meanCovid(indicesBrain(i)) = anovaTable{i}.groupmeans('group').Mean(2);
    end
    i
end
diff_means = meanControls - meanCovid;
niftiwrite(single(h_anova), fullfile(aslImagesPath, 'h_anova_voxels.nii'), imageInfo, 'Compressed', 1);
niftiwrite(single(h_hypo_anova), fullfile(aslImagesPath, 'h_anova_hypo_voxels.nii'), imageInfo, 'Compressed', 1);
niftiwrite(single(h_hyper_anova), fullfile(aslImagesPath, 'h_anova_hyper_voxels.nii'), imageInfo, 'Compressed', 1);
niftiwrite(single(diff_means), fullfile(aslImagesPath, 'diff_control-covid.nii'), imageInfo, 'Compressed', 1);
niftiwrite(single(meanControls), fullfile(aslImagesPath, 'meanControls_sign.nii'), imageInfo, 'Compressed', 1);
niftiwrite(single(meanCovid), fullfile(aslImagesPath, 'meanCovid_sign.nii'), imageInfo, 'Compressed', 1);
niftiwrite(single(p_anova), fullfile(aslImagesPath, 'p_anova_voxels.nii'), imageInfo, 'Compressed', 1);
niftiwrite(single(p_age), fullfile(aslImagesPath, 'p_age_voxels.nii'), imageInfo, 'Compressed', 1);
niftiwrite(single(p_gender), fullfile(aslImagesPath, 'p_gender_voxels.nii'), imageInfo, 'Compressed', 1);
% Boxlplots regions:
%h_hyper_anova

%% MULTIPLE COMPARISONS
numPermutations = 500;
fprintf('Permutations:\n');
for j = 1 : numPermutations
    reshuffled_groups = group(randperm(numel(group)));
    indicesCOVIDShuffled = strcmp(reshuffled_groups, 'Covid Prolongado');
    indicesControlShuffled = ~indicesCOVIDShuffled;
    [h_aux, p_aux, ci, stats] = ttest2(permute(imageAllAslSmoothed_8mm(:,:,:,indicesControlShuffled), [4 1 2 3]), ...
        permute(imageAllAslSmoothed_8mm(:,:,:,indicesCOVIDShuffled), [4 1 2 3]), 'tail', 'right');
    p_aux = permute(p_aux, [2 3 4 1]);
    tstats_aux = permute(stats.tstat, [2 3 4 1]);
    tstats(:,:,:,j) = tstats_aux;
    pstats(:,:,:,j) = p_aux;
    fprintf('%d\n', j);
end
%%
tstats_max_perm = max(reshape(tstats,[], numPermutations) ,[],1); % max tstats for each permutation
tstats_min_perm = min(reshape(tstats,[], numPermutations),[],1); % min tstats for each permutation
% Get number of values lower than the tstats_max_perm for each of the
% original values:
num_greater = zeros(size(tstats_orig));
num_lower = zeros(size(tstats_orig));
for i = 1 : numPermutations
    num_greater = num_greater + (tstats_max_perm(i) > tstats_orig);
    num_lower = num_lower + (tstats_min_perm(i) < tstats_orig);
end
p_greater = single(num_greater./ numPermutations);
p_lower = single(num_lower ./ numPermutations);

niftiwrite(p_greater, fullfile(aslImagesPath, 'p_greater_perm_voxels.nii'), imageInfo, 'Compressed', 1);
niftiwrite(p_lower, fullfile(aslImagesPath, 'p_lower_perm_voxels.nii'), imageInfo, 'Compressed', 1);
%% SIMPLE KRUSKAL WALLIS TESTS
for i = 1 : size(cbfDK_regions, 2)
    p_cbfDK(i) = kruskalwallis(cbfDK_regions(:,i), cbfDK.Grupo, 'off');
end
for i = 1 : size(cbfHammers_regions, 2)
    p_cbfHammers(i) = kruskalwallis(cbfHammers_regions(:,i), cbfHammers.Grupo, 'off');
end
% Identify regions with signficant p values
indicesSign_cbfDK = find(p_cbfDK < 0.05);
indicesSign_cbfHammers = find(p_cbfHammers < 0.05);
% Get the names:
disp('Regions with signifcant differences:')
for i = 1 : numel(indicesSign_cbfDK)
    fprintf("DK. Region %s. Medians cbf controls %.1f, covid %.1f. Num NaN: %d\n", cbfDK.Properties.VariableNames{indexFirstRegion+indicesSign_cbfDK(i)-1}, ...
        median(cbfDK_regions(indicesControls,indicesSign_cbfDK(i)), "omitnan"), median(cbfDK_regions(indicesCovid, indicesSign_cbfDK(i)), "omitnan"),...
        sum(isnan(cbfDK_regions(:, indicesSign_cbfDK(i)))));
end

for i = 1 : numel(indicesSign_cbfHammers)
    fprintf("Hammers. Region %s. Medians cbf controls %.1f, covid %.1f. Num NaN: %d\n", cbfHammers.Properties.VariableNames{indexFirstRegion+indicesSign_cbfHammers(i)-1}, ...
        median(cbfHammers_regions(indicesControls,indicesSign_cbfHammers(i)), "omitnan"), median(cbfHammers_regions(indicesCovid, indicesSign_cbfHammers(i)), "omitnan"),...
        sum(isnan(cbfHammers_regions(:, indicesSign_cbfHammers(i)))));
end
%%  WHITE MATTER
wmhVolumesPath = '/home/martin/data/UNSAM/CovidProject2/DataAnalysis/AnatomicalMRI/WMH/';
wmhTable = readtable([wmhVolumesPath, 'bianca_data_thr_0_7.0_cluster_5_deep.csv']);
% Leave only the data of the subjects with ASL
indicesAsl = [];
for i = 1 : numel(cbfDK.ID)
    ind=find(strcmp(cbfDK.ID{i}, wmhTable.ID)>0);
    if ~isempty(ind)
       indicesAsl = [indicesAsl ind'];
       wmhMatchedTable(i,:) = wmhTable(ind,:);
    else
        warning(sprintf('Subject %s without data.\n', cbfDK.ID{i}));
        wmhMatchedTable(i,1) = cbfDK.ID(i);
        wmhMatchedTable(i,3) = array2table(NaN);
        wmhMatchedTable(i,4) = array2table(NaN);
    end
end

wmhTable = wmhTable(indicesAsl,:);


%% EXPLORE THE DATA
X = [allCognitiveTests_perc cbfHammers_regions(:, indicesSign_cbfHammers)];
varNames = {'TMT_A'; 'TMT_B'; 'DSF'; 'DSB'; 'Stroop-W'; 'Stroop-C';...
    'Stroop-WC'; 'Stroop-interf'; 'MOCA'};
varNames = cat(1,varNames,cbfDK.Properties.VariableNames(indexFirstRegion+indicesSign_cbfDK-1)');
%cbfDK.Properties.VariableNames{indexFirstRegion+indicesSign_cbfDK(i)-1}
%X(find(sum(isnan(X),2)>0),:) = [];
gplotmatrix(X,[], cbfHammers.Grupo)
figure;
set(gcf, 'Position', [100 100 2600 2400]);
[H,AX,BigAx] = gplotmatrix(X,[],categorical(cbfHammers.Grupo),['b' 'r' ],[],[20],true, 'grpbars', varNames);%'m' 'g' 'r'

%% TRAIN A MODEL WITH RELEVANT DATA FOR REGRESSION
featuresPredictorsOfPerfusion = double([age_years strcmp(gender, 'F') symptom uint8(smoker) uint8(alcoholConsumption) strcmp(healthIssues, 'SI') fas eqolDimensiones allCognitiveTests_perc ...
    wmhMatchedTable.NumberOfClustersOfWMH wmhMatchedTable.TotalVolumeOfClusters]);

output_1 = [cbfHammers_regions(:, indicesSign_cbfHammers(1))];
output_2 = [cbfHammers_regions(:, indicesSign_cbfHammers(2))];
output_5 = [cbfHammers_regions(:, indicesSign_cbfHammers(5))];
output_5 = [cbfHammers_regions(:, indicesSign_cbfHammers(6))];
%% TRAIN A MODEL WITH RELEVANT DATA FOR CLASSIFICATION
featuresPredictorsOfLongCovid = double([age_years strcmp(gender, 'F') symptom uint8(smoker) uint8(alcoholConsumption) strcmp(healthIssues, 'SI') fas eqolDimensiones allCognitiveTests_perc ...
    wmhMatchedTable.NumberOfClustersOfWMH wmhMatchedTable.TotalVolumeOfClusters cbfHammers_regions(:, indicesSign_cbfHammers)]);
featuresPredictorsOfLongCovid = double([cbfHammers_regions wmhMatchedTable.NumberOfClustersOfWMH wmhMatchedTable.TotalVolumeOfClusters]);
group = cbfHammers.Grupo;
% Even data sets
indicesCovidMatched = indicesCovid(randperm(numel(indicesCovid), numel(indicesControls)));
featuresPredictorsOfLongCovidMatched = [featuresPredictorsOfLongCovid(indicesCovidMatched,:); featuresPredictorsOfLongCovid(indicesControls,:)];
groupMatched = [group(indicesCovidMatched); group(indicesControls)];


%% CORRELATION OF PERFUSION 
corrPerf = corr(cbfHammers_regions, featuresPredictorsOfPerfusion);

%% CORRELATION MATRIX
X = [table2array(freesurferAllCognitiveTests_perc) gmVolume wmVolume brainVolume ageFreesurfer_years groupBinary];
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
saveas(gca, [fullFilename], 'png');
close all
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
%% USING PCA TO THE REGIONS
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


%%
mdlStep = stepwiselm(tbl, 'linear', 'Verbose',2);%'MPG ~ Weight'

mdl1 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical');
variables = mdl1.Coefficients.Variables;

mdl2 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical+mocaFreesurfer_score');
variables = mdl2.Coefficients.Variables;

mdl3 = fitlm(tbl, 'gmVolume~ageSienax_years+groupCategorical+mocaFreesurfer_perc');
variables = mdl3.Coefficients.Variables;
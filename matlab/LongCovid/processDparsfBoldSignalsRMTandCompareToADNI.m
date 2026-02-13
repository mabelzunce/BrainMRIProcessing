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
addpath('../fMRI/')
%% DATA PATHS
dataPath = [imagingPartitionPath '/CovidProject/Estudio/PreprocessedMRI/DPARSF/'];
roiSignalsPath = [dataPath '/ResultsAAL/ROISignals_FunImgARWSDCFN/'];
% Subjects data:
studyDataPath = [dataPartitionPath '/UNSAM/CovidProject2/'];
resultsPath = '../DataAnalysis/';
summaryExcelFilename = fullfile(studyDataPath, 'Respuestas.xlsx');
outputPath = [dataPartitionPath '/CovidProject/Estudio/fMRIAnalyisis/RMT/'];
if ~isdir(outputPath)
    mkdir(outputPath)
end
%% ADNI PATH
% ADNI DATA PATHS
adniDataPath = [imagingPartitionPath '/ADNIdata/Trini/'];
adniRoiSignalsPath = [adniDataPath '/ResultsAAL/ROISignals_NiftiPreprocessedAllBatchesNorm/'];
%roiSignalsPath = '/home/martin/data/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/ADNI3_Advanced_MB_fMRI/AAL/';
adniCsvFilename = [adniDataPath 'SubjectsDataAndTests.csv'];
% READ CSV
subjectsInfoAdni = readtable(adniCsvFilename);

% READ ALL THE BOLD SIGNALS, THE DEMOGRAPHICS DATA AND COMPUTE THE CORRELATION MATRICES
subjectsToExclude = {'114_S_6039', '035_S_6953', '128_S_2002', '031_S_4021', '130_S_5231'};
%% SUBJECTS TO EXCLUDE:
% Apply some filters for unwanted volunteers:
subjectsToExclude = {'CP0011', 'CP0015', 'CP0144'};%{'CP0011'}; % 11 (foreign language) and 15 (elder and probaly with some form of dementia). Remove empty names.

%% READ ALL THE BOLD SIGNALS, THE DEMOGRAPHICS DATA AND COMPUTE THE CORRELATION MATRICES
roiSignalsDir = dir([roiSignalsPath,'ROISignals_*.mat']);
% Get filenames, removing dirs
roiSignalsFilenames = {roiSignalsDir(~[roiSignalsDir.isdir]).name};
% group them by clinical group:
indexSubject = 1;
for i = 1 : numel(roiSignalsFilenames)
    roiSignalsTemp = load([roiSignalsPath roiSignalsFilenames{i}]);
    [filepath,name,ext] = fileparts(roiSignalsFilenames{i});
    splitName = split(name,'ROISignals_');
    name = splitName{2};
    
    if sum(strcmp(name, subjectsToExclude)) == 0
        % if nan in the signals, remove them:
        if sum(isnan(roiSignalsTemp.signals)) > 0
            warning(sprintf('Subject %s removed because of nan ROI signals', name));
        else
            subjectNames{indexSubject} = name;
            roiSignalsOriginal{indexSubject} = roiSignalsTemp.signals;    
            lengthSignals(indexSubject) = size(roiSignalsOriginal{indexSubject}, 1);
            % Normalize the signals:
            roiSignals{indexSubject} = zscore(roiSignalsOriginal{indexSubject}); % Works in columns.
            roiSignals{indexSubject}(isnan(roiSignals{indexSubject})) = 0; % Regions with zero signal
            crossCorr(:,:,indexSubject) = corr(roiSignals{indexSubject});
            indexSubject = indexSubject+1;
            %fig = check_fMRI_bold_signals(roiSignals{i});   
        end  
    end
    %fig = check_fMRI_bold_signals(roiSignals{i});       
end
% Correlations between zero signals are NaN, force them to 0.
crossCorr(isnan(crossCorr)) = 0; % Regions with zero signal

%% LOAD MRI REPORTS EXCEL
%responsesTable = readtable(summaryExcelFilename, 'NumHeaderLines', 3);
responsesTable = readtable(summaryExcelFilename,'Sheet', 'respuestasCuestionario', 'Range', 'A4:EJ191', 'ReadVariableNames', false, 'ReadRowNames', false);
% COLUMNS ORDER
% Groups:
colsDemographicData = 4:7;
colsCovidData = 8:36;
colsLongCovidSymptoms = 25:36;
colsHealthHistory = 38:55;

% Leave only the data of the subjects with fMRI
indicesFMRI = [];
for i = 1 : numel(subjectNames)
    ind=find(strcmp(subjectNames{i}, responsesTable.Var1)>0);
    if ~isempty(ind)
       indicesFMRI = [indicesFMRI ind'];
       responsesMatchedTable(i,:) = responsesTable(ind, :);
    else
        warning(sprintf('Subject %s without data.\n', subjectNames{i}));
        responsesMatchedTable(i,1) = subjectNames(i);
        %responsesMatchedTable(i,2:end) = NaN;
        for j = 2 : size(responsesMatchedTable,2)
            if isnumeric(table2array(responsesMatchedTable(i,j)))
                responsesMatchedTable(i,j) = array2table(NaN);
            end
        end
    end
end

responsesTable = responsesTable(indicesFMRI,:);

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
severityCovid = categorical(responsesMatchedTable.Var22, severityNames);
% Symptoms:
symptom = table2array(responsesMatchedTable(:,colsLongCovidSymptoms));
symptomsNames = {'Headaches', 'Fatigue', 'Anosmia', 'Loss of Taste', ...
    'Dyspnoea', 'Muscle weakness', 'Muscle ache', 'Brain fog', 'Communication problems',...
    'Sleeping problems', 'Memory problems', 'Problems with attention'};

labelsSymptoms = {'Yes, and they persist', 'Yes, but no anymore', 'No'};
educationLevels = {'Primario incompleto', 'Primario completo','Secundario incompleto', ...
    'Secundario completo', 'Terciario incompleto', 'Universitario incompleto', 'Terciario completo', ...
    'Universitario completo', 'Posgrado'};
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

%% LOAD ADNI DATA AND SIGNALS
roiSignalsDir = dir(fullfile(adniRoiSignalsPath, '*.mat'));
% Get filenames, removing dirs
roiSignalsFilenames = {roiSignalsDir(~[roiSignalsDir.isdir]).name};
% group them by clinical group:
indexSubject = 1;
for i = 1 : numel(roiSignalsFilenames)
    roiSignalsTemp = load([roiSignalsPath roiSignalsFilenames{i}]);
    [filepath,name,ext] = fileparts(roiSignalsFilenames{i});
    splitName = split(name,'ROISignals_');
    name = splitName{2};
    % Leave only AD and CN
    indexInTable = find(strcmp(name, subjectsInfo.SubjectID));
    if sum(strcmp(name, subjectsToExclude)) == 0
        % if nan in the signals, remove them:
        if sum(isnan(roiSignalsTemp.signals)) > 0
            warning(sprintf('Subject %s removed because of nan ROI signals', name));
        elseif strcmp(subjectsInfoAdni.ResearchGroup(indexInTable), 'AD') || strcmp(subjectsInfoAdni.ResearchGroup(indexInTable), 'CN')
            subjectNamesAdni{indexSubject} = name;
            roiSignalsOriginalAdni{indexSubject} = roiSignalsTemp.signals;    
            lengthSignals(indexSubject) = size(roiSignalsOriginalAdni{indexSubject}, 1);
            % Normalize the signals:
            roiSignalsAdni{indexSubject} = zscore(roiSignalsOriginalAdni{indexSubject}); % Works in columns.
            roiSignalsAdni{indexSubject}(isnan(roiSignalsAdni{indexSubject})) = 0; % Regions with zero signal
            crossCorrAdni(:,:,indexSubject) = corr(roiSignalsAdni{indexSubject});
            indexSubject = indexSubject+1;
            %fig = check_fMRI_bold_signals(roiSignals{i});   
        end  
    end
end
% Correlations between zero signals are NaN, force them to 0.
crossCorrAdni(isnan(crossCorrAdni)) = 0; % Regions with zero signal

% Get demographics of valid data:
for i = 1 : numel(subjectNamesAdni)
% Get demographic and other info:
    indexInTable = find(strcmp(subjectNamesAdni{i}, subjectsInfoAdni.SubjectID));
    sexAdni(i) = categorical(subjectsInfoAdni.Sex(indexInTable));
    groupAdni(i) = categorical(subjectsInfoAdni.ResearchGroup(indexInTable));
    ageAdni(i) = subjectsInfoAdni.Age(indexInTable);
    eduAdni_years(i) = subjectsInfoAdni.PTEDUCAT(indexInTable);
    cdrsbAdni(i) = subjectsInfoAdni.CDRSB(indexInTable);
    mmseAdni(i) = subjectsInfoAdni.MMSE(indexInTable);
    ventriclesAdni(i) = subjectsInfoAdni.Ventricles(indexInTable);
    % site:
    siteAdni{i} = subjectNamesAdni{i}(1:findstr(subjectNamesAdni{i}, '_S')-1);
    % Get the scanner:
    imProtocolAdni = strsplit(subjectsInfoAdni.ImagingProtocol{indexInTable},';');
    aux = strsplit(imProtocol{3},'=');
    scannerAdni{i} = aux{2};
end
% Column vectors:
sexAdni = vec(sexAdni);
groupAdni = vec(groupAdni);
ageAdni = vec(ageAdni);
eduAdni_years = vec(eduAdni_years);
cdrsbAdni = vec(cdrsbAdni);
mmseAdni = vec(mmseAdni);
ventriclesAdni = vec(ventriclesAdni);
siteAdni = vec(categorical(siteAdni));
scannerAdni = vec(categorical(scannerAdni));
%% MERGE DATA
%% VIEW ALL THE CORRELATION MATRICES
% AND SAVE IT IN AN ANIMATED GIF.
filenameGif = [outputPath 'CorrelationMatrices'];
figure;
for i = 1 : numel(subjectNames)
    him = imagesc(crossCorr(:,:,i));
    %cm = him.Parent.Colormap; 
    text(20,10,sprintf('Subject %d. Group %s', i, group{i}), 'FontSize', 14, 'FontWeight', 'bold', 'Color','w')
    f=getframe(gca);
    [im,cm] = rgb2ind(f.cdata,256); 
    % Write to the GIF File 
    if i == 1 
      imwrite(im,cm,filenameGif,'gif', 'Loopcount',inf); 
    else 
      imwrite(im,cm,filenameGif,'gif','WriteMode','append','DelayTime',0.2);
    end
end
%% LINEAR REGRESSION MODEL TO REGRESS OUT COVARIABLES
ageFMRI = vec(age_years);
sexFMRI = vec(categorical(gender));
groupFMRI = vec(categorical(group));
variablesInTheRegression = {'Intercept','Age', 'Sex', 'Group'}; % In the same orhter as the variables are included in the regression.
for i = 1 : size(crossCorr,1)
    for j = 1 : size(crossCorr,1)
        % Fit linear model
        cc_ij = squeeze(crossCorr(i,j,:));
        tbl = table(cc_ij, ageFMRI, sexFMRI, groupFMRI);
        mdl = fitlm(tbl,'cc_ij ~ ageFMRI + sexFMRI + groupFMRI');
        models{i,j} = mdl;
        pvalues(i,j,:) = mdl.Coefficients.pValue; % In the same order as in the model (starting with the intercept);
        coefficients(i,j,:) = mdl.Coefficients.Estimate;
        % Manually regress out variables:
        cc_ij_reg = cc_ij - mdl.Coefficients.Estimate(2)*ageFMRI - mdl.Coefficients.Estimate(3)*(sexFMRI=='M') - mdl.Coefficients.Estimate(4)*(groupFMRI=='COVID');
        % The cc coefficient after regressing out the variables.
        crossCorrRegOut(i,j,:) = cc_ij_reg;
        % IT should be the same as geting the residuals from  the
        % regression.
        crossCorrResiduals(i,j,:) = mdl.Coefficients.Estimate(1) + mdl.Residuals.Raw;
    end
end

%% CHECK THE MODEL
alpha = 0.05;
% p-values
figure;
set(gcf, 'Position', [100 100 2400 1200])
t = tiledlayout('flow','TileSpacing','compact');
for i = 1 : size(coefficients,3)
    nexttile
    imagesc(pvalues(:,:,i), [0 0.5]);
    colorbar
    title(variablesInTheRegression{i})
end
saveas(gca, [outputPath 'pValuesRegression'], 'png');
% Significance
figure;
set(gcf, 'Position', [100 100 2400 1200])
t = tiledlayout('flow','TileSpacing','compact');
for i = 1 : size(coefficients,3)
    nexttile
    imagesc(pvalues(:,:,i)<alpha);
    title(variablesInTheRegression{i})
end
saveas(gca, [outputPath 'significantVariablesRegression'], 'png');
%% COMPUTE MEAN VALUES FOR EACH RANDOM MATRIX
meanCC =  mean(crossCorrResiduals, 3); %mean(crossCorrRegOut, 3);
stdCC = std(crossCorrResiduals, 0, 3);
figure;
set(gcf, 'Position', [100 100 2400 1200])
t = tiledlayout('flow','TileSpacing','compact');
nexttile
imagesc(meanCC);
saveas(gca, [outputPath 'meanCovMatrix'], 'png');
%% FIND BRAIN NETWORKS USING RMT WITH THE MEAN CC MATRIX
% 1) Get the eigenvalues and eigenvectores outside the Marchenk-Pastur
% Compute eigen values and vectors:
[v,lambda]=eig(meanCC);
lambda = diag(lambda);
% Ratio length/ROIs:
N = size(meanCC,1);
L = mean(lengthSignals);
Q=N/L;
% Std dev of the data:
stdDevSignals = 1; % In theory, the singals are z normalised.
% Probability Density Function
lambdaMin = stdDevSignals*(1-sqrt(Q))^2;     % Boundaries -
lambdaMax = stdDevSignals*(1+sqrt(Q))^2;    % Boundaries +

% eigenvalues
indicesBrainNetworks =  abs(lambda) > lambdaMax;
lambdaBrainNetworks = lambda(indicesBrainNetworks);
vBrainNetworks = v(:,indicesBrainNetworks);
numBrainNetworks = numel(lambdaBrainNetworks);
fprintf('Number of brain networks discovered: %d.\n', numBrainNetworks);
fprintf('Eigenvalues of the brain networks: '); disp(lambdaBrainNetworks');


% Plot the distributions and eigenvalues.
numBinsHist = 100;
[histLambdas,lambdaBins] = hist(abs(lambda),[linspace(0,lambdaMax,numBinsHist), [lambdaMax+1:max(lambda)]]);
histLambdas=histLambdas/sum(histLambdas);     % Normalization la densidad
% Theoretical pdf evaluated in the n points between lambda -a y lambda + b
binsMP = linspace(lambdaMin,lambdaMax,numBinsHist);
ft = @(lambda,a,b,c, s) (1./(2*pi*lambda*c*s^(2))).*sqrt((b-lambda).*(lambda-a));
distMarchenkPastur = ft(binsMP, lambdaMin, lambdaMax, Q, stdDevSignals);
distMarchenkPastur = distMarchenkPastur./sum(distMarchenkPastur);
% Plot Marchenko-Pastur
figure;
bar(lambdaBins, histLambdas);
hold on;
plot(binsMP, distMarchenkPastur, 'LineWidth', 3);
ylabel('P(lambda)')


plot(abs(lambda), zeros(1, numel(lambda)), 'x', 'LineWidth', 3, 'MarkerSize',10);
legend('Eigenvalues histogram', sprintf('Distr Marchenko-Pastur Q=N/L=%.1f',Q), 'Eigenvalues')
saveas(gca, [outputPath 'eigenvalues'], 'png');
xlim([0 10])
saveas(gca, [outputPath 'eigenvaluesZoomed'], 'png');
%%
% 2) Get only the most important components for each brain network using
% the IPR.
brainNetworks = zeros(size(vBrainNetworks)); % Brain networks is the autovector with only the 1:IPR components.
iprNetworks = round(1./sum(vBrainNetworks.^4));
fprintf('IPR for each brain network: '); disp(iprNetworks);
% Plot the contrbiutions:
figure;
for i = 1 : numBrainNetworks
    [sortedComponents, indexSortedComponents] = sort(abs(vBrainNetworks(:,i)), 'descend');
    % Get ROIS of this brain networkL
    roisBrainNetworks{i} = indexSortedComponents(1:iprNetworks(i));
    % Define a brainnetowrk keeping the components of the autovector that
    % contribute most using the IPR as threshold.
    brainNetworks(roisBrainNetworks{i},i) = vBrainNetworks(roisBrainNetworks{i},i);
    % Estimate the contribution:
    cumContribution = cumsum(sortedComponents)./sum(sortedComponents);
    cumContributionIpr(i) = cumContribution(iprNetworks(i));
    plot([1:numel(cumContribution)], cumContribution, ':', 'LineWidth', 2);
    hold on;
    legendNetworks{i} = sprintf('Brain network %d (Eigenvalue=%.1f, IPR=%d)', i, lambdaBrainNetworks(i), iprNetworks(i));
end
ylim([0 1])
% Plot the IPR
for i = 1 : numBrainNetworks
    plot(iprNetworks(i), cumContributionIpr(i), 'xk', 'LineWidth', 2, 'MarkerSize', 10);
end
legend(legendNetworks, 'Location', 'SouthEast');

saveas(gca, [outputPath 'ComponentContributionsAndIPR'], 'png');
%% SHOW BRAIN NETWORKS

%% GET THE "STRENGTH OF SIGNAL" OF EACH BRAIN NETWORK
% Whitening matrix
W = diag(lambdaBrainNetworks.^-(1/2))*brainNetworks'; 

% Apply the prwhitenning matrix:
for i = 1 : numel(subjectNames)
    roiSignalsOriginal{i}(isnan(roiSignalsOriginal{i})) = 0;
    S{i} = W*roiSignalsOriginal{i}'; % We use a cell array just in case the signals don't have the same length for different subjects.
    %S{i} = W*roiSignals{i}';
    % Get the rms (normalized by number of time points):
    for j = 1 : numBrainNetworks
        rmsBrainNetowrksSubjects(i,j) = rms(S{i}(j,:))./size(S{1},2);
    end
end

% Compare the rms between the two groups for each brain network using a
% kruskal wallis statistical test
for i = 1 : numBrainNetworks
    p_kw(i) = kruskalwallis(rmsBrainNetowrksSubjects(:,i),groupFMRI, 'off');
    p_anova(i) = anova1(rmsBrainNetowrksSubjects(:,i),groupFMRI, 'off');
end
%% VIOLIN PLOTS
% Get the code for violin plots from: https://github.com/bastibe/Violinplot-Matlab.git
violinPlotLibraryPath = '/home/martin/data/UNSAM/Tools/Violinplot-Matlab/'; % Where you have the library.
addpath(violinPlotLibraryPath);
for i = 1 : numBrainNetworks
    figure;
    violinplot(rmsBrainNetowrksSubjects(:,i),groupFMRI);
    ylabel('RMS Brain Netowrk Signal');
    saveas(gca, [outputPath sprintf('ViolinPlotsRmsBrainNetwork_%d',i)], 'png');
end
%% SHOW BRAIN NETWORKS
% Get AAL Rois first:
aal = niftiread('../fMRI/ROImask.nii');
info=niftiinfo('../fMRI/ROImask.nii');
aal = double(aal);
% MNI 152 to have a structural reference:
mni152 = niftiread('../fMRI/T1.nii');
mni152 = double(mni152);
infoMni152 = niftiinfo('../fMRI/T1.nii');
infoMni152.Datatype='double';
% Now go through the networks:
for i = 1 : numBrainNetworks
    % Create an image for this brain network:
    imageRoisBrainNetowrk{i} = zeros(size(aal));
    imageRoisBrainNetowrk{i}(roisBrainNetworks{i}) = 1;
    % Another with the contribution of each component:
    imageBrainNetowrk{i} = zeros(size(aal));
    for j = 1 : numel(roisBrainNetworks{i})
        imageBrainNetowrk{i}(imageRoisBrainNetowrk{i} == roisBrainNetworks{i}(1)) = brainNetworks(roisBrainNetworks{i}(j),i);
    end
    
    % Show the image slices using a reblue colormap (you nned to have the
    % function in the same path):
    figure;
    cmapRedBlues = redblue();
    % Range of values
    maxValues = max(imageBrainNetowrk{i}(:));
    minValues = min(imageBrainNetowrk{i}(:));

    % Scale the image for visualization with labeloverlay:
    scaledImage = uint8(round((imageBrainNetowrk{i}-minValues)./(maxValues-minValues).*63)+1);

    %labels = labels .* mask;

    for l = 10 : 5 : size(imageBrainNetowrk{i},3)-10
        structMri = normalize(rot90(mni152(:,:,l)), 'range', [0 1]);
        scaledSlice = rot90(scaledImage(:,:,l));

        subplot(3,3,l/5-1)
        imshow(labeloverlay(structMri, scaledSlice,'Colormap',cmapRedBlues));
        
    end
    
    
    aux=round(minValues:((maxValues-minValues)/10):maxValues,2);
    tiks=cell(1,11);
    for k=1:11
        tiks{k}=num2str(aux(k));
    end
    %tiks{6}=num2str(0);
    
    colorbar('Colormap',cmapRedBlues,'Position',[0.89 0.18 0.021 0.7],'TickLabels',tiks,'Limits',[0 max(labels(:))/64]);
    
    saveas(gcf,Slicename);

    close all

end
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
dataPath = [dataPartitionPath '/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/DataBaseDrive-20231031T195909Z-001/DataBase/ADNI_fMRI_screening/AAL/'];
roiSignalsPath = [dataPath '/ROISignals/'];
%roiSignalsPath = '/home/martin/data/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/ADNI3_Advanced_MB_fMRI/AAL/';
csvFilename = [dataPath 'DataBaseSubjects.csv'];
outputPath = [dataPartitionPath '/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/RMT_Analysis/'];
if ~isdir(outputPath)
    mkdir(outputPath)
end
%% READ CSV
subjectsInfo = readtable(csvFilename);

%% READ ALL THE BOLD SIGNALS
roiSignalsDir = dir(roiSignalsPath);
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
    if strcmp(subjectsInfo.ResearchGroup(indexInTable), 'AD') || strcmp(subjectsInfo.ResearchGroup(indexInTable), 'CN')
        subjectNames{indexSubject} = name;
        roiSignals{indexSubject} = roiSignalsTemp.ROISignals;    
        lengthSignals(indexSubject) = size(roiSignals{indexSubject}, 1);
        % Normalize the signals:
        roiSignals{indexSubject} = zscore(roiSignals{indexSubject}); % Works in columns.
        crossCorr(:,:,indexSubject) = corr(roiSignals{indexSubject});
        indexSubject = indexSubject+1;
        %fig = check_fMRI_bold_signals(roiSignals{i});   
    end
    
end

% Get demographics of valid data:
for i = 1 : numel(subjectNames)
% Get demographic and other info:
    indexInTable = find(strcmp(subjectNames{i}, subjectsInfo.SubjectID));
    sex(i) = categorical(subjectsInfo.Sex(indexInTable));
    group(i) = categorical(subjectsInfo.ResearchGroup(indexInTable));
    age(i) = subjectsInfo.Age(indexInTable);
    edu_years(i) = subjectsInfo.PTEDUCAT(indexInTable);
    cdrsb(i) = subjectsInfo.CDRSB(indexInTable);
    mmse(i) = subjectsInfo.MMSE(indexInTable);
    ventricles(i) = subjectsInfo.Ventricles(indexInTable);
    % site:
    site{i} = subjectNames{i}(1:findstr(subjectNames{i}, '_S')-1);
    % Get the scanner:
    imProtocol = strsplit(subjectsInfo.ImagingProtocol{indexInTable},';');
    aux = strsplit(imProtocol{3},'=');
    scanner{i} = aux{2};
end
% Column vectors:
sex = vec(sex);
group = vec(group);
age = vec(age);
edu_years = vec(edu_years);
cdrsb = vec(cdrsb);
mmse = vec(mmse);
ventricles = vec(ventricles);
site = vec(categorical(site));
scanner = vec(categorical(scanner));
%% LINEAR REGRESSION MODEL TO REGRESS OUT COVARIABLES
variablesInTheRegression = {'Intercept','Age', 'Sex', 'Group'}; % In the same orhter as the variables are included in the regression.
for i = 1 : size(crossCorr,1)
    for j = 1 : size(crossCorr,1)
        % Fit linear model
        cc_ij = squeeze(crossCorr(i,j,:));
        tbl = table(cc_ij, age, sex, group);
        mdl = fitlm(tbl,'cc_ij ~ age + sex + group');
        models{i,j} = mdl;
        pvalues(i,j,:) = mdl.Coefficients.pValue; % In the same order as in the model (starting with the intercept);
        coefficients(i,j,:) = mdl.Coefficients.Estimate;
        % Manually regress out variables:
        cc_ij_reg = cc_ij - mdl.Coefficients.Estimate(2)*age - mdl.Coefficients.Estimate(3)*(sex=='F') - mdl.Coefficients.Estimate(4)*(group=='AD');
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
meanCC = mean(crossCorrRegOut, 3);
stdCC = std(crossCorrRegOut, 0, 3);
meanCCRes = mean(crossCorrResiduals, 3);
figure;
set(gcf, 'Position', [100 100 2400 1200])
t = tiledlayout('flow','TileSpacing','compact');
nexttile
imagesc(meanCC);
nexttile
imagesc(meanCCRes);
saveas(gca, [outputPath 'meanCovMatrix'], 'png');
%% FIND BRAIN NETWORKS USING RMT WITH THE MEAN CC MATRIX
% 1) Get the eigenvalues and eigenvectores outside the Marchenk-Pastur
% Compute eigen values and vectors:
[v,lambda]=eig(meanCC);
lambda = diag(lambda);
% Ratio length/ROIs:
N = size(meanCC,1);
L = min(lengthSignals);
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
numBinsHist = 50;
[histLambdas,lambdaBins] = hist(lambda,linspace(lambdaMin,lambdaMax,numBinsHist));
histLambdas=histLambdas/sum(histLambdas);     % Normalization la densidad
% Theoretical pdf evaluated in the n points between lambda -a y lambda + b
ft = @(lambda,a,b,c, s) (1./(2*pi*lambda*c*s^(2))).*sqrt((b-lambda).*(lambda-a));
distMarchenkPastur = ft(lambdaBins, lambdaMin, lambdaMax, Q, stdDevSignals);
distMarchenkPastur = distMarchenkPastur./sum(distMarchenkPastur);
% Plot Marchenko-Pastur
figure;
plot(lambdaBins, distMarchenkPastur);
hold on;
plot(abs(lambda), 0, 'x');
%%
% 2) Get only the most important components for each brain network using
% the IPR:
iprNetworks = round(1./sum(vBrainNetworks.^4));
fprintf('IPR for each brain network: '); disp(iprNetworks);
% Plot the contrbiutions:
figure;
for i = 1 : size(vBrainNetworks,2)
    [sortedComponents, indexSortedComponents] = sort(abs(vBrainNetworks(:,i)), 'descend');
    cumContribution = cumsum(sortedComponents)./sum(sortedComponents);
    cumContributionIpr(i) = cumContribution(iprNetworks(i));
    plot([1:numel(cumContribution)], cumContribution, ':', 'LineWidth', 2);
    hold on;
    legendNetworks{i} = sprintf('Brain network %d (Eigenvalue=%.1f)', i, lambdaBrainNetworks(i));
end
ylim([0 1])
    % Plot the IPR
for i = 1 : size(vBrainNetworks,2)
    plot(iprNetworks(i), cumContributionIpr(i), 'xk', 'LineWidth', 2, 'MarkerSize', 10);
end
legend(legendNetworks, 'Location', 'SouthEast');

saveas(gca, [outputPath 'ComponentContributionsAndIPR'], 'png');
%% SHOW BRAIN NETWORKS

%% GET THE "STRENGTH OF SIGNAL" OF EACH BRAIN NETWORK


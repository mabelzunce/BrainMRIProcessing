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

%% READ ALL THE BOLD SIGNALS, THE DEMOGRAPHICS DATA AND COMPUTE THE CORRELATION MATRICES
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
%% VIEW ALL THE CORRELATION MATRICES
% AND SAVE IT IN AN ANIMATED GIF.
filenameGif = [outputPath 'CorrelationMatrices'];
figure;
for i = 1 : numel(subjectNames)
    him = imagesc(crossCorr(:,:,i));
    %cm = him.Parent.Colormap; 
    text(20,10,sprintf('Subject %d. Group %s', i, group(i)), 'FontSize', 14, 'FontWeight', 'bold', 'Color','w')
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
L = min(lengthSignals);
Q=L/N;
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
    S{i} = W*roiSignals{i}'; % We use a cell array just in case the signals don't have the same length for different subjects.
    % Get the rms (normalized by number of time points):
    for j = 1 : numBrainNetworks
        rmsBrainNetowrksSubjects(i,j) = rms(S{i}(j,:))./size(S{1},2);
    end
end

% Compare the rms between the two groups for each brain network using a
% kruskal wallis statistical test
for i = 1 : numBrainNetworks
    p_kw(i) = kruskalwallis(rmsBrainNetowrksSubjects(:,i),group);
    p_anova(i) = anova1(rmsBrainNetowrksSubjects(:,i),group);
end
%% VIOLIN PLOTS
% Get the code for violin plots from: https://github.com/bastibe/Violinplot-Matlab.git
violinPlotLibraryPath = '/home/martin/data/UNSAM/Tools/Violinplot-Matlab/'; % Where you have the library.
addpath(violinPlotLibraryPath);
for i = 1 : numBrainNetworks
    figure;
    violinplot(rmsBrainNetowrksSubjects(:,i),group);
    ylabel('RMS Brain Netowrk Signal');
    saveas(gca, [outputPath sprintf('ViolinPlotsRmsBrainNetwork_%d',i)], 'png');
end
%% SHOW BRAIN NETWORKS
% Get AAL Rois first:
aal = niftiread('./ROImask.nii');
info=niftiinfo('./ROImask.nii');
aal = double(aal);
% MNI 152 to have a structural reference:
mni152 = niftiread('./T1.nii');
mni152 = double(mni152);
infoMni152 = niftiinfo('./T1.nii');
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
    maxValues = max(imageBrainNetowrk{i});
    minValues = min(imageBrainNetowrk{i});

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
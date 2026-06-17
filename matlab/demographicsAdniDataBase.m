adniDatabaseCsv = '/home/martin/data/UNSAM/CEMSC3/ProcesamientoADNI/DataBase/fMRI_screening_initalvisit/fMRI_screening_initalvisit_2023_04_27.csv';
tableAdniDatabase = readtable(adniDatabaseCsv);
%% Stats
numSubjects = numel(tableAdniDatabase.SubjectID);
figure;
subplot(1,2,2)
hist(categorical(tableAdniDatabase.Visit));
subplot(1,2,1)
hist(categorical(tableAdniDatabase.Phase));
%%
figure;
pie(categorical(tableAdniDatabase.Sex));
figure;
pie(categorical(tableAdniDatabase.ResearchGroup));
%% Age by group
[histGroups, groupCat]=hist(categorical(tableAdniDatabase.ResearchGroup));
for i = 1 : numel(groupCat)
    indexThisGroup = strcmp(tableAdniDatabase.ResearchGroup,groupCat{i});
    meanAge(i) = mean(tableAdniDatabase.Age(indexThisGroup));
    stdAge(i) = std(tableAdniDatabase.Age(indexThisGroup));
end

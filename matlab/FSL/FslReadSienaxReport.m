function [volumes, labels] = FslReadSienaxReport(filenameSienax)
    % Newer matlab versions:
    %sienaxReport = readcell(filenameSienax,  "FileType" ,'text');
    % Older Matlab versions:
    fid = fopen(filenameSienax);
    sienaxReport = textscan(fid,'%s','delimiter','\t');
    sienaxReport = sienaxReport{1};
    fclose(fid);
    
    % greyMatterStrings = splitstring(sienaxReport{14});
    % whiteMatterStrings = splitstring(sienaxReport{15});
    % brainStrings = splitstring(sienaxReport{16});
    % 
    % % Volumes:
    % normGreyMatterVolume = str2num(greyMatterStrings(2,:));
    % normWhitteMatterVolume = str2num(whiteMatterStrings(2,:));
    % normBrainVolume = str2num(brainStrings(2,:));
    % greyMatterVolume = str2num(greyMatterStrings(3,:));
    % whitteMatterVolume = str2num(whiteMatterStrings(3,:));
    % brainVolume = str2num(brainStrings(3,:));
    
    for i = 1:3
        volStrings = strsplit(sienaxReport{13+i});
        for j = 1:2            
            volumes(i,j) = str2num(volStrings{1+j});
        end
        labels{i} = volStrings{1};
    end
end
function FslBet(filenameMri, outputFilename, dTE)
    [filepath,name,ext] = fileparts(filenameMri);
    filenameMri_brain = [filepath '\' name '_brain' ext];
    command = sprintf('bet %s %s -f 0.5 -g 0', filenameMri, filenameMag1_brain);
    if out < 0
        filenameMag1_brain = [];
    end
    return filenameMag1_brain;
end
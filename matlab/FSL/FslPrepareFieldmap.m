function outputFilename = FslPrepareFieldmap(filenamePhase, filenameMag1, outputFilename, dTE)
    % Call BET to extract brain from magnitude:
    [filepath,name,ext] = fileparts(filenameMag1);
    filenameMag1_brain = [filepath '\' name '_brain' ext];
    command = sprintf('bet %s %s -f 0.5 -g 0', filenameMag1, filenameMag1_brain);
    out = system(command);
    % Call command:
    command = sprintf('fsl_prepare_fieldmap SIEMENS %s %s %s %f', ...
        filenamePhase, filenameMag1_brain, outputFilename, dTE);
    out = system(command);
    if out < 0
        outputFilename = [];
    end
    return 
end
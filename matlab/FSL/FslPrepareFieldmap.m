function outputFilename = FslPrepareFieldmap(filenamePhase, filenameMag1, outputPathFieldmap, dTE)
    outputFilename = [outputPathFieldmap 'fmap_rads.nii.gz'];
    % Call BET to extract brain from magnitude:
    [filepath,name,ext] = fileparts(filenameMag1);
    if contains(name, '.')
        [p, name, ext2] = fileparts(name);
    else
        ext2 = '';
    end
    % Use the same outputpath as for the fieldmaps:
    filenameMag1_brain = [outputPathFieldmap '/' name '_brain' ext2 ext];
    command = sprintf('bet %s %s -f 0.55 -g 0', filenameMag1, filenameMag1_brain); % Here is preffered a more agressive Bet
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
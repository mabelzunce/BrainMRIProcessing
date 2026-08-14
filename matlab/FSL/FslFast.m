% FSL Bet
function filenameSeg = FslFast(filenameT1, outputPath, parameters)
    filenameSeg = '';
    [filepath,name,ext] = fileparts(filenameT1);
    if contains(name, '.')
        [p, name, ext2] = fileparts(name);
    end
    % The output using the output path:
    basefilenameOutput = fullfile(outputPath, name);
    command = sprintf('fast %s %s -o %s', parameters, filenameT1, basefilenameOutput);
    out=system(command);
    if exist([basefilenameOutput '_seg.nii.gz'])
        filenameSeg = [basefilenameOutput '_seg.nii.gz'];
    end
end
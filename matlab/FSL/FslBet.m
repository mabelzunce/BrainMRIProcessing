% FSL Bet
function filenameMri_brain = FslBet(filenameMri, outputPath, parameters, centre)
    if nargin == 4
        parameters = [parameters sprintf(' -c %d %d %d', centre(1), centre(2), centre(3))];
    end
    [filepath,name,ext] = fileparts(filenameMri);
    if contains(name, '.')
        [p, name, ext2] = fileparts(name);
    end
    % The output using the output path:
    filenameMri_brain = fullfile(outputPath, [name '_brain']);
    command = sprintf('bet %s %s %s', filenameMri, filenameMri_brain, parameters);
    out=system(command);
    if out < 0
        filenameMri_brain = [];
    end
end
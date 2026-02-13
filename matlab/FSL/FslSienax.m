% FSL Bet
function out = FslSienax(filenameMri, outputSubdir, betOptions, ...
    fastOptions, filenameLesionsMask)
    if nargin == 4
        lesionsMask = [];
        filenameLesionsMask = [];
    end
    [filepath,name,ext] = fileparts(filenameMri);
    if contains(name, '.')
        [p, name, ext2] = fileparts(name);
    end
    % Siena needs to run where the output
    currentPath = pwd;
    cd(filepath);
    if isempty(fastOptions)
        command = sprintf('sienax %s -o %s -B "%s"', name, ...
            outputSubdir, betOptions);
    elseif isempty(filenameLesionsMask)
        command = sprintf('sienax %s -o %s -B "%s" -S "%s"', name, ...
            outputSubdir, betOptions, fastOptions);
    else
        command = sprintf('sienax %s -o %s -B "%s" -S "%s" -lm %s', name, ...
            outputSubdir, betOptions, fastOptions, filenameLesionsMask);
    end
    out=system(command);
    cd(currentPath)
end
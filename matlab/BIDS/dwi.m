root_dir = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Delfina/BIDS/BIDS_MRI/';
subj_folders = dir(fullfile(root_dir, 'sub-CP*'));

sin_cambios = {};
log_cambios = {};

for i = 1:length(subj_folders)
    subj = subj_folders(i).name;
    dwi_path = fullfile(root_dir, subj, 'dwi');
    if ~isfolder(dwi_path), continue; end
    
    dwi_files = dir(fullfile(dwi_path, '*_dwi.nii.gz'));
    nombres_finales = {}; 
    
    for j = 1:length(dwi_files)
        basefile = dwi_files(j).name;
        [~, name_nii, ~] = fileparts(basefile);
        if endsWith(name_nii, '.nii')
            name_nii = extractBefore(name_nii, '.nii');
        end

        json_file = fullfile(dwi_path, [name_nii '.json']);
        bval_file = fullfile(dwi_path, [name_nii '.bval']);
        bvec_file = fullfile(dwi_path, [name_nii '.bvec']);
        nii_file  = fullfile(dwi_path, [name_nii '.nii.gz']);

        if ~isfile(json_file)
            sin_cambios{end+1} = basefile; continue;
        end
        
        txt = fileread(json_file);
        datos = jsondecode(txt);

        % (1) Phase Encoding
        dir_label = '';
        if isfield(datos, 'PhaseEncodingDirection')
            ped = datos.PhaseEncodingDirection;
            switch ped
                case 'i',  dir_label = 'dir-RL';
                case 'j',  dir_label = 'dir-PA';
                case 'j-', dir_label = 'dir-AP';
            end
        end

        % (2) bval
        acq_label = '';
        if isfile(bval_file)
            bvals = load(bval_file);
            uniq_b = unique(bvals);
            if isequal(uniq_b, 0)
                acq_label = 'acq-b0';
            elseif isequal(uniq_b, [0 1000])
                acq_label = 'acq-b1000';
            elseif isequal(uniq_b, 2000)
                acq_label = 'acq-b2000';
            elseif isequal(uniq_b, [0 1000 2000])
                acq_label = 'acq-b1000b2000';
            end
        end

        if isempty(dir_label) && isempty(acq_label)
            sin_cambios{end+1} = basefile;
            continue;
        end

        % (3) Reconstruir nuevo nombre
        partes = split(name_nii, '_');
        nuevas = {};
        for k = 1:length(partes)
            if startsWith(partes{k}, {'acq-', 'dir-'})
                continue; 
            end
            nuevas{end+1} = partes{k};
        end

        if ~isempty(acq_label)
            nuevas{end+1} = acq_label;
        end
        if ~isempty(dir_label)
            nuevas{end+1} = dir_label;
        end
        nuevas{end+1} = 'dwi';

        nuevo_base = strjoin(nuevas, '_');
        nuevo_base_sin_run = regexprep(nuevo_base, '_run-[0-9]+', '');

        if ~any(strcmp(nombres_finales, nuevo_base_sin_run))
            final_base = nuevo_base_sin_run;
        else
            final_base = nuevo_base;  % conservar run-xx
        end

        nombres_finales{end+1} = final_base;

        rename_file(dwi_path, [name_nii '.nii.gz'], [final_base '.nii.gz']);
        rename_file(dwi_path, [name_nii '.json'],   [final_base '.json']);
        rename_file(dwi_path, [name_nii '.bval'],   [final_base '.bval']);
        rename_file(dwi_path, [name_nii '.bvec'],   [final_base '.bvec']);

        log_cambios{end+1} = sprintf('%s → %s', name_nii, final_base);
    end
end

if ~isempty(log_cambios)
    fprintf('Archivos renombrados:\n');
    fprintf('%s\n', log_cambios{:});
end

if ~isempty(sin_cambios)
    fprintf('\n⚠ Archivos sin cambios:\n');
    fprintf('%s\n', sin_cambios{:});
end

function rename_file(folder, oldname, newname)
    oldpath = fullfile(folder, oldname);
    newpath = fullfile(folder, newname);
    if isfile(oldpath)
        movefile(oldpath, newpath);
    end
end



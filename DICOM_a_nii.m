%% convertir todo el dicom a nii

Files = dir('/home/triniibar/Documents/I_parametros_AD/FunRaw/');
% aca se escribe el directorio donde estan las DICOM de TODOS los subjects
% este dir tiene que tener subcarpetas para cada sub

for sub=1:length(Files)
    ext = string(Files(sub).name) ; % extension del archivo
    if (ext~= '.') && (ext~= '..')
        
        path_DICOM = strcat( '/home/triniibar/Documents/I_parametros_AD/FunRaw/', string(ext));
        % el path de la carpeta con los archivos DICOM
        
        path_nii = strcat( '/home/triniibar/Documents/I_parametros_AD/FunImg/', string(ext));
        % el path donde se van a guardar las iamgebes nii
        
        varargout = dicm2nii(path_DICOM, path_nii, 0) ;
        % el ultimo numero deterina el tipo de archivo exportado
        
    end
    
end

% output file format:
%      0 or '.nii'           for single nii uncompressed.
%      1 or '.nii.gz'        for single nii compressed (default).
%      2 or '.hdr'           for hdr/img pair uncompressed.
%      3 or '.hdr.gz'        for hdr/img pair compressed.
%      4 or '.nii 3D'        for 3D nii uncompressed (SPM12).
%      5 or '.nii.gz 3D'     for 3D nii compressed.
%      6 or '.hdr 3D'        for 3D hdr/img pair uncompressed (SPM8).
%      7 or '.hdr.gz 3D'     for 3D hdr/img pair compressed.
%      'bids'                for bids data structure 
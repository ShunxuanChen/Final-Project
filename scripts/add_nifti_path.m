% add_nifti_path.m
% Add NIfTI tools to MATLAB path
folder_path = 'C:\Users\UChi\2025Winter\Neuroimaging\Final\NIfTI_20140122';
if exist(folder_path, 'dir')
    fprintf('Found NIfTI tools folder at: %s\n', folder_path);
    addpath(genpath(folder_path));
    fprintf('Added folder and subfolders to path\n');
    
    % Check if load_nii.m exists
    if exist('load_nii', 'file')
        fprintf('load_nii.m is now on the MATLAB path\n');
    else
        fprintf('WARNING: load_nii.m was not found on the path!\n');
        
        % Search for load_nii.m in the folder
        files = dir(fullfile(folder_path, '**', 'load_nii.m'));
        if ~isempty(files)
            fprintf('Found load_nii.m at: %s\n', fullfile(files(1).folder, files(1).name));
            addpath(files(1).folder);
            fprintf('Added folder containing load_nii.m to path\n');
        else
            fprintf('ERROR: load_nii.m not found in the directory!\n');
        end
    end
else
    fprintf('ERROR: NIfTI tools folder not found at: %s\n', folder_path);
end

% batch_fc_analysis.m
% Script to run functional connectivity analysis for all subjects

% Add NIfTI tools to path
nifti_path = 'C:\Users\UChi\2025Winter\Neuroimaging\Final\NIfTI_20140122';
if exist(nifti_path, 'dir')
    addpath(genpath(nifti_path));
    fprintf('Added NIfTI tools to MATLAB path\n');
else
    error('NIfTI tools folder not found. Please update the path in this script.');
end

% Define paths
data_path_controls = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_controls/ds171_R1.0.0';
data_path_mdd = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_mdd/ds171_R1.0.0';
output_dir = 'functional_connectivity_results';

% Create output directory if needed
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% Find all subjects
control_dirs = dir(fullfile(data_path_controls, 'sub-control*'));
mdd_dirs = dir(fullfile(data_path_mdd, 'sub-mdd*'));

control_subjects = {control_dirs.name};
mdd_subjects = {mdd_dirs.name};

fprintf('Found %d control subjects and %d MDD subjects\n', ...
    length(control_subjects), length(mdd_subjects));

% Process control subjects
for i = 1:length(control_subjects)
    subject = control_subjects{i};
    subject_id = strrep(subject, 'sub-', '');
    fprintf('Processing control subject %s (%d of %d)...\n', subject_id, i, length(control_subjects));
    
    % Get functional files
    func_dir = fullfile(data_path_controls, subject, 'func');
    
    % Try multiple patterns to find music files
    music_files = dir(fullfile(func_dir, '*task-music*run-1*bold.nii.gz'));
    
    % If no files found, try a more flexible pattern
    if isempty(music_files)
        music_files = dir(fullfile(func_dir, '*task-music*bold.nii.gz'));
        if ~isempty(music_files)
            fprintf('  Found music files with alternative pattern\n');
        end
    end
    
    % If still no files found, check directory for all nifti files
    if isempty(music_files)
        all_files = dir(fullfile(func_dir, '*.nii.gz'));
        fprintf('  Warning: No music files found. Directory contains %d .nii.gz files\n', length(all_files));
        continue;
    end
    
    % Process the first found music file
    func_file = fullfile(func_dir, music_files(1).name);
    fprintf('  Using file: %s\n', music_files(1).name);
    
    % Extract time series and compute functional connectivity
    process_subject_fc(func_file, subject_id, output_dir);
end

% Process MDD subjects
for i = 1:length(mdd_subjects)
    subject = mdd_subjects{i};
    subject_id = strrep(subject, 'sub-', '');
    fprintf('Processing MDD subject %s (%d of %d)...\n', subject_id, i, length(mdd_subjects));
    
    % Get functional files
    func_dir = fullfile(data_path_mdd, subject, 'func');
    
    % Try multiple patterns to find music files
    music_files = dir(fullfile(func_dir, '*task-music*run-1*bold.nii.gz'));
    
    % If no files found, try a more flexible pattern
    if isempty(music_files)
        music_files = dir(fullfile(func_dir, '*task-music*bold.nii.gz'));
        if ~isempty(music_files)
            fprintf('  Found music files with alternative pattern\n');
        end
    end
    
    % If still no files found, list all files for debugging
    if isempty(music_files)
        all_files = dir(fullfile(func_dir, '*.nii.gz'));
        fprintf('  Warning: No music files found. Directory contains %d .nii.gz files\n', length(all_files));
        
        % If any files exist, print their names for debugging
        if ~isempty(all_files)
            fprintf('  Files in directory:\n');
            for j = 1:min(5, length(all_files))  % Print first 5 files max
                fprintf('    %s\n', all_files(j).name);
            end
            if length(all_files) > 5
                fprintf('    ... and %d more files\n', length(all_files) - 5);
            end
        end
        continue;
    end
    
    % Process the first found music file
    func_file = fullfile(func_dir, music_files(1).name);
    fprintf('  Using file: %s\n', music_files(1).name);
    
    % Extract time series and compute functional connectivity
    process_subject_fc(func_file, subject_id, output_dir);
end

fprintf('Batch functional connectivity analysis complete.\n');

% Function to process functional connectivity for a single subject
function process_subject_fc(func_file, subject_id, output_dir)
    % Decompress .nii.gz file if needed
    [filepath, filename, ext] = fileparts(func_file);
    if strcmp(ext, '.gz')
        fprintf('  Decompressing file...\n');
        gunzip(func_file, filepath);
        nii_file = fullfile(filepath, filename);
    else
        nii_file = func_file;
    end
    
    % Load NIFTI file
    fprintf('  Loading NIFTI file...\n');
    try
        nii = load_untouch_nii(nii_file);
        data = double(nii.img);
        [nx, ny, nz, nt] = size(data);
        fprintf('  Image dimensions: %d x %d x %d x %d\n', nx, ny, nz, nt);
    catch e
        fprintf('  Error loading NIFTI: %s\n', e.message);
        if strcmp(ext, '.gz') && exist(nii_file, 'file')
            delete(nii_file);
        end
        return;
    end
    
    % Create ROI mask at center of brain (simplified approach)
    center_x = round(nx/2);
    center_y = round(ny/2);
    center_z = round(nz/2);
    radius_vox = 5;
    
    % Create mask
    mask = zeros(nx, ny, nz);
    [x, y, z] = ndgrid(1:nx, 1:ny, 1:nz);
    dist = sqrt((x - center_x).^2 + (y - center_y).^2 + (z - center_z).^2);
    mask(dist <= radius_vox) = 1;
    
    % Extract ROI time series
    roi_indices = find(mask);
    roi_data = zeros(length(roi_indices), nt);
    
    for i = 1:length(roi_indices)
        [x_i, y_i, z_i] = ind2sub([nx, ny, nz], roi_indices(i));
        roi_data(i, :) = squeeze(data(x_i, y_i, z_i, :));
    end
    
    % Calculate mean ROI time series
    roi_ts = mean(roi_data, 1);
    
    % Compute functional connectivity
    fprintf('  Computing functional connectivity...\n');
    
    % Reshape data for computation
    reshaped_data = reshape(data, nx*ny*nz, nt);
    valid_voxels = ~isnan(sum(reshaped_data, 2));
    valid_data = reshaped_data(valid_voxels, :);
    
    % Ensure time series is a row vector
    roi_ts = reshape(roi_ts, 1, numel(roi_ts));
    
    % Compute correlation in batches
    fc_map_flat = zeros(nx*ny*nz, 1);
    batch_size = 5000;
    num_batches = ceil(size(valid_data, 1) / batch_size);
    
    for batch = 1:num_batches
        start_idx = (batch-1) * batch_size + 1;
        end_idx = min(batch * batch_size, size(valid_data, 1));
        batch_indices = start_idx:end_idx;
        
        batch_data = valid_data(batch_indices, :);
        
        % Compute correlations
        for j = 1:size(batch_data, 1)
            voxel_ts = batch_data(j, :);
            
            % Ensure voxel time series is a row vector
            voxel_ts = reshape(voxel_ts, 1, numel(voxel_ts));
            
            % Compute correlation
            try
                c = corrcoef(roi_ts, voxel_ts);
                fc_map_flat(valid_voxels(start_idx + j - 1)) = c(1, 2);
            catch
                % If corrcoef fails, try manual calculation
                roi_norm = (roi_ts - mean(roi_ts)) / std(roi_ts);
                voxel_norm = (voxel_ts - mean(voxel_ts)) / std(voxel_ts);
                fc_map_flat(valid_voxels(start_idx + j - 1)) = sum(roi_norm .* voxel_norm) / length(roi_norm);
            end
        end
    end
    
    % Reshape to 3D
    fc_map = reshape(fc_map_flat, nx, ny, nz);
    
    % Save results
    save(fullfile(output_dir, sprintf('fc_map_%s_music_run1.mat', subject_id)), 'fc_map');
    
    % Clean up temporary files
    if strcmp(ext, '.gz') && exist(nii_file, 'file')
        delete(nii_file);
    end
    
    fprintf('  Functional connectivity analysis for subject %s complete.\n', subject_id);
end
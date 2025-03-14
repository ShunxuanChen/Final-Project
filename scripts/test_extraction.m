% test_extraction.m
% Script to test time series extraction on a single subject

% Add NIfTI tools to path directly
nifti_path = 'C:\Users\UChi\2025Winter\Neuroimaging\Final\NIfTI_20140122';
if exist(nifti_path, 'dir')
    addpath(genpath(nifti_path));
    fprintf('Added NIfTI tools to MATLAB path\n');
else
    error('NIfTI tools folder not found. Please update the path in this script.');
end

% Load ROI information
if exist('roi_masks/roi_info.mat', 'file')
    load('roi_masks/roi_info.mat');
    fprintf('Loaded ROI information from file\n');
else
    fprintf('ERROR: ROI information not found. Please run create_roi_masks.m first.\n');
    return;
end

% Define data paths
data_path_controls = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_controls/ds171_R1.0.0';

% Select first control subject for testing
control_dirs = dir(fullfile(data_path_controls, 'sub-control*'));
subject = control_dirs(1).name;
subject_id = strrep(subject, 'sub-', '');
fprintf('Testing extraction on subject: %s\n', subject_id);

% Get functional file paths
func_dir = fullfile(data_path_controls, subject, 'func');
func_files = dir(fullfile(func_dir, '*bold.nii.gz'));

% Use first functional run for testing
test_file = fullfile(func_dir, func_files(1).name);
fprintf('Using functional file: %s\n', func_files(1).name);

% Decompress file
[filepath, filename, ext] = fileparts(test_file);
if strcmp(ext, '.gz')
    fprintf('Decompressing file...\n');
    gunzip(test_file, filepath);
    nii_file = fullfile(filepath, filename);
else
    nii_file = test_file;
end

% Load the functional data
fprintf('Loading NIFTI file...\n');
nii = load_untouch_nii(nii_file);
data = double(nii.img);
[nx, ny, nz, nt] = size(data);
fprintf('Functional data dimensions: %d x %d x %d x %d\n', nx, ny, nz, nt);

% Create directory for results
results_dir = 'extraction_test_results';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('Created results directory: %s\n', results_dir);
end

% Extract time series for each ROI
roi_ts = cell(length(roi_info.names), 1);

for r = 1:length(roi_info.names)
    roi_name = roi_info.names{r};
    roi_center = roi_info.centers{r};
    
    fprintf('Extracting time series for %s...\n', roi_name);
    
    % Load the ROI mask
    mask_file = fullfile('roi_masks', sprintf('%s_mask.nii', roi_name));
    if exist(mask_file, 'file')
        mask_nii = load_untouch_nii(mask_file);
        mask = mask_nii.img;
    else
        % Create mask on-the-fly if file not found
        fprintf('  Mask file not found, creating mask on-the-fly\n');
        mask = zeros(nx, ny, nz);
        [x, y, z] = ndgrid(1:nx, 1:ny, 1:nz);
        dist = sqrt((x - roi_center(1)).^2 + (y - roi_center(2)).^2 + (z - roi_center(3)).^2);
        mask(dist <= roi_info.radius) = 1;
    end
    
    % Find voxels in the ROI
    roi_voxels = find(mask > 0);
    num_voxels = length(roi_voxels);
    fprintf('  Number of voxels in ROI: %d\n', num_voxels);
    
    % Extract time series from each voxel
    voxel_ts = zeros(num_voxels, nt);
    for i = 1:num_voxels
        [x, y, z] = ind2sub([nx, ny, nz], roi_voxels(i));
        voxel_ts(i, :) = squeeze(data(x, y, z, :));
    end
    
    % Calculate mean time series
    mean_ts = mean(voxel_ts, 1);
    roi_ts{r} = mean_ts;
    
    % Plot the time series
    figure;
    plot(mean_ts);
    title(sprintf('Time Series: %s', roi_name));
    xlabel('Time (TR)');
    ylabel('BOLD Signal');
    saveas(gcf, fullfile(results_dir, sprintf('%s_timeseries.png', roi_name)));
    
    % Save the time series data
    ts_data = struct();
    ts_data.roi_name = roi_name;
    ts_data.mean_ts = mean_ts;
    ts_data.voxel_ts = voxel_ts;
    ts_data.num_voxels = num_voxels;
    save(fullfile(results_dir, sprintf('%s_timeseries.mat', roi_name)), 'ts_data');
end

% Clean up temporary files
if strcmp(ext, '.gz') && exist(nii_file, 'file')
    delete(nii_file);
end

fprintf('Time series extraction test complete. Results saved to %s\n', results_dir);
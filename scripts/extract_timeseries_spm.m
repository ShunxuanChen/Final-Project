% extract_time_series.m
% Simple standalone script to extract time series from a brain region

% Add NIfTI tools to path
addpath(genpath('C:\Users\UChi\2025Winter\Neuroimaging\Final\NIfTI_20140122'));

% Define data path
data_path = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_controls/ds171_R1.0.0';

% Select first subject
subjects = dir(fullfile(data_path, 'sub-control*'));
subject = subjects(1).name;
subject_id = strrep(subject, 'sub-', '');
fprintf('Processing subject: %s\n', subject_id);

% Get first functional file
func_dir = fullfile(data_path, subject, 'func');
func_files = dir(fullfile(func_dir, '*bold.nii.gz'));
func_file = fullfile(func_dir, func_files(1).name);
fprintf('Using file: %s\n', func_files(1).name);

% Create output folder
output_dir = 'time_series_results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Load the data
fprintf('Loading NIFTI file...\n');

% Decompress file if needed
[filepath, filename, ext] = fileparts(func_file);
if strcmp(ext, '.gz')
    fprintf('Decompressing file...\n');
    gunzip(func_file, filepath);
    nii_file = fullfile(filepath, filename);
else
    nii_file = func_file;
end

% Load NIFTI
nii = load_untouch_nii(nii_file);
data = double(nii.img);
[nx, ny, nz, nt] = size(data);
fprintf('Image dimensions: %d x %d x %d x %d\n', nx, ny, nz, nt);

% Create a simple ROI at the center of the brain
center_x = round(nx/2);
center_y = round(ny/2);
center_z = round(nz/2);
radius = 5;  % Radius in voxels

fprintf('Creating ROI centered at [%d, %d, %d] with radius %d voxels\n', ...
    center_x, center_y, center_z, radius);

% Create mask
mask = zeros(nx, ny, nz);
[x, y, z] = ndgrid(1:nx, 1:ny, 1:nz);
dist = sqrt((x - center_x).^2 + (y - center_y).^2 + (z - center_z).^2);
mask(dist <= radius) = 1;

% Count voxels in ROI
roi_voxels = find(mask > 0);
num_voxels = length(roi_voxels);
fprintf('ROI contains %d voxels\n', num_voxels);

% Extract time series from each voxel
voxel_ts = zeros(num_voxels, nt);
for i = 1:num_voxels
    [x_i, y_i, z_i] = ind2sub([nx, ny, nz], roi_voxels(i));
    voxel_ts(i, :) = squeeze(data(x_i, y_i, z_i, :));
end

% Calculate mean time series
mean_ts = mean(voxel_ts, 1);

% Plot the time series
figure;
plot(mean_ts);
title(sprintf('Mean Time Series (%d voxels)', num_voxels));
xlabel('Time (TR)');
ylabel('BOLD Signal');
saveas(gcf, fullfile(output_dir, 'mean_time_series.png'));

% Save the results
results = struct();
results.mean_ts = mean_ts;
results.num_voxels = num_voxels;
results.roi_center = [center_x, center_y, center_z];
results.roi_radius = radius;

save(fullfile(output_dir, 'time_series_results.mat'), 'results');

% Clean up temporary files
if strcmp(ext, '.gz') && exist(nii_file, 'file')
    delete(nii_file);
end

fprintf('Extraction complete. Results saved to %s\n', output_dir);
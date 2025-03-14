% main_extraction_debug.m
% Script to extract time series from ROIs with debug information

% Add NIfTI tools to path
addpath(genpath('C:\Users\UChi\2025Winter\Neuroimaging\Final\NIfTI_20140122'));

% Define paths to your data
data_path_controls = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_controls/ds171_R1.0.0';
data_path_mdd = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_mdd/ds171_R1.0.0';

% Find all subjects
control_dirs = dir(fullfile(data_path_controls, 'sub-control*'));
mdd_dirs = dir(fullfile(data_path_mdd, 'sub-mdd*'));

control_subjects = {control_dirs.name};
mdd_subjects = {mdd_dirs.name};

fprintf('Found %d control subjects and %d MDD subjects\n', ...
    length(control_subjects), length(mdd_subjects));

% Process a single subject to test
subject = control_subjects{1};
subject_id = strrep(subject, 'sub-', '');
fprintf('Processing subject %s...\n', subject_id);

% Get functional files
func_dir = fullfile(data_path_controls, subject, 'func');
func_files = dir(fullfile(func_dir, '*bold.nii.gz'));

% Process one file for testing
func_file = fullfile(func_dir, func_files(1).name);
fprintf('Using functional file: %s\n', func_files(1).name);

% Extract time series
fprintf('Extracting time series...\n');

% Decompress the file
[filepath, filename, ext] = fileparts(func_file);
if strcmp(ext, '.gz')
    fprintf('Decompressing file...\n');
    gunzip(func_file, filepath);
    nii_file = fullfile(filepath, filename);
else
    nii_file = func_file;
end

% Load the data
fprintf('Loading NIFTI file...\n');
try
    nii = load_untouch_nii(nii_file);
    data = double(nii.img);
    [nx, ny, nz, nt] = size(data);
    fprintf('Image dimensions: %d x %d x %d x %d\n', nx, ny, nz, nt);
catch e
    fprintf('Error loading NIFTI file: %s\n', e.message);
    return;
end

% Try a simpler approach - just use a central sphere in the brain
fprintf('Creating ROI mask at center of the brain...\n');

% Define center of the brain
center_x = round(nx/2);
center_y = round(ny/2);
center_z = round(nz/2);
fprintf('Brain center in voxel coordinates: [%d, %d, %d]\n', center_x, center_y, center_z);

% Create mask for spherical ROI (using voxel coordinates)
mask = zeros(nx, ny, nz);
[x, y, z] = ndgrid(1:nx, 1:ny, 1:nz);
radius_vox = 5; % Radius in voxels instead of mm
dist = sqrt((x - center_x).^2 + (y - center_y).^2 + (z - center_z).^2);
mask(dist <= radius_vox) = 1;

% Count voxels in the ROI
num_voxels = sum(mask(:));
fprintf('Number of voxels in ROI: %d\n', num_voxels);

% Check if mask contains any voxels
if num_voxels == 0
    fprintf('ERROR: ROI mask contains no voxels. Try adjusting the radius or center.\n');
    return;
end

% Visualize the mask
slice_x = center_x;
slice_y = center_y;
slice_z = center_z;

figure;
subplot(1,3,1);
imagesc(squeeze(mask(slice_x,:,:))'); title('Sagittal View');
subplot(1,3,2);
imagesc(squeeze(mask(:,slice_y,:))'); title('Coronal View');
subplot(1,3,3);
imagesc(squeeze(mask(:,:,slice_z))); title('Axial View');
sgtitle('ROI Mask');
saveas(gcf, 'roi_mask_visualization.png');

% Extract time series
roi_voxels = find(mask > 0);
voxel_ts = zeros(num_voxels, nt);

fprintf('Extracting time series from %d voxels...\n', num_voxels);
for i = 1:num_voxels
    [x_i, y_i, z_i] = ind2sub([nx, ny, nz], roi_voxels(i));
    voxel_ts(i, :) = squeeze(data(x_i, y_i, z_i, :));
end

% Compute average time series
ts = mean(voxel_ts, 1);

% Plot the time series
figure;
plot(ts);
title('ROI Time Series');
xlabel('Time (TR)');
ylabel('BOLD Signal');
saveas(gcf, 'roi_timeseries_plot.png');

% Save all results for debugging
test_results = struct();
test_results.timeseries = ts;
test_results.voxel_timeseries = voxel_ts;
test_results.mask = mask;
test_results.num_voxels = num_voxels;
test_results.center = [center_x, center_y, center_z];
test_results.radius = radius_vox;

save('test_extraction.mat', 'test_results');

% Clean up
if strcmp(ext, '.gz') && exist(nii_file, 'file')
    delete(nii_file);
end

fprintf('Test extraction complete. Results saved to test_extraction.mat\n');
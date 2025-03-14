% create_roi_masks.m
% Script to create and save ROI mask files

% First ensure NIfTI tools are on the path
run('add_nifti_path.m');

% Define output directory for masks
mask_dir = 'roi_masks';
if ~exist(mask_dir, 'dir')
    mkdir(mask_dir);
    fprintf('Created mask directory: %s\n', mask_dir);
end

% Define data path to get a reference image
data_path_controls = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_controls/ds171_R1.0.0';
ref_subject = dir(fullfile(data_path_controls, 'sub-control*'));
ref_subject = ref_subject(1).name;
func_dir = fullfile(data_path_controls, ref_subject, 'func');
func_files = dir(fullfile(func_dir, '*bold.nii.gz'));
ref_file = fullfile(func_dir, func_files(1).name);

% Decompress reference file
[filepath, filename, ext] = fileparts(ref_file);
if strcmp(ext, '.gz')
    fprintf('Decompressing reference file...\n');
    gunzip(ref_file, filepath);
    nii_file = fullfile(filepath, filename);
else
    nii_file = ref_file;
end

% Load reference image to get dimensions
fprintf('Loading reference image...\n');
ref_nii = load_untouch_nii(nii_file);
[nx, ny, nz, ~] = size(ref_nii.img);
fprintf('Reference image dimensions: %d x %d x %d\n', nx, ny, nz);

% Define ROIs consistently for the entire pipeline
roi_names = {'HeschlsGyrus', 'STG', 'MTG', 'Amygdala'};
roi_centers = {
    [40, 30, 25], % HeschlsGyrus approx location
    [50, 35, 25], % STG approx location
    [55, 40, 20], % MTG approx location
    [35, 20, 15]  % Amygdala approx location
};
roi_radius = 5;

% Create and save each ROI mask
for r = 1:length(roi_names)
    roi_name = roi_names{r};
    roi_center = roi_centers{r};
    
    fprintf('Creating mask for %s...\n', roi_name);
    
    % Create ROI mask
    mask = zeros(nx, ny, nz);
    [x, y, z] = ndgrid(1:nx, 1:ny, 1:nz);
    dist = sqrt((x - roi_center(1)).^2 + (y - roi_center(2)).^2 + (z - roi_center(3)).^2);
    mask(dist <= roi_radius) = 1;
    
    % Count voxels in the ROI
    num_voxels = sum(mask(:));
    fprintf('  Number of voxels in ROI: %d\n', num_voxels);
    
    % Create a new NIfTI structure for the mask
    mask_nii = ref_nii;
    mask_nii.img = mask;
    mask_nii.hdr.dime.dim(1) = 3; % 3D volume
    mask_nii.hdr.dime.dim(5) = 1; % Just one volume
    
    % Save the mask
    mask_file = fullfile(mask_dir, sprintf('%s_mask.nii', roi_name));
    save_untouch_nii(mask_nii, mask_file);
    fprintf('  Saved mask to: %s\n', mask_file);
    
    % Create a figure to visualize the mask
    figure;
    
    % Display three orthogonal slices through the center of the ROI
    subplot(1, 3, 1);
    imagesc(squeeze(mask(:, roi_center(2), :))');
    title(sprintf('%s: Coronal view', roi_name));
    axis equal tight;
    
    subplot(1, 3, 2);
    imagesc(squeeze(mask(roi_center(1), :, :))');
    title(sprintf('%s: Sagittal view', roi_name));
    axis equal tight;
    
    subplot(1, 3, 3);
    imagesc(squeeze(mask(:, :, roi_center(3))));
    title(sprintf('%s: Axial view', roi_name));
    axis equal tight;
    
    % Save the figure
    saveas(gcf, fullfile(mask_dir, sprintf('%s_visualization.png', roi_name)));
end

% Save ROI information for later scripts
roi_info = struct();
roi_info.names = roi_names;
roi_info.centers = roi_centers;
roi_info.radius = roi_radius;
save(fullfile(mask_dir, 'roi_info.mat'), 'roi_info');

% Clean up temporary files
if strcmp(ext, '.gz') && exist(nii_file, 'file')
    delete(nii_file);
end

fprintf('All ROI masks created and saved to %s\n', mask_dir);
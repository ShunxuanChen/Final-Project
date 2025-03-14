% representational_similarity_analysis.m
% Script for performing RSA on musical and non-musical valence processing

% Add NIfTI tools to path
addpath(genpath('C:\Users\UChi\2025Winter\Neuroimaging\Final\NIfTI_20140122'));

% Define output directory
output_dir = 'rsa_results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% Define paths to your data
data_path_controls = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_controls/ds171_R1.0.0';
data_path_mdd = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_mdd/ds171_R1.0.0';

% Define ROIs for RSA analysis
roi_names = {'HeschlsGyrus', 'STG', 'MTG', 'Amygdala'};
roi_sizes = [5, 5, 5, 4]; % Radius of each ROI in voxels

% Define ROI center coordinates (simplified - using voxel coordinates)
roi_centers = {
    [40, 30, 25], % HeschlsGyrus approx location
    [50, 35, 25], % STG approx location
    [55, 40, 20], % MTG approx location
    [35, 20, 15]  % Amygdala approx location
};

% Define model RDMs
% 1. Modality model: music vs. non-music
modality_model = [
    0 0 0 1 1; % music-positive similar to music-positive
    0 0 0 1 1; % music-positive similar to music-negative
    0 0 0 1 1; % music-negative similar to music-negative
    1 1 1 0 0; % non-music-positive similar to non-music-positive
    1 1 1 0 0  % non-music-positive similar to non-music-negative
];

% 2. Valence model: positive vs. negative
valence_model = [
    0 1 1 0 1; % positive similar to positive
    1 0 0 1 0; % positive different from negative
    1 0 0 1 0; % negative similar to negative
    0 1 1 0 1; % positive similar to positive
    1 0 0 1 0  % negative similar to negative
];

% Define experimental conditions
conditions = {'music_positive', 'music_negative', 'nonmusic_positive', 'nonmusic_negative'};

% For this example, just process one control subject
control_dirs = dir(fullfile(data_path_controls, 'sub-control*'));
subject = control_dirs(1).name;
subject_id = strrep(subject, 'sub-', '');
fprintf('Processing subject %s for RSA analysis...\n', subject_id);

% Get functional run files
func_dir = fullfile(data_path_controls, subject, 'func');
func_files = dir(fullfile(func_dir, '*bold.nii.gz'));

% Map runs to conditions - note: this is an assumption about your data organization
% In a real analysis, you would need to know which runs correspond to which conditions
condition_files = {
    fullfile(func_dir, func_files(1).name), % music_positive
    fullfile(func_dir, func_files(2).name), % music_negative
    fullfile(func_dir, func_files(3).name), % nonmusic_positive
    fullfile(func_dir, func_files(4).name)  % nonmusic_negative
};

% Pre-allocate arrays for results
subject_rdms = cell(length(roi_names), 1);
modality_cors = zeros(length(roi_names), 1);
valence_cors = zeros(length(roi_names), 1);

% Process each ROI
for r = 1:length(roi_names)
    roi_name = roi_names{r};
    roi_center = roi_centers{r};
    roi_radius = roi_sizes(r);
    
    fprintf('Processing ROI: %s\n', roi_name);
    
    % Pre-allocate matrix for ROI patterns across conditions
    roi_patterns = cell(length(conditions), 1);
    
    % Process each condition
    for c = 1:length(conditions)
        condition = conditions{c};
        func_file = condition_files{c};
        
        fprintf('  Processing condition: %s\n', condition);
        
        % Decompress .nii.gz file
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
        nii = load_untouch_nii(nii_file);
        data = double(nii.img);
        [nx, ny, nz, nt] = size(data);
        
        % Create ROI mask
        mask = zeros(nx, ny, nz);
        [x, y, z] = ndgrid(1:nx, 1:ny, 1:nz);
        dist = sqrt((x - roi_center(1)).^2 + (y - roi_center(2)).^2 + (z - roi_center(3)).^2);
        mask(dist <= roi_radius) = 1;
        
        % Extract voxel-wise patterns within ROI
        roi_indices = find(mask);
        roi_pattern = zeros(length(roi_indices), nt);
        
        for i = 1:length(roi_indices)
            [x_i, y_i, z_i] = ind2sub([nx, ny, nz], roi_indices(i));
            roi_pattern(i, :) = squeeze(data(x_i, y_i, z_i, :));
        end
        
        % Store the pattern
        roi_patterns{c} = roi_pattern;
        
        % Clean up
        if strcmp(ext, '.gz') && exist(nii_file, 'file')
            delete(nii_file);
        end
    end
    
    % Compute neural RDM (correlation distance between patterns)
    n_conditions = length(conditions);
    neural_rdm = zeros(n_conditions, n_conditions);
    
    for i = 1:n_conditions
        for j = 1:n_conditions
            if i == j
                neural_rdm(i, j) = 0; % Zero distance to self
            else
                % Average patterns over time
                pattern_i = mean(roi_patterns{i}, 2);
                pattern_j = mean(roi_patterns{j}, 2);
                
                % Compute correlation
                corr_ij = corr(pattern_i, pattern_j);
                
                % Convert to distance (1 - correlation)
                neural_rdm(i, j) = 1 - corr_ij;
            end
        end
    end
    
    % Store the neural RDM
    subject_rdms{r} = neural_rdm;
    
    % Compare with model RDMs
    % Flatten RDMs to vectors for correlation
    neural_rdm_vec = neural_rdm(triu(true(n_conditions), 1));
    modality_model_vec = modality_model(triu(true(n_conditions), 1));
    valence_model_vec = valence_model(triu(true(n_conditions), 1));
    
    % Compute Spearman's rank correlation
    modality_cors(r) = corr(neural_rdm_vec, modality_model_vec, 'type', 'Spearman');
    valence_cors(r) = corr(neural_rdm_vec, valence_model_vec, 'type', 'Spearman');
    
    fprintf('  Results for %s:\n', roi_name);
    fprintf('    Correlation with modality model: %.4f\n', modality_cors(r));
    fprintf('    Correlation with valence model: %.4f\n', valence_cors(r));
    
    % Visualize the neural RDM
    figure;
    imagesc(neural_rdm);
    colorbar;
    title(sprintf('Neural RDM: %s', roi_name));
    xlabel('Condition');
    ylabel('Condition');
    set(gca, 'XTick', 1:n_conditions, 'XTickLabel', conditions);
    set(gca, 'YTick', 1:n_conditions, 'YTickLabel', conditions);
    colormap('jet');
    saveas(gcf, fullfile(output_dir, sprintf('rdm_%s_%s.png', subject_id, roi_name)));
end

% Save all results
rsa_results = struct();
rsa_results.subject_id = subject_id;
rsa_results.roi_names = roi_names;
rsa_results.conditions = conditions;
rsa_results.neural_rdms = subject_rdms;
rsa_results.modality_correlation = modality_cors;
rsa_results.valence_correlation = valence_cors;

save(fullfile(output_dir, sprintf('rsa_results_%s.mat', subject_id)), 'rsa_results');

fprintf('RSA analysis complete. Results saved to %s\n', output_dir);
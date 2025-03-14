% batch_rsa_analysis.m
% Script to perform RSA for all subjects

% Add NIfTI tools to path directly instead of using a separate script
nifti_path = 'C:\Users\UChi\2025Winter\Neuroimaging\Final\NIfTI_20140122';
if exist(nifti_path, 'dir')
    addpath(genpath(nifti_path));
    fprintf('Added NIfTI tools to MATLAB path\n');
else
    error('NIfTI tools folder not found. Please update the path in this script.');
end

% Define data paths
data_path_controls = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_controls/ds171_R1.0.0';
data_path_mdd = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_mdd/ds171_R1.0.0';

% Define output directory
output_dir = 'rsa_results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% Load ROI information if available, otherwise define default ROIs
if exist('roi_masks/roi_info.mat', 'file')
    load('roi_masks/roi_info.mat');
    fprintf('Loaded ROI information from file\n');
    roi_names = roi_info.names;
    roi_centers = roi_info.centers;
    roi_radius = roi_info.radius;
else
    fprintf('ROI information file not found. Using default ROIs.\n');
    roi_names = {'HeschlsGyrus', 'STG', 'MTG', 'Amygdala'};
    roi_centers = {
        [40, 30, 25], % HeschlsGyrus approx location
        [50, 35, 25], % STG approx location
        [55, 40, 20], % MTG approx location
        [35, 20, 15]  % Amygdala approx location
    };
    roi_radius = 5;
end

% Define model RDMs
% 1. Modality model: music vs. non-music
modality_model = [
    0, 0, 1, 1;   % music-positive similar to music-positive and music-negative
    0, 0, 1, 1;   % music-negative similar to music-positive and music-negative
    1, 1, 0, 0;   % nonmusic-positive similar to nonmusic-positive and nonmusic-negative
    1, 1, 0, 0    % nonmusic-negative similar to nonmusic-positive and nonmusic-negative
];

% 2. Valence model: positive vs. negative
valence_model = [
    0, 1, 0, 1;   % positive similar to positive, different from negative
    1, 0, 1, 0;   % negative similar to negative, different from positive
    0, 1, 0, 1;   % positive similar to positive, different from negative
    1, 0, 1, 0    % negative similar to negative, different from positive
];

% Define condition labels
conditions = {'music_positive', 'music_negative', 'nonmusic_positive', 'nonmusic_negative'};

% Get subjects
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
    
    % Get functional run files
    func_dir = fullfile(data_path_controls, subject, 'func');
    
    % Map conditions to runs - this is critical for RSA
    % Assuming a standard pattern in the dataset
    music_pos_files = dir(fullfile(func_dir, '*task-music*run-1*bold.nii.gz'));
    music_neg_files = dir(fullfile(func_dir, '*task-music*run-2*bold.nii.gz'));
    nonmusic_pos_files = dir(fullfile(func_dir, '*task-nonmusic*run-4*bold.nii.gz'));
    nonmusic_neg_files = dir(fullfile(func_dir, '*task-nonmusic*run-5*bold.nii.gz'));
    
    % Make sure we have all the needed files
    if isempty(music_pos_files) || isempty(music_neg_files) || ...
       isempty(nonmusic_pos_files) || isempty(nonmusic_neg_files)
        fprintf('  Warning: Missing condition files for subject %s\n', subject_id);
        continue;
    end
    
    % Collect all condition files
    condition_files = {
        fullfile(func_dir, music_pos_files(1).name),     % music_positive
        fullfile(func_dir, music_neg_files(1).name),     % music_negative
        fullfile(func_dir, nonmusic_pos_files(1).name),  % nonmusic_positive
        fullfile(func_dir, nonmusic_neg_files(1).name)   % nonmusic_negative
    };
    
    % Process each ROI
    subject_results = struct();
    subject_results.subject_id = subject_id;
    subject_results.group = 'control';
    subject_results.conditions = conditions;
    subject_results.roi_names = roi_names;
    subject_results.neural_rdms = cell(length(roi_names), 1);
    subject_results.modality_correlation = zeros(length(roi_names), 1);
    subject_results.valence_correlation = zeros(length(roi_names), 1);
    
    for r = 1:length(roi_names)
        roi_name = roi_names{r};
        roi_center = roi_centers{r};
        
        fprintf('  Processing ROI: %s\n', roi_name);
        
        % Check if ROI mask exists
        mask_file = fullfile('roi_masks', sprintf('%s_mask.nii', roi_name));
        
        % Create conditions pattern array
        condition_patterns = cell(length(conditions), 1);
        
        % Process each condition
        for c = 1:length(conditions)
            condition = conditions{c};
            func_file = condition_files{c};
            
            fprintf('    Processing condition: %s\n', condition);
            
            % Extract voxel patterns for this condition
            % Decompress file if needed
            [filepath, filename, ext] = fileparts(func_file);
            if strcmp(ext, '.gz')
                fprintf('    Decompressing file...\n');
                gunzip(func_file, filepath);
                nii_file = fullfile(filepath, filename);
            else
                nii_file = func_file;
            end
            
            % Load NIFTI
            fprintf('    Loading NIFTI file...\n');
            try
                nii = load_untouch_nii(nii_file);
                data = double(nii.img);
                [nx, ny, nz, nt] = size(data);
            catch e
                fprintf('    Error loading NIFTI: %s\n', e.message);
                continue;
            end
            
            % Create ROI mask (either load from file or create on-the-fly)
            if exist(mask_file, 'file')
                mask_nii = load_untouch_nii(mask_file);
                mask = mask_nii.img;
            else
                % Create mask on-the-fly
                mask = zeros(nx, ny, nz);
                [x, y, z] = ndgrid(1:nx, 1:ny, 1:nz);
                dist = sqrt((x - roi_center(1)).^2 + (y - roi_center(2)).^2 + (z - roi_center(3)).^2);
                mask(dist <= roi_radius) = 1;
            end
            
            % Find voxels in the ROI
            roi_voxels = find(mask > 0);
            if isempty(roi_voxels)
                fprintf('    Warning: No voxels found in ROI.\n');
                continue;
            end
            
            % Extract voxel patterns (average over time)
            voxel_patterns = zeros(length(roi_voxels), 1);
            for v = 1:length(roi_voxels)
                [x_i, y_i, z_i] = ind2sub([nx, ny, nz], roi_voxels(v));
                voxel_ts = squeeze(data(x_i, y_i, z_i, :));
                voxel_patterns(v) = mean(voxel_ts);
            end
            
            % Store the pattern
            condition_patterns{c} = voxel_patterns;
            
            % Clean up
            if strcmp(ext, '.gz') && exist(nii_file, 'file')
                delete(nii_file);
            end
        end
        
        % Compute neural RDM
        n_conditions = length(conditions);
        neural_rdm = zeros(n_conditions, n_conditions);
        
        for c1 = 1:n_conditions
            for c2 = 1:n_conditions
                if c1 == c2
                    neural_rdm(c1, c2) = 0; % No distance to self
                else
                    % Check if patterns exist
                    if isempty(condition_patterns{c1}) || isempty(condition_patterns{c2})
                        neural_rdm(c1, c2) = NaN;
                        continue;
                    end
                    
                    % Compute correlation
                    try
                        % Ensure patterns are column vectors
                        pattern1 = reshape(condition_patterns{c1}, [], 1);
                        pattern2 = reshape(condition_patterns{c2}, [], 1);
                        
                        c_val = corrcoef(pattern1, pattern2);
                        neural_rdm(c1, c2) = 1 - c_val(1, 2); % Convert to distance
                    catch
                        neural_rdm(c1, c2) = NaN;
                    end
                end
            end
        end
        
        % Store neural RDM
        subject_results.neural_rdms{r} = neural_rdm;
        
        % Compare with model RDMs
        % Convert RDMs to vector form (upper triangle)
        mask = triu(true(n_conditions), 1);
        neural_vec = neural_rdm(mask);
        modality_vec = modality_model(mask);
        valence_vec = valence_model(mask);
        
        % Remove NaN values if any
        valid = ~isnan(neural_vec);
        if sum(valid) > 0
            neural_vec = neural_vec(valid);
            modality_vec = modality_vec(valid);
            valence_vec = valence_vec(valid);
            
            % Calculate Spearman correlation with models
            subject_results.modality_correlation(r) = corr(neural_vec, modality_vec, 'type', 'Spearman');
            subject_results.valence_correlation(r) = corr(neural_vec, valence_vec, 'type', 'Spearman');
        else
            subject_results.modality_correlation(r) = NaN;
            subject_results.valence_correlation(r) = NaN;
        end
        
        fprintf('    Results for %s:\n', roi_name);
        fprintf('      Correlation with modality model: %.4f\n', subject_results.modality_correlation(r));
        fprintf('      Correlation with valence model: %.4f\n', subject_results.valence_correlation(r));
    end
    
    % Save subject results
    save(fullfile(output_dir, sprintf('rsa_results_%s.mat', subject_id)), 'subject_results');
end

% Process MDD subjects using the same approach
for i = 1:length(mdd_subjects)
    subject = mdd_subjects{i};
    subject_id = strrep(subject, 'sub-', '');
    fprintf('Processing MDD subject %s (%d of %d)...\n', subject_id, i, length(mdd_subjects));
    
    % Get functional run files
    func_dir = fullfile(data_path_mdd, subject, 'func');
    
    % Map conditions to runs
    music_pos_files = dir(fullfile(func_dir, '*task-music*run-1*bold.nii.gz'));
    music_neg_files = dir(fullfile(func_dir, '*task-music*run-2*bold.nii.gz'));
    nonmusic_pos_files = dir(fullfile(func_dir, '*task-nonmusic*run-4*bold.nii.gz'));
    nonmusic_neg_files = dir(fullfile(func_dir, '*task-nonmusic*run-5*bold.nii.gz'));
    
    % Make sure we have all the needed files
    if isempty(music_pos_files) || isempty(music_neg_files) || ...
       isempty(nonmusic_pos_files) || isempty(nonmusic_neg_files)
        fprintf('  Warning: Missing condition files for subject %s\n', subject_id);
        continue;
    end
    
    % Collect all condition files
    condition_files = {
        fullfile(func_dir, music_pos_files(1).name),     % music_positive
        fullfile(func_dir, music_neg_files(1).name),     % music_negative
        fullfile(func_dir, nonmusic_pos_files(1).name),  % nonmusic_positive
        fullfile(func_dir, nonmusic_neg_files(1).name)   % nonmusic_negative
    };
    
    % Process each ROI
    subject_results = struct();
    subject_results.subject_id = subject_id;
    subject_results.group = 'mdd';
    subject_results.conditions = conditions;
    subject_results.roi_names = roi_names;
    subject_results.neural_rdms = cell(length(roi_names), 1);
    subject_results.modality_correlation = zeros(length(roi_names), 1);
    subject_results.valence_correlation = zeros(length(roi_names), 1);
    
    for r = 1:length(roi_names)
        roi_name = roi_names{r};
        roi_center = roi_centers{r};
        
        fprintf('  Processing ROI: %s\n', roi_name);
        
        % Check if ROI mask exists
        mask_file = fullfile('roi_masks', sprintf('%s_mask.nii', roi_name));
        
        % Create conditions pattern array
        condition_patterns = cell(length(conditions), 1);
        
        % Process each condition
        for c = 1:length(conditions)
            condition = conditions{c};
            func_file = condition_files{c};
            
            fprintf('    Processing condition: %s\n', condition);
            
            % Extract voxel patterns for this condition
            % Decompress file if needed
            [filepath, filename, ext] = fileparts(func_file);
            if strcmp(ext, '.gz')
                fprintf('    Decompressing file...\n');
                gunzip(func_file, filepath);
                nii_file = fullfile(filepath, filename);
            else
                nii_file = func_file;
            end
            
            % Load NIFTI
            fprintf('    Loading NIFTI file...\n');
            try
                nii = load_untouch_nii(nii_file);
                data = double(nii.img);
                [nx, ny, nz, nt] = size(data);
            catch e
                fprintf('    Error loading NIFTI: %s\n', e.message);
                continue;
            end
            
            % Create ROI mask (either load from file or create on-the-fly)
            if exist(mask_file, 'file')
                mask_nii = load_untouch_nii(mask_file);
                mask = mask_nii.img;
            else
                % Create mask on-the-fly
                mask = zeros(nx, ny, nz);
                [x, y, z] = ndgrid(1:nx, 1:ny, 1:nz);
                dist = sqrt((x - roi_center(1)).^2 + (y - roi_center(2)).^2 + (z - roi_center(3)).^2);
                mask(dist <= roi_radius) = 1;
            end
            
            % Find voxels in the ROI
            roi_voxels = find(mask > 0);
            if isempty(roi_voxels)
                fprintf('    Warning: No voxels found in ROI.\n');
                continue;
            end
            
            % Extract voxel patterns (average over time)
            voxel_patterns = zeros(length(roi_voxels), 1);
            for v = 1:length(roi_voxels)
                [x_i, y_i, z_i] = ind2sub([nx, ny, nz], roi_voxels(v));
                voxel_ts = squeeze(data(x_i, y_i, z_i, :));
                voxel_patterns(v) = mean(voxel_ts);
            end
            
            % Store the pattern
            condition_patterns{c} = voxel_patterns;
            
            % Clean up
            if strcmp(ext, '.gz') && exist(nii_file, 'file')
                delete(nii_file);
            end
        end
        
        % Compute neural RDM
        n_conditions = length(conditions);
        neural_rdm = zeros(n_conditions, n_conditions);
        
        for c1 = 1:n_conditions
            for c2 = 1:n_conditions
                if c1 == c2
                    neural_rdm(c1, c2) = 0; % No distance to self
                else
                    % Check if patterns exist
                    if isempty(condition_patterns{c1}) || isempty(condition_patterns{c2})
                        neural_rdm(c1, c2) = NaN;
                        continue;
                    end
                    
                    % Compute correlation
                    try
                        % Ensure patterns are column vectors
                        pattern1 = reshape(condition_patterns{c1}, [], 1);
                        pattern2 = reshape(condition_patterns{c2}, [], 1);
                        
                        c_val = corrcoef(pattern1, pattern2);
                        neural_rdm(c1, c2) = 1 - c_val(1, 2); % Convert to distance
                    catch
                        neural_rdm(c1, c2) = NaN;
                    end
                end
            end
        end
        
        % Store neural RDM
        subject_results.neural_rdms{r} = neural_rdm;
        
        % Compare with model RDMs
        % Convert RDMs to vector form (upper triangle)
        mask = triu(true(n_conditions), 1);
        neural_vec = neural_rdm(mask);
        modality_vec = modality_model(mask);
        valence_vec = valence_model(mask);
        
        % Remove NaN values if any
        valid = ~isnan(neural_vec);
        if sum(valid) > 0
            neural_vec = neural_vec(valid);
            modality_vec = modality_vec(valid);
            valence_vec = valence_vec(valid);
            
            % Calculate Spearman correlation with models
            subject_results.modality_correlation(r) = corr(neural_vec, modality_vec, 'type', 'Spearman');
            subject_results.valence_correlation(r) = corr(neural_vec, valence_vec, 'type', 'Spearman');
        else
            subject_results.modality_correlation(r) = NaN;
            subject_results.valence_correlation(r) = NaN;
        end
        
        fprintf('    Results for %s:\n', roi_name);
        fprintf('      Correlation with modality model: %.4f\n', subject_results.modality_correlation(r));
        fprintf('      Correlation with valence model: %.4f\n', subject_results.valence_correlation(r));
    end
    
    % Save subject results
    save(fullfile(output_dir, sprintf('rsa_results_%s.mat', subject_id)), 'subject_results');
end

fprintf('RSA analysis complete. Results saved to %s\n', output_dir);
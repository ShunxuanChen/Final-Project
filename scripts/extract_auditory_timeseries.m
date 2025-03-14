% extract_auditory_timeseries.m
% Script to extract time series from auditory cortex ROIs

% Define paths to your data
data_path_controls = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_controls/ds171_R1.0.0';
data_path_mdd = 'C:/Users/UChi/2025Winter/Neuroimaging/Final/ds171_R1.0.0_mdd/ds171_R1.0.0';

% Load subject directories
control_dirs = dir(fullfile(data_path_controls, 'sub-control*'));
mdd_dirs = dir(fullfile(data_path_mdd, 'sub-mdd*'));

control_subjects = {control_dirs.name};
mdd_subjects = {mdd_dirs.name};

fprintf('Found %d control subjects and %d MDD subjects\n', ...
    length(control_subjects), length(mdd_subjects));

% Check if SPM is available
try
    spm('Ver');
    fprintf('SPM is available.\n');
    use_spm = true;
catch
    fprintf('SPM not found. Will use alternative methods.\n');
    use_spm = false;
end

% Define ROIs for auditory processing
% We'll use coordinates from the literature for key auditory regions
roi_names = {'HeschlsGyrus_L', 'HeschlsGyrus_R', 'STG_L', 'STG_R'};
roi_coords = {
    [-42, -22, 7],  % Left Heschl's Gyrus
    [46, -14, 8],   % Right Heschl's Gyrus
    [-58, -22, 6],  % Left Superior Temporal Gyrus
    [58, -10, 4]    % Right Superior Temporal Gyrus
};
roi_radius = 8; % 8mm sphere

% Create structure to store all time series data
all_data = struct();
all_data.control = struct();
all_data.mdd = struct();

% Create temporary directory for ROI masks
mask_dir = fullfile(pwd, 'roi_masks');
if ~exist(mask_dir, 'dir')
    mkdir(mask_dir);
end

% Extract time series for each subject and condition
% First, let's process control subjects
for i = 1:length(control_subjects)
    subject = control_subjects{i};
    subject_id = strrep(subject, 'sub-', '');
    fprintf('Processing control subject %s (%d of %d)...\n', subject_id, i, length(control_subjects));
    
    % Create subject data structure
    all_data.control.(subject_id) = struct();
    
    % Get functional files for this subject
    func_dir = fullfile(data_path_controls, subject, 'func');
    music_files = dir(fullfile(func_dir, '*task-music*bold.nii.gz'));
    nonmusic_files = dir(fullfile(func_dir, '*task-nonmusic*bold.nii.gz'));
    
    % Process music runs
    all_data.control.(subject_id).music = process_runs(func_dir, music_files, roi_names, roi_coords, roi_radius, use_spm);
    
    % Process non-music runs
    all_data.control.(subject_id).nonmusic = process_runs(func_dir, nonmusic_files, roi_names, roi_coords, roi_radius, use_spm);
end

% Process MDD subjects
for i = 1:length(mdd_subjects)
    subject = mdd_subjects{i};
    subject_id = strrep(subject, 'sub-', '');
    fprintf('Processing MDD subject %s (%d of %d)...\n', subject_id, i, length(mdd_subjects));
    
    % Create subject data structure
    all_data.mdd.(subject_id) = struct();
    
    % Get functional files for this subject
    func_dir = fullfile(data_path_mdd, subject, 'func');
    music_files = dir(fullfile(func_dir, '*task-music*bold.nii.gz'));
    nonmusic_files = dir(fullfile(func_dir, '*task-nonmusic*bold.nii.gz'));
    
    % Process music runs
    all_data.mdd.(subject_id).music = process_runs(func_dir, music_files, roi_names, roi_coords, roi_radius, use_spm);
    
    % Process non-music runs
    all_data.mdd.(subject_id).nonmusic = process_runs(func_dir, nonmusic_files, roi_names, roi_coords, roi_radius, use_spm);
end

% Save the results
save('auditory_timeseries_data.mat', 'all_data');
fprintf('Analysis complete. Data saved to auditory_timeseries_data.mat\n');

% Function to process runs and extract time series
function run_data = process_runs(func_dir, run_files, roi_names, roi_coords, roi_radius, use_spm)
    run_data = struct();
    
    for r = 1:length(run_files)
        run_file = fullfile(run_files(r).folder, run_files(r).name);
        run_id = sprintf('run%d', r);
        
        % Extract run name from filename for better identification
        [~, run_name, ~] = fileparts(run_files(r).name);
        run_name = strrep(run_name, '.nii', ''); % Remove extension
        
        run_data.(run_id) = struct('name', run_name);
        
        % Extract time series for each ROI
        for roi_idx = 1:length(roi_names)
            roi_name = roi_names{roi_idx};
            roi_coord = roi_coords{roi_idx};
            
            % Extract time series using appropriate method
            if use_spm
                ts = extract_timeseries_spm(run_file, roi_coord, roi_radius);
            else
                ts = extract_timeseries_basic(run_file, roi_coord, roi_radius);
            end
            
            run_data.(run_id).(roi_name) = ts;
        end
    end
end

% Simple placeholder function - if you have SPM, I'll provide the actual implementation
function ts = extract_timeseries_spm(func_file, roi_coord, roi_radius)
    % Placeholder - This will be implemented based on your available tools
    ts = zeros(1, 100); % Dummy data
    fprintf('  Extracting time series using SPM from %s\n', func_file);
end

% Simple placeholder function - for basic extraction without SPM
function ts = extract_timeseries_basic(func_file, roi_coord, roi_radius)
    % Placeholder - This will be implemented based on your available tools
    ts = zeros(1, 100); % Dummy data
    fprintf('  Extracting time series using basic method from %s\n', func_file);
end
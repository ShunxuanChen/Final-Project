% compute_functional_connectivity.m
% Function to compute functional connectivity from ROI time series

function compute_functional_connectivity(timeseries_file, output_dir)
    % Check inputs
    if ~exist(timeseries_file, 'file')
        error('Time series file not found: %s', timeseries_file);
    end
    
    % Create output directory if needed
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
        fprintf('Created output directory: %s\n', output_dir);
    end
    
    % Load time series data
    fprintf('Loading time series data...\n');
    load(timeseries_file);
    
    % Create structure to store results
    fc_results = struct();
    
    % Process control group
    fprintf('Computing functional connectivity for control group...\n');
    control_subjects = fieldnames(all_data.control);
    
    for i = 1:length(control_subjects)
        subject = control_subjects{i};
        fprintf('  Processing subject %s (%d of %d)...\n', subject, i, length(control_subjects));
        
        % Create subject results structure
        fc_results.control.(subject) = struct();
        
        % Process music condition
        fc_results.control.(subject).music = struct();
        
        if isfield(all_data.control.(subject), 'music')
            music_runs = fieldnames(all_data.control.(subject).music);
            
            for j = 1:length(music_runs)
                run = music_runs{j};
                fc_results.control.(subject).music.(run) = struct();
                
                % Get ROIs
                roi_names = fieldnames(all_data.control.(subject).music.(run));
                
                % Compute FC matrix (correlation between ROIs)
                n_rois = length(roi_names);
                fc_matrix = zeros(n_rois, n_rois);
                
                for roi1 = 1:n_rois
                    for roi2 = 1:n_rois
                        ts1 = all_data.control.(subject).music.(run).(roi_names{roi1});
                        ts2 = all_data.control.(subject).music.(run).(roi_names{roi2});
                        
                        % Ensure time series are vectors
                        ts1 = reshape(ts1, 1, numel(ts1));
                        ts2 = reshape(ts2, 1, numel(ts2));
                        
                        % Compute correlation
                        try
                            c = corrcoef(ts1, ts2);
                            fc_matrix(roi1, roi2) = c(1, 2);
                        catch
                            % If corrcoef fails, use a manual calculation
                            try
                                % Make sure both time series have the same length
                                min_length = min(length(ts1), length(ts2));
                                ts1 = ts1(1:min_length);
                                ts2 = ts2(1:min_length);
                                
                                % Normalize both time series
                                ts1_norm = (ts1 - mean(ts1)) / std(ts1);
                                ts2_norm = (ts2 - mean(ts2)) / std(ts2);
                                
                                % Calculate correlation
                                fc_matrix(roi1, roi2) = sum(ts1_norm .* ts2_norm) / min_length;
                            catch
                                % If everything fails, just use zero
                                fc_matrix(roi1, roi2) = 0;
                                fprintf('    Warning: Could not compute correlation for ROIs %s and %s\n', ...
                                    roi_names{roi1}, roi_names{roi2});
                            end
                        end
                    end
                end
                
                % Store FC matrix
                fc_results.control.(subject).music.(run).fc_matrix = fc_matrix;
                fc_results.control.(subject).music.(run).roi_names = roi_names;
            end
        end
        
        % Process non-music condition
        fc_results.control.(subject).nonmusic = struct();
        
        if isfield(all_data.control.(subject), 'nonmusic')
            nonmusic_runs = fieldnames(all_data.control.(subject).nonmusic);
            
            for j = 1:length(nonmusic_runs)
                run = nonmusic_runs{j};
                fc_results.control.(subject).nonmusic.(run) = struct();
                
                % Get ROIs
                roi_names = fieldnames(all_data.control.(subject).nonmusic.(run));
                
                % Compute FC matrix (correlation between ROIs)
                n_rois = length(roi_names);
                fc_matrix = zeros(n_rois, n_rois);
                
                for roi1 = 1:n_rois
                    for roi2 = 1:n_rois
                        ts1 = all_data.control.(subject).nonmusic.(run).(roi_names{roi1});
                        ts2 = all_data.control.(subject).nonmusic.(run).(roi_names{roi2});
                        
                        % Ensure time series are vectors
                        ts1 = reshape(ts1, 1, numel(ts1));
                        ts2 = reshape(ts2, 1, numel(ts2));
                        
                        % Compute correlation
                        try
                            c = corrcoef(ts1, ts2);
                            fc_matrix(roi1, roi2) = c(1, 2);
                        catch
                            % If corrcoef fails, use a manual calculation
                            try
                                % Make sure both time series have the same length
                                min_length = min(length(ts1), length(ts2));
                                ts1 = ts1(1:min_length);
                                ts2 = ts2(1:min_length);
                                
                                % Normalize both time series
                                ts1_norm = (ts1 - mean(ts1)) / std(ts1);
                                ts2_norm = (ts2 - mean(ts2)) / std(ts2);
                                
                                % Calculate correlation
                                fc_matrix(roi1, roi2) = sum(ts1_norm .* ts2_norm) / min_length;
                            catch
                                % If everything fails, just use zero
                                fc_matrix(roi1, roi2) = 0;
                                fprintf('    Warning: Could not compute correlation for ROIs %s and %s\n', ...
                                    roi_names{roi1}, roi_names{roi2});
                            end
                        end
                    end
                end
                
                % Store FC matrix
                fc_results.control.(subject).nonmusic.(run).fc_matrix = fc_matrix;
                fc_results.control.(subject).nonmusic.(run).roi_names = roi_names;
            end
        end
    end
    
    % Process MDD group with similar logic
    fprintf('Computing functional connectivity for MDD group...\n');
    mdd_subjects = fieldnames(all_data.mdd);
    
    for i = 1:length(mdd_subjects)
        subject = mdd_subjects{i};
        fprintf('  Processing subject %s (%d of %d)...\n', subject, i, length(mdd_subjects));
        
        % Create subject results structure
        fc_results.mdd.(subject) = struct();
        
        % Process music condition
        fc_results.mdd.(subject).music = struct();
        
        if isfield(all_data.mdd.(subject), 'music')
            music_runs = fieldnames(all_data.mdd.(subject).music);
            
            for j = 1:length(music_runs)
                run = music_runs{j};
                fc_results.mdd.(subject).music.(run) = struct();
                
                % Get ROIs
                roi_names = fieldnames(all_data.mdd.(subject).music.(run));
                
                % Compute FC matrix (correlation between ROIs)
                n_rois = length(roi_names);
                fc_matrix = zeros(n_rois, n_rois);
                
                for roi1 = 1:n_rois
                    for roi2 = 1:n_rois
                        ts1 = all_data.mdd.(subject).music.(run).(roi_names{roi1});
                        ts2 = all_data.mdd.(subject).music.(run).(roi_names{roi2});
                        
                        % Ensure time series are vectors
                        ts1 = reshape(ts1, 1, numel(ts1));
                        ts2 = reshape(ts2, 1, numel(ts2));
                        
                        % Compute correlation
                        try
                            c = corrcoef(ts1, ts2);
                            fc_matrix(roi1, roi2) = c(1, 2);
                        catch
                            % If corrcoef fails, use a manual calculation
                            try
                                % Make sure both time series have the same length
                                min_length = min(length(ts1), length(ts2));
                                ts1 = ts1(1:min_length);
                                ts2 = ts2(1:min_length);
                                
                                % Normalize both time series
                                ts1_norm = (ts1 - mean(ts1)) / std(ts1);
                                ts2_norm = (ts2 - mean(ts2)) / std(ts2);
                                
                                % Calculate correlation
                                fc_matrix(roi1, roi2) = sum(ts1_norm .* ts2_norm) / min_length;
                            catch
                                % If everything fails, just use zero
                                fc_matrix(roi1, roi2) = 0;
                                fprintf('    Warning: Could not compute correlation for ROIs %s and %s\n', ...
                                    roi_names{roi1}, roi_names{roi2});
                            end
                        end
                    end
                end
                
                % Store FC matrix
                fc_results.mdd.(subject).music.(run).fc_matrix = fc_matrix;
                fc_results.mdd.(subject).music.(run).roi_names = roi_names;
            end
        end
        
        % Process non-music condition
        fc_results.mdd.(subject).nonmusic = struct();
        
        if isfield(all_data.mdd.(subject), 'nonmusic')
            nonmusic_runs = fieldnames(all_data.mdd.(subject).nonmusic);
            
            for j = 1:length(nonmusic_runs)
                run = nonmusic_runs{j};
                fc_results.mdd.(subject).nonmusic.(run) = struct();
                
                % Get ROIs
                roi_names = fieldnames(all_data.mdd.(subject).nonmusic.(run));
                
                % Compute FC matrix (correlation between ROIs)
                n_rois = length(roi_names);
                fc_matrix = zeros(n_rois, n_rois);
                
                for roi1 = 1:n_rois
                    for roi2 = 1:n_rois
                        ts1 = all_data.mdd.(subject).nonmusic.(run).(roi_names{roi1});
                        ts2 = all_data.mdd.(subject).nonmusic.(run).(roi_names{roi2});
                        
                        % Ensure time series are vectors
                        ts1 = reshape(ts1, 1, numel(ts1));
                        ts2 = reshape(ts2, 1, numel(ts2));
                        
                        % Compute correlation
                        try
                            c = corrcoef(ts1, ts2);
                            fc_matrix(roi1, roi2) = c(1, 2);
                        catch
                            % If corrcoef fails, use a manual calculation
                            try
                                % Make sure both time series have the same length
                                min_length = min(length(ts1), length(ts2));
                                ts1 = ts1(1:min_length);
                                ts2 = ts2(1:min_length);
                                
                                % Normalize both time series
                                ts1_norm = (ts1 - mean(ts1)) / std(ts1);
                                ts2_norm = (ts2 - mean(ts2)) / std(ts2);
                                
                                % Calculate correlation
                                fc_matrix(roi1, roi2) = sum(ts1_norm .* ts2_norm) / min_length;
                            catch
                                % If everything fails, just use zero
                                fc_matrix(roi1, roi2) = 0;
                                fprintf('    Warning: Could not compute correlation for ROIs %s and %s\n', ...
                                    roi_names{roi1}, roi_names{roi2});
                            end
                        end
                    end
                end
                
                % Store FC matrix
                fc_results.mdd.(subject).nonmusic.(run).fc_matrix = fc_matrix;
                fc_results.mdd.(subject).nonmusic.(run).roi_names = roi_names;
            end
        end
    end
    
    % Save results
    save(fullfile(output_dir, 'fc_results.mat'), 'fc_results');
    fprintf('Functional connectivity analysis complete. Results saved to %s\n', ...
        fullfile(output_dir, 'fc_results.mat'));
end
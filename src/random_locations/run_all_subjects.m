function run_all_subjects(subjects, data_dir, n_runs)

spm('defaults','eeg');
rng('shuffle');

for idx=1:n_runs
    
    % Run localization and sliding time window inversion
    for i=1:length(subjects)
        invert_burst_tc_subject(subjects(i),idx,'base_dir',...
            fullfile('../../output/data',data_dir));
    end
end
function run_all_subjects(subjects, data_dir, n_runs)

spm('defaults','eeg');
rng('shuffle');

% Run preprocessing
for run_idx=1:n_runs
    for i=1:length(subjects)
        subj_id=subjects(i).subj_id;
        % Shuffle times
        preprocess(fullfile('..\..\output\data',data_dir,subj_id));

        % Run localization and sliding time window inversion
        invert_burst_tc_subject(subjects(i),run_idx,'base_dir',...
            fullfile('../../output/data',data_dir));
    end
end

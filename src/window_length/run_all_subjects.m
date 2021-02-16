function run_all_subjects(subjects, data_dir)

spm('defaults','eeg');

window_size=[2 5 15];

% Run localization and sliding time window inversion
for j=1:length(window_size)
    for i=1:length(subjects)
        invert_burst_tc_subject(subjects(i),'base_dir',...
            fullfile('../../output/data',data_dir), 'win_size', window_size(j));
    end
end
function run_all_subjects(subjects, data_dir)

spm('defaults','eeg');
thresholds=[.7 .75 .85 ];
for j=1:length(thresholds)
    % Run localization and sliding time window inversion
    for i=1:length(subjects)
        invert_burst_tc_subject(subjects(i),thresholds(j), 'base_dir',...
            fullfile('../../output/data',data_dir));
    end
end
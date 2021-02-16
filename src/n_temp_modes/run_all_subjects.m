function run_all_subjects(subjects, data_dir)

spm('defaults','eeg');
n_temporal_modes=[2 6];

for j=1:length(n_temporal_modes)
    % Run localization and sliding time window inversion
    for i=1:length(subjects)
       invert_burst_tc_subject(subjects(i),'base_dir',...
            fullfile('../../output/data',data_dir), 'n_temp_modes', n_temporal_modes(j));
    end
end

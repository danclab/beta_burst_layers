function run_all_subjects(subjects, data_dir)

spm('defaults','eeg');
patch_size=[2.5 7.5 10];

for j=1:length(patch_size)
    % Run localization and sliding time window inversion
    for i=1:length(subjects)
       invert_burst_tc_subject(subjects(i),'base_dir',...
            fullfile('../../output/data',data_dir), 'patch_size', patch_size(j));
    end
end

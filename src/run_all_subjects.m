function run_all_subjects(subjects, data_dir)

spm('defaults','eeg');

addpath('./preprocessing');
% Run preprocessing
for i=1:length(subjects)
    subj_id=subjects(i).subj_id;
    preprocess(fullfile('..\data',data_dir,subj_id),...
        fullfile('..\output\data',data_dir,subj_id));
end
rmpath('./preprocessing');

% Run localization and sliding time window inversion
for i=1:length(subjects)
    invert_burst_tc_subject(subjects(i),'base_dir',fullfile('../output/data',data_dir));
end

% Invert combined surface
for i=1:length(subjects)
    subj_dir=fullfile('../output/data',data_dir,subjects(i).subj_id);
    % Preprocessed data file
    data_file=fullfile(subj_dir, 'mrcresp_TafdfC.mat');
    coreg_fname=fullfile(subj_dir, 'grey_mrcresp_TafdfC.mat');
    mesh_fname=fullfile('../data/surf/',...
        sprintf('%s-synth',subjects(i).subj_id),'surf',...
        'white.ds-pial.ds.link_vector.gii');
    invert_ebb(subjects(i),data_file,coreg_fname,mesh_fname,'surf_dir',...
        '../data/surf','mri_dir', '../data/mri','patch_size',5,...
        'n_temp_modes',4);
end

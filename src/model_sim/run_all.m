function run_all(subj_info, data_dir, type)

rng('default');
rng('shuffle');
spm('defaults','eeg');

g=gifti(fullfile('../../output/data', data_dir, subj_info.subj_id,'pial_mrcresp_TafdfC_1_t-2000_-1800_f_1.gii'));
[max_val,pial_vertex]=max(g.cdata(:));

for idx=1:10
    simulate_model_bursts(subj_info, pial_vertex, idx, type,...
        'base_dir', fullfile('../../output/data',data_dir),...
        'surf_dir', '../../data/surf', 'mri_dir', '../../data/mri');

    %Run localization and sliding time window inversion
    invert_burst_tc_subject(subj_info,idx,type,'base_dir',...
        fullfile('../../output/data',data_dir));

    addpath('..');
    % Invert combined surface
    subj_dir=fullfile('../../output/data',data_dir,subj_info.subj_id,'model_sim');
    % Preprocessed data file
    data_file=fullfile(subj_dir, sprintf('msim_%s_%d_grey_rcresp_TafdfC.mat',type,idx));
    coreg_fname=fullfile(subj_dir, sprintf('grey_msim_%s_%d_grey_rcresp_TafdfC.mat',type,idx));
    mesh_fname=fullfile('../../data/surf/',...
        sprintf('%s-synth',subj_info.subj_id),'surf',...
        'white.ds-pial.ds.link_vector.gii');
    invert_ebb(subj_info,data_file,coreg_fname,mesh_fname,'surf_dir',...
        '../../data/surf','mri_dir', '../../data/mri','patch_size',5,...
        'n_temp_modes',4);
    rmpath('..');
end

% Plot results
plot_subject_tc(subj_info, type, 'base_dir',fullfile('../../output/data',data_dir));
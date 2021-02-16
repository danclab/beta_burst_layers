function run_all(subj_info, data_dir, type)

rng('default');
rng('shuffle');
spm('defaults','eeg');

g=gifti(fullfile('../../../output/data', data_dir, subj_info.subj_id,'pial_mrcresp_TafdfC_1_t-2000_-1800_f_1.gii'));
[max_val,pial_vertex]=max(g.cdata(:));

coregerrs=[0 .25 .5 .75 1 1.5 2 3 4 5];

orig_fid=[subj_info.nas; subj_info.lpa; subj_info.rpa];
meanfid=mean(orig_fid);
zeromeanfid=[orig_fid-repmat(meanfid,3,1) ones(3,1)];

for c_idx=1:length(coregerrs)
    coregerr=coregerrs(c_idx);
    subj_dir=fullfile('../../../output/data',data_dir,subj_info.subj_id,'simple_sim',sprintf('coregerr_%d',coregerr));
    mkdir(subj_dir);
    for idx=1:50
        
        % Translation vector
        translation=coregerr*.1;
        shift_vec=randn(1,3);
        shift_vec=shift_vec./sqrt(sum(shift_vec.^2)).*translation;
        
        % Rotation vector
        rotation_rad=coregerr*pi/180.0;
        rot_vec=randn(1,3);
        rot_vec=rot_vec./sqrt(sum(rot_vec.^2)).*rotation_rad;
        
        % Apply transformation to fiducial locations
        P=[shift_vec rot_vec];
        [A]=spm_matrix(P);
        fid=(A*zeromeanfid')';
        fid=fid(:,1:3)+repmat(meanfid,3,1);        
        
        subj_info.nas=orig_fid(1,:);
        subj_info.lpa=orig_fid(2,:);
        subj_info.rpa=orig_fid(3,:);
        
        simulate_simple_bursts(subj_info, pial_vertex, idx, type, coregerr,...
            'base_dir', fullfile('../../../output/data',data_dir),...
            'surf_dir', '../../../data/surf', 'mri_dir', '../../../data/mri');

        subj_info.nas=fid(1,:);
        subj_info.lpa=fid(2,:);
        subj_info.rpa=fid(3,:);
        
        %Run localization and sliding time window inversion
        invert_burst_tc_subject(subj_info,idx,type,coregerr,'base_dir',...
            fullfile('../../../output/data',data_dir));

%         % Invert combined surface
%         addpath('../..');
%         % Preprocessed data file
%         data_file=fullfile(subj_dir, sprintf('msim_%s_%d_grey_rcresp_TafdfC.mat',type,idx));
%         coreg_fname=fullfile(subj_dir, sprintf('grey_msim_%s_%d_grey_rcresp_TafdfC.mat',type,idx));
%         mesh_fname=fullfile('../../../data/surf/',...
%             sprintf('%s-synth',subj_info.subj_id),'surf',...
%             'white.ds-pial.ds.link_vector.gii');
%         invert_ebb(subj_info,data_file,coreg_fname,mesh_fname,'surf_dir',...
%             '../../../data/surf','mri_dir', '../../../data/mri','patch_size',5,...
%             'n_temp_modes',4);
        rmpath('../..');
        
        delete(fullfile(subj_dir,'msim_standard*'));
        delete(fullfile(subj_dir,'sim_standard*'));
        delete(fullfile(subj_dir,'SPMgainmatrix_*'));
        delete(fullfile(subj_dir,'pial*'));
        delete(fullfile(subj_dir,'white*'));
    end    
end
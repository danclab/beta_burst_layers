function invert_burst_tc_subject(subj_info, sim_idx, type, coregerr, varargin)

defaults = struct('base_dir', '../../../data/JB_BUTTON_LOCKED_d3_ers',...
    'surf_dir', '../../../data/surf', 'mri_dir', '../../../data/mri', 'patch_size', 5,...
    'win_size', 10, 'win_overlap',true,'n_temp_modes', 4,...
    'data_type', 'mean_evoked', 'plot_dir', '');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');
spm_jobman('initcfg');

base_dir_parts=strsplit(params.base_dir,filesep);
data_dir=fullfile('../../../output/data',base_dir_parts{end},subj_info.subj_id,'simple_sim',sprintf('coregerr_%d',coregerr));

% Where to save plots
if length(params.plot_dir)==0
    params.plot_dir=fullfile('../../../output/figures',base_dir_parts{end},...
        subj_info.subj_id,params.data_type,...
        sprintf('%d_temp_modes', params.n_temp_modes),'simple_sim',sprintf('coregerr_%d',coregerr));
end
if exist(params.plot_dir,'dir')~=7
    mkdir(params.plot_dir);
end

pial_clusters=invert_burst_subject(subj_info,...
    sim_idx, type, coregerr, 'base_dir',params.base_dir,'surf_dir', params.surf_dir ,...
    'mri_dir', params.mri_dir,'patch_size', params.patch_size,...
    'data_type', params.data_type, 'plot_dir', params.plot_dir, 'plot', false);

fname=sprintf('msim_%s_%d_grey_rcresp_TafdfC.mat',type,sim_idx);
if strcmp(params.data_type,'evoked')
    fname=sprintf('sim_%s_%d_grey_rcresp_TafdfC.mat',type,sim_idx);
end
data_file=fullfile(data_dir, fname);
D=spm_eeg_load(data_file);
times=D.time;
zero_time=times(round(length(times)/2))+.5*(times(round(length(times)/2+1))-times(round(length(times)/2)));

subj_surf_dir=fullfile(params.surf_dir, sprintf('%s-synth',...
    subj_info.subj_id),'surf');
pial_fname=fullfile(subj_surf_dir,'pial.ds.link_vector.gii');
pial_surf=gifti(pial_fname);
wm_fname=fullfile(subj_surf_dir,'white.ds.link_vector.gii');

% Setup mesh lists
mesh_fnames={wm_fname, pial_fname};
mesh_names={'white','pial'};

for cluster_idx=1:length(pial_clusters)
    cluster=pial_clusters(cluster_idx);
    max_vert=cluster.vertices(cluster.max_idx);
    
    norm_diffs=[];
    for i=1:length(cluster.vertices)
        norm_diffs(i)=atan2(norm(cross(pial_surf.normals(max_vert,:),pial_surf.normals(cluster.vertices(i),:))),...
            dot(pial_surf.normals(max_vert,:),pial_surf.normals(cluster.vertices(i),:)));
    end
    cluster.inv_verts=cluster.vertices(find(abs(norm_diffs<.1)));

    cluster.tc_fvals=zeros(length(cluster.inv_verts),length(mesh_names),length(times));
    for v_idx=1:length(cluster.inv_verts)
        pial_max_vert=cluster.inv_verts(v_idx);

        priors={pial_max_vert, pial_max_vert};
        recompute_lgain=(v_idx==1);
        for m_idx=1:length(mesh_names)
            [cluster.tc_fvals(v_idx,m_idx,:),wois]=invert_sliding_window(subj_info,...
                priors{m_idx}, mesh_names{m_idx}, mesh_fnames{m_idx},...
                sim_idx, type, coregerr, 'base_dir', params.base_dir,'mri_dir', params.mri_dir,...
                'patch_size', 5, 'n_temp_modes', 4,...
                'data_type', params.data_type, 'win_size', 10,...
                'win_overlap',params.win_overlap,...
                'recompute_lgain',recompute_lgain);
        end
        centered_wois=wois-zero_time*1000;
        inner_times=centered_wois(:,1)+.5*(centered_wois(:,2)-centered_wois(:,1));
        if params.win_overlap
            left_idx=round((params.win_size-1)/2)+1;
            right_idx=length(times)-round((params.win_size-1)/2);
        else
            left_idx=1;
            right_idx=length(centered_wois);
        end
    end

    cluster.f_diff=squeeze(cluster.tc_fvals(:,2,left_idx:right_idx)-cluster.tc_fvals(:,1,left_idx:right_idx));
    cluster.f_diff=reshape(cluster.f_diff,[size(cluster.tc_fvals,1) (right_idx-left_idx+1)]); 
    new_pial_clusters(cluster_idx)=cluster;
end

invert_burst_tc_results=[];
invert_burst_tc_results.params=params;
invert_burst_tc_results.times=inner_times;
invert_burst_tc_results.left_idx=left_idx;
invert_burst_tc_results.right_idx=right_idx;
invert_burst_tc_results.clusters=new_pial_clusters;   

save(fullfile(data_dir, sprintf('invert_burst_tc_results_%s_%d.mat',type,sim_idx)), 'invert_burst_tc_results');

fig=figure();
subplot(2,1,1);
hold all
for cluster_idx=1:length(new_pial_clusters)
    plot(inner_times(left_idx:right_idx),mean(new_pial_clusters(cluster_idx).f_diff,1));
end
plot([inner_times(left_idx) inner_times(right_idx)],[0 0],'k');
plot([inner_times(left_idx) inner_times(right_idx)],[3 3],'k--');
plot([inner_times(left_idx) inner_times(right_idx)],[-3 -3],'k--');
xlim([inner_times(left_idx) inner_times(right_idx)]);
xlabel('Time (ms)')
ylabel('BC Fpial-Fwhite');
subplot(2,1,2);
hold all
for cluster_idx=1:length(new_pial_clusters)
    plot(inner_times(left_idx:right_idx),mean(new_pial_clusters(cluster_idx).f_diff,1));
end
plot([inner_times(left_idx) inner_times(right_idx)],[0 0],'k');
plot([inner_times(left_idx) inner_times(right_idx)],[3 3],'k--');
plot([inner_times(left_idx) inner_times(right_idx)],[-3 -3],'k--');
xlim([inner_times(left_idx) inner_times(right_idx)]);
xlabel('Time (ms)')
ylabel('Fpial-Fwhite');
% saveas(fig, fullfile(params.plot_dir, sprintf('pial-wm_%s_%d_free_energy_burst_tc_pial.png',type,sim_idx)), 'png');
% saveas(fig, fullfile(params.plot_dir, sprintf('pial-wm_%s_%d_free_energy_burst_tc_pial.eps',type,sim_idx)), 'eps');
% saveas(fig, fullfile(params.plot_dir, sprintf('pial-wm_%s_%d_free_energy_burst_tc_pial.fig',type,sim_idx)), 'fig');

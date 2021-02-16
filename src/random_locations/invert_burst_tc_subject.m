function invert_burst_tc_subject(subj_info, idx, varargin)

defaults = struct('base_dir', '../../data/JB_BUTTON_LOCKED_d3_ers',...
    'surf_dir', '../../data/surf', 'mri_dir', '../../data/mri', 'patch_size', 5,...
    'win_size', 10, 'win_overlap',true,'n_temp_modes', 4,...
    'data_type', 'mean_evoked');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');
spm_jobman('initcfg');

base_dir_parts=strsplit(params.base_dir,filesep);
data_dir=fullfile('../../output/data',base_dir_parts{end},subj_info.subj_id);

% Load localizer results
load(fullfile(data_dir, 'invert_burst_subject_results.mat'));

% Meshes to use
subj_surf_dir=fullfile(params.surf_dir, sprintf('%s-synth',...
    subj_info.subj_id),'surf');
pial_fname=fullfile(subj_surf_dir,'pial.ds.link_vector.gii');
pial_surf=gifti(pial_fname);
pial_inflated_fname=fullfile(subj_surf_dir, 'pial.ds.inflated.gii');
pial_inflated_surf=gifti(pial_inflated_fname);
wm_fname=fullfile(subj_surf_dir,'white.ds.link_vector.gii');

% For each cluster found by the localizer (so we have same number of random
% clusters
for cluster_idx=1:length(invert_burst_subject_results.clusters)
    % Pick random vertex
    v=randi([1,size(pial_surf.vertices,1)]);
    % Find vertices within 10mm
    pial_dist=sqrt(sum((pial_inflated_surf.vertices-repmat(pial_inflated_surf.vertices(v,:),size(pial_surf.vertices,1),1)).^2,2));    
    vs=find((pial_dist<=10) & (pial_dist>0));
    % Create new cluster
    pial_clusters(cluster_idx).vertices=[v; vs];
    pial_clusters(cluster_idx).max_idx=1;
    pial_clusters(cluster_idx).coords=pial_surf.vertices([v; vs],:);
end

% Load averaged burst file
fname='mrcresp_TafdfC.mat';
if strcmp(params.data_type,'evoked')
    fname='rcresp_TafdfC.mat';
end
data_file=fullfile(data_dir, fname);
D=spm_eeg_load(data_file);

% Get times and zero time
times=D.time;
zero_time=times(round(length(times)/2))+.5*(times(round(length(times)/2+1))-times(round(length(times)/2)));

% Setup mesh lists
mesh_fnames={wm_fname, pial_fname};
mesh_names={'white','pial'};

% For each random cluster
for cluster_idx=1:length(pial_clusters)
    cluster=pial_clusters(cluster_idx);
    max_vert=cluster.vertices(cluster.max_idx);
    
    % Difference in link vector angle at max vertex of cluster and all
    % other vertices in cluster
    norm_diffs=[];
    for i=1:length(cluster.vertices)
        norm_diffs(i)=atan2(norm(cross(pial_surf.normals(max_vert,:),pial_surf.normals(cluster.vertices(i),:))),...
            dot(pial_surf.normals(max_vert,:),pial_surf.normals(cluster.vertices(i),:)));
    end
    % Invert at each vertex with link vector angle within .1 radians of the
    % max vertex
    cluster.inv_verts=cluster.vertices(find(abs(norm_diffs<.1)));

    % Run sliding window inversion at each selected vertex
    cluster.tc_fvals=zeros(length(cluster.inv_verts),length(mesh_names),length(times));
    for v_idx=1:length(cluster.inv_verts)
        pial_prior=cluster.inv_verts(v_idx);

        % Because of our cool downsampling method, the vertices in the
        % white matter and pial surface are cooresponding
        priors={pial_prior, pial_prior};
        % Only have to recompute the lead field gain matrix once
        recompute_lgain=(cluster_idx==1) & (v_idx==1);
        
        % Run sliding window inversion for each mesh
        for m_idx=1:length(mesh_names)
            [cluster.tc_fvals(v_idx,m_idx,:),wois]=invert_sliding_window(subj_info,...
                priors{m_idx}, mesh_names{m_idx}, mesh_fnames{m_idx},...
                'base_dir', params.base_dir,'mri_dir', params.mri_dir,...
                'patch_size', 5, 'n_temp_modes', 4,...
                'data_type', params.data_type, 'win_size', 10,...
                'win_overlap',params.win_overlap,...
                'recompute_lgain',recompute_lgain);
        end
        % Create centered WOI time
        centered_wois=wois-zero_time*1000;
        % Figure out center of each sliding time window
        inner_times=centered_wois(:,1)+.5*(centered_wois(:,2)-centered_wois(:,1));
        % If overlapping windows, cut off the edges
        if params.win_overlap
            left_idx=round((params.win_size-1)/2)+1;
            right_idx=length(times)-round((params.win_size-1)/2);
        else
            left_idx=1;
            right_idx=length(centered_wois);
        end
    end

    % Compute free energy difference
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

save(fullfile(data_dir, sprintf('invert_burst_tc_results-random_locations_%d.mat',idx)), 'invert_burst_tc_results');

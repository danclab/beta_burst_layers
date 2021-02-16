function clusters=invert_burst_subject(subj_info, sim_idx, type, varargin)

defaults = struct('base_dir','../../data/JB_BUTTON_LOCKED_d3_ers',...
    'surf_dir', '../../data/surf', 'mri_dir', '../../data/mri',...
    'patch_size', 5, 'plot', true, 'n_temp_modes', 4,...
    'data_type', 'mean_evoked', 'plot_dir', '');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Threshold percentage
percent_thresh=.80;

% Where to put output data
base_dir_parts=strsplit(params.base_dir,filesep);
data_dir=fullfile('../../output/data',base_dir_parts{end},subj_info.subj_id,'simple_sim');
mkdir(data_dir);

% Data file to load
data_file=fullfile(data_dir, sprintf('msim_%s_%d_grey_rcresp_TafdfC.mat', type, sim_idx));
if strcmp(params.data_type,'evoked')
    data_file=fullfile(data_dir, sprintf('sim_%s_%d_grey_rcresp_TafdfC.mat', type, sim_idx));
end

% Where to save plots
if length(params.plot_dir)==0
    params.plot_dir=fullfile('../../output/figures',base_dir_parts{end},...
        subj_info.subj_id,params.data_type,...
        sprintf('%d_temp_modes', params.n_temp_modes),'simple_sim');
end
if exist(params.plot_dir,'dir')~=7
    mkdir(params.plot_dir);
end
addpath('..');
% Subject surfaces
subj_surf_dir=fullfile(params.surf_dir, sprintf('%s-synth',...
    subj_info.subj_id),'surf');
pial_original_fname=fullfile(subj_surf_dir,'pial.gii');
pial_fname=fullfile(subj_surf_dir,'pial.ds.link_vector.gii');
pial_inflated_fname=fullfile(subj_surf_dir, 'pial.ds.inflated.gii');
wm_fname=fullfile(subj_surf_dir,'white.ds.link_vector.gii');
wm_inflated_fname=fullfile(subj_surf_dir, 'white.ds.inflated.gii');

norm=compute_surface_normals(params.surf_dir, sprintf('%s-synth',subj_info.subj_id), 'pial', 'link_vector');
pial=gifti(pial_fname);
pial.normals=norm;
save(pial,pial_fname);
 
pial_inflated=gifti(pial_inflated_fname);    

% Which vertices belong to which hemisphere
pial_hemisphere_map=get_hemisphere_map(pial_fname, pial_original_fname,'recompute',false);

% Get LH motor cortex mask
[anat_pial_mask,anat_wm_mask]=get_anatomical_mask('lh_motor', pial_fname,...
    pial_inflated_fname, wm_fname, wm_inflated_fname,...
    subj_info.coords('lh_motor'), 50, 'recompute', false);

spm('defaults','eeg');
spm_jobman('initcfg');

% Create smoothed meshes
[smoothkern]=spm_eeg_smoothmesh_mm(pial_fname, params.patch_size);

% Coregistered filename
coreg_fname=fullfile(data_dir, sprintf('pial_msim_%s_%d_grey_rcresp_TafdfC.mat', type, sim_idx));
    
invert_ebb(subj_info, data_file, coreg_fname, pial_fname,...
    'surf_dir', params.surf_dir, 'mri_dir', params.mri_dir,...
    'patch_size', params.patch_size, 'n_temp_modes', params.n_temp_modes);        

% Load mesh results
mesh_results=gifti(fullfile(data_dir, sprintf('pial_msim_%s_%d_grey_rcresp_TafdfC_1_t-2000_-1800_f_1.gii',type,sim_idx)));
    
% Create mask
mask=find(pial_hemisphere_map==1);
mask=intersect(mask, anat_pial_mask);
        
% Find vertices in mask where results within percentile threshold
cluster_mask=intersect(mask,find(mesh_results.cdata(:)>=max(mesh_results.cdata(mask))*percent_thresh));

cluster_mesh=struct();
% Get vertices and all faces containing a vertex from the ROI
cluster_mesh.vertices=pial.vertices;
% Find faces with vertices in p_mask
[rows,cols]=find(ismember(pial.faces,cluster_mask));
% Only include faces with more than one vertex in mask
cluster_mesh.faces=pial.faces(unique(rows),:);
% Split into clusters
cluster_fv=splitFV(cluster_mesh);
    
clusters=[];
for c_idx=1:length(cluster_fv)
    vs=find(ismember(pial.vertices,cluster_fv(c_idx).vertices,'rows')>0);
    clusters(c_idx).vertices=vs;
    clusters(c_idx).max_idx=find(mesh_results.cdata(vs)==max(mesh_results.cdata(vs)));
    clusters(c_idx).coords=pial.vertices(vs,:);
end
    
if params.plot    
    data=zeros(size(mesh_results.cdata));
    nclusters=length(clusters);
    for cluster_idx=1:length(clusters)
        cluster=clusters(cluster_idx);
        data(cluster.vertices)=cluster_idx*1.0/nclusters;                
    end
    [ax,metric_data]=plot_surface_metric(pial_inflated, data,...
        'clip_vals',false, 'threshold', 0.1, 'limits', [0 1], 'custom_cm', false,...
        'specular_strength', 0.0, 'ambient_strength', 0.8, 'face_lighting', '');
    set(ax,'CameraViewAngle',6.028);
    set(ax,'CameraUpVector',subj_info.camera_up_vector('topdown'));
    set(ax,'CameraPosition',subj_info.camera_position('topdown'));
    fig=get(ax,'Parent');
    saveas(fig, fullfile(params.plot_dir, sprintf('pial_%s_%d_evoked_burst_t-Inf_Inf_clusters.png',type,sim_idx)), 'png');

    [ax,metric_data]=plot_surface_metric(pial_inflated,mesh_results.cdata(:),...
        'clip_vals',false,'custom_cm',false);
    set(ax,'CameraViewAngle',6.028);
    set(ax,'CameraUpVector',subj_info.camera_up_vector('topdown'));
    set(ax,'CameraPosition',subj_info.camera_position('topdown'));
    fig=get(ax,'Parent');
    saveas(fig, fullfile(params.plot_dir, sprintf('pial_%s_%d_evoked_burst_t-Inf_Inf.png',type,sim_idx)), 'png');

    [ax,metric_data]=plot_surface_metric(pial_inflated,mesh_results.cdata(:),...
        'clip_vals',false,'custom_cm',false, 'mask', cluster_mask);
    set(ax,'CameraViewAngle',6.028);
    set(ax,'CameraUpVector',subj_info.camera_up_vector('topdown'));
    set(ax,'CameraPosition',subj_info.camera_position('topdown'));
    fig=get(ax,'Parent');
    saveas(fig, fullfile(params.plot_dir, sprintf('pial_%s_%d_evoked_burst_t-Inf_Inf_mask.png',type,sim_idx)), 'png');
end       

invert_burst_subject_results=[];
invert_burst_subject_results.params=params;
invert_burst_subject_results.percent_thresh=percent_thresh;
invert_burst_subject_results.data_file=data_file;
invert_burst_subject_results.clusters=clusters;
save(fullfile(data_dir, sprintf('invert_burst_subject_results_%s_%d.mat',type,sim_idx)), 'invert_burst_subject_results');

rmpath('..');

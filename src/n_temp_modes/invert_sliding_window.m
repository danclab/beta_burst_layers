function [f_vals,wois]=invert_sliding_window(subj_info, prior, mesh_name, ...
    mesh_fname, varargin)

defaults = struct('base_dir', '../../data/JB_BUTTON_LOCKED_d3_ers',...
    'mri_dir', '../../data/mri', 'patch_size', 5, 'n_temp_modes', 1,...
    'data_type', 'mean_evoked', 'win_size', 10, 'win_overlap',true,...
    'recompute_lgain',true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Where to put output data
base_dir_parts=strsplit(params.base_dir,filesep);
data_dir=fullfile('../../output/data',base_dir_parts{end},subj_info.subj_id);

% Data file to load
fname='mrcresp_TafdfC.mat';
if strcmp(params.data_type,'evoked')
    fname='rcresp_TafdfC.mat';
end
data_file=fullfile(data_dir, fname);
D=spm_eeg_load(data_file);


% Create wois
times=D.time;
wois=[];
if params.win_overlap
    for t_idx=1:length(times)
        win_l=max(1,ceil(t_idx-params.win_size/2));
        win_r=min(length(times),floor(t_idx+params.win_size/2));
        woi=[times(win_l) times(win_r)].*1000;
        wois(t_idx,:)=woi;
    end
else
    ts=linspace(times(1),times(end),(times(end)-times(1))./(params.win_size/1000)).*1000;
    wois=[];
    for i=2:length(ts)
        wois(end+1,:)=[ts(i-1) ts(i)];
    end
end

spm('defaults','eeg');
spm_jobman('initcfg');

% Create smoothed meshes
[smoothkern]=spm_eeg_smoothmesh_mm(mesh_fname,params.patch_size);

% Coregistered filename
coreg_fname=fullfile(data_dir, sprintf('%s_mrcresp_TafdfC.mat',mesh_name));

% Coregister to mesh
clear jobs
matlabbatch={};
batch_idx=1;

% Copy datafile
matlabbatch{batch_idx}.spm.meeg.other.copy.D = {data_file};
matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = coreg_fname;
batch_idx=batch_idx+1;

% Coregister simulated dataset to reconstruction mesh
matlabbatch{batch_idx}.spm.meeg.source.headmodel.D = {coreg_fname};
matlabbatch{batch_idx}.spm.meeg.source.headmodel.val = 1;
matlabbatch{batch_idx}.spm.meeg.source.headmodel.comment = '';
matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir, subj_info.subj_id, [subj_info.headcast_t1 ',1'])};
matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {mesh_fname};
matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = subj_info.nas;
matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = subj_info.lpa;
matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = subj_info.rpa;
matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{batch_idx}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{batch_idx}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';            
spm_jobman('run', matlabbatch);    

if ~params.recompute_lgain
    D=spm_eeg_load(coreg_fname);
    [path,fname,ext]=fileparts(coreg_fname);
    D.inv{1}.gainmat = ['SPMgainmatrix_' fname '_1.mat'];
    save(D);
end

% Setup spatial modes for cross validation
spatialmodesname=fullfile(data_dir, sprintf('%s_testmodes.mat',mesh_name));    
[spatialmodesname,Nmodes,pctest]=spm_eeg_inv_prep_modes_xval(coreg_fname, [], spatialmodesname, 1, 0);

% so use all vertices that will be simulated on (plus a few more) as MSP priors
Ip=[prior];
% Save priors
patchfilename=fullfile(data_dir, 'patch.mat');
save(patchfilename,'Ip');

clear jobs
matlabbatch={};
batch_idx=1;

% Source reconstruction
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.D = {coreg_fname};
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.val = 1;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.whatconditions.all = 1;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.invfunc = 'Classic';
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.invtype = 'MSP'; %;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.wois = wois;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.foi = [1 256];
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.hanning = 1;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.isfixedpatch.fixedpatch.fixedfile = {patchfilename}; % '<UNDEFINED>';
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.isfixedpatch.fixedpatch.fixedrows = [1 Inf]; %'<UNDEFINED>';
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.patchfwhm = -params.patch_size; %% NB A fiddle here- need to properly quantify
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.mselect = 0;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.nsmodes = Nmodes;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.umodes = {spatialmodesname};
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.ntmodes = params.n_temp_modes;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.priors.priorsmask = {''};
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.priors.space = 0;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.restrict.locs = zeros(0, 3);
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.restrict.radius = 32;
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.isstandard.custom.outinv = '';
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.modality = {'All'};
matlabbatch{batch_idx}.spm.meeg.source.invertiter_slidingwindow.crossval = [pctest 1];                                
batch_idx=batch_idx+1;

[a,b]=spm_jobman('run', matlabbatch);
% Get F-values for inversion
Drecon=spm_eeg_load(a{1}.D{1});                
f_vals=Drecon.inv{1}.inverse.crossF;                  

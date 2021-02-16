function simulate_simple_bursts(subj_info, pial_vertex, sim_idx, type, coregerr, varargin)

defaults = struct('base_dir','../../../output/data/JB_BUTTON_LOCKED_d3_ers',...
    'surf_dir', '../../../data/surf', 'mri_dir', '../../../data/mri');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults', 'EEG');
spm_jobman('initcfg'); 

% Setup directories and files
data_dir=fullfile(params.base_dir,subj_info.subj_id);
data_file=fullfile(data_dir, 'rcresp_TafdfC.mat');
output_data_dir=fullfile(params.base_dir,subj_info.subj_id,'simple_sim',sprintf('coregerr_%d',coregerr));
mkdir(output_data_dir);

% Load meshes
subj_surf_dir=fullfile(params.surf_dir, sprintf('%s-synth',...
    subj_info.subj_id),'surf');
pial_fname=fullfile(subj_surf_dir,'pial.ds.link_vector.gii');
white_fname=fullfile(subj_surf_dir,'white.ds.link_vector.gii');
grey_fname=fullfile(subj_surf_dir,'white.ds-pial.ds.link_vector.gii');
    
white=gifti(white_fname);
pial=gifti(pial_fname);

% Convert pial and white matter vertices to gray matter mesh vertices
grey_pial_vertex=size(white.vertices,1)+pial_vertex;
grey_white_vertex=pial_vertex;

% Copy data file    
reg_file=fullfile(output_data_dir, 'grey_rcresp_TafdfC.mat');

%if sim_idx==1
    clear jobs
    matlabbatch=[];
    matlabbatch{1}.spm.meeg.other.copy.D = {data_file};
    matlabbatch{1}.spm.meeg.other.copy.outfile = reg_file;
    spm_jobman('run', matlabbatch);

    % Coregister to gray mesh
    matlabbatch=[];
    matlabbatch{1}.spm.meeg.source.headmodel.D = {reg_file};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,subj_info.subj_id, [subj_info.headcast_t1 ',1'])};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {grey_fname};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = subj_info.nas;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = subj_info.lpa;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = subj_info.rpa;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    spm_jobman('run', matlabbatch);

    % Setup spatial modes for cross validation
    spatialmodesname=fullfile(output_data_dir, 'testmodes.mat');    
    [spatialmodesname,Nmodes,pctest]=spm_eeg_inv_prep_modes_xval(reg_file, [], spatialmodesname, 1, 0);

    clear jobs
    matlabbatch={};
    batch_idx=1;

    % Run inversion (required for simulation)
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.D = {reg_file};
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.val = 1;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.whatconditions.all = 1;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invtype = 'EBB'; %;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.woi = [-Inf Inf];
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.foi = [0 256];
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.hanning = 0;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.npatches = 512;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.niter = 1;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = -5; %% NB A fiddle here- need to properly quantify
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = Nmodes;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.umodes = {spatialmodesname};
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = 4;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 0;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.modality = {'All'};
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.crossval = [pctest 1];      
    spm_jobman('run', matlabbatch);
%end

% Load inverted data file
D=spm_eeg_load(reg_file);   

% Create simulated signal
woi=[D.time(1)-(D.time(2)-D.time(1)) D.time(end)];
% time zero is midpoint of WOI
zero_time=D.time((length(D.time)-1)/2+1);

% Get dipole orientations
pial_unit_norm=pial.normals(pial_vertex,:);
white_unit_norm=-1.*white.normals(pial_vertex,:);

% Signal in deep layers
deep_signal_width=.025; % 25ms
deep_signal=exp(-((D.time-zero_time).^2)/(2*deep_signal_width^2));

% Signal in superficial layers
superficial_signal_width=.01; % 10ms
superficial_signal=exp(-((D.time-zero_time).^2)/(2*superficial_signal_width^2));

% Simulate source
if strcmp(type,'standard')
    % White wide - pial narrow
    % Coordinate of each signal
    simpos=[D.inv{1}.mesh.tess_mni.vert(grey_white_vertex,:)
        D.inv{1}.mesh.tess_mni.vert(grey_pial_vertex,:)];
    % Orientation of each dipole
    ormni=[white_unit_norm; pial_unit_norm];
    % Dipole momemnts
    nAmdipmom=[6 8];
    % Combined signal
    simsignal=[deep_signal; superficial_signal];
    % Width of patch
    dipfwhm=[5 5];
end

% Set sim signal to have unit variance
% simsignal=simsignal./repmat(std(simsignal,[],2),1,size(simsignal,2),1);
[D,meshsourceind]=spm_eeg_simulate({D}, sprintf('sim_%s_%d_',type,sim_idx), simpos,...
   simsignal, ormni, woi, [], -20, [], [], dipfwhm,...
   nAmdipmom, []);
save(D);

clear jobs
matlabbatch={};
batch_idx=1;    

% Average
matlabbatch{batch_idx}.spm.meeg.averaging.average.D = {fullfile(output_data_dir, sprintf('sim_%s_%d_grey_rcresp_TafdfC.mat',type,sim_idx))};
matlabbatch{batch_idx}.spm.meeg.averaging.average.userobust.standard = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.plv = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.prefix = fullfile(output_data_dir, 'm');
batch_idx=batch_idx+1;

spm_jobman('run', matlabbatch); 

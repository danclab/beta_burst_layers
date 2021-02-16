function simsignal=simulate_model_bursts(subj_info, pial_vertex, sim_idx, type, varargin)

defaults = struct('base_dir','../../output/data/JB_BUTTON_LOCKED_d3_ers',...
    'surf_dir', '../../data/surf', 'mri_dir', '../../data/mri');  %define default values
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
output_data_dir=fullfile(params.base_dir,subj_info.subj_id,'model_sim');
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

all_dist=[];
all_prox=[];
all_cum=[];
for i=0:49
    fid=fopen(fullfile('../../data/model/18feb16_beta_ongoing_dist_only_F/data',sprintf('dpl_%d.txt',i)),'r');
    C=textscan(fid,'%f');
    trial_data=reshape(C{1},4,length(C{1})/4);
    t_idx=find(trial_data(1,:)>=50);
    all_dist(i+1,:)=trial_data(2,t_idx);
    fclose(fid);
    
    fid=fopen(fullfile('../../data/model/18feb16_beta_ongoing_prox_only_E/data',sprintf('dpl_%d.txt',i)),'r');
    C=textscan(fid,'%f');
    trial_data=reshape(C{1},4,length(C{1})/4);
    t_idx=find(trial_data(1,:)>=50);
    all_prox(i+1,:)=trial_data(2,t_idx);
    fclose(fid);

    fid=fopen(fullfile('../../data/model/18feb5_beta_ongoing_D',sprintf('dpl_%d.txt',i)),'r');
    C=textscan(fid,'%f');
    trial_data=reshape(C{1},4,length(C{1})/4);
    t_idx=find(trial_data(1,:)>=50);
    all_cum(i+1,:)=trial_data(2,t_idx);
    fclose(fid);
end

ds_dist=ft_preproc_resample(all_dist,4.0000e+04,250,'resample');
ds_prox=ft_preproc_resample(all_prox,4.0000e+04,250,'resample');
ds_cum=ft_preproc_resample(all_cum,4.0000e+04,250,'resample');

% Copy data file
reg_file=fullfile(output_data_dir, 'grey_rcresp_TafdfC.mat');

if sim_idx==1
    clear jobs
    matlabbatch=[];
    matlabbatch{1}.spm.meeg.other.copy.D = {data_file};
    matlabbatch{1}.spm.meeg.other.copy.outfile = fullfile(output_data_dir, 'grey_rcresp_TafdfC.mat');
    spm_jobman('run', matlabbatch);

    % Remove extra trials
    D=spm_eeg_load(fullfile(output_data_dir, 'grey_rcresp_TafdfC.mat'));
    D1=D.badtrials(1:size(D,3)-size(ds_prox,1),true);
    S=[];
    S.D=D1;
    S.prefix='r';
    D2=spm_eeg_remove_bad_trials(S);
    
    clear jobs
    matlabbatch=[];
    matlabbatch{1}.spm.meeg.other.copy.D = {fullfile(output_data_dir, 'rgrey_rcresp_TafdfC.mat')};
    matlabbatch{1}.spm.meeg.other.copy.outfile = reg_file;
    spm_jobman('run', matlabbatch);

    delete(fullfile(output_data_dir, 'rgrey_rcresp_TafdfC.mat'));
    delete(fullfile(output_data_dir, 'rgrey_rcresp_TafdfC.dat'));
    
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
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.modality = {'All'};
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.crossval = [pctest 1];      
    spm_jobman('run', matlabbatch);
end

% Load inverted data file
D=spm_eeg_load(reg_file);   

% Create simulated signal
woi=[D.time(1)-(D.time(2)-D.time(1)) D.time(end)];

% Get dipole orientations
pial_unit_norm=pial.normals(pial_vertex,:);
white_unit_norm=-1.*white.normals(pial_vertex,:);

% Simulate source 
deep_signal = zeros(1,length(D.time),size(ds_prox,1));
deep_signal(1,round(length(D.time)/2-size(ds_prox,2)/2)+1:round(length(D.time)/2+size(ds_prox,2)/2),:)=ds_prox';
superficial_signal = zeros(1,length(D.time),size(ds_dist,1));
superficial_signal(1,round(length(D.time)/2-size(ds_dist,2)/2)+1:round(length(D.time)/2+size(ds_dist,2)/2),:)=ds_dist';
cum_signal = zeros(1,length(D.time),size(ds_cum,1));
cum_signal(1,round(length(D.time)/2-size(ds_cum,2)/2)+1:round(length(D.time)/2+size(ds_cum,2)/2),:)=ds_cum';

% Normalize
deep_mag=max(abs(squeeze(mean(deep_signal,3))));
deep_signal=deep_signal./deep_mag;
superficial_mag=max(abs(squeeze(mean(superficial_signal,3))));
superficial_signal=superficial_signal./superficial_mag;

% Flip superficial signal
superficial_signal=-1.*superficial_signal;

SNRdB=-20;
%SNRdB=-10;
noiseFt=[];

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
    %nAmdipmom=[31.7187 67.6499];
    % Combined signal
    simsignal=[deep_signal; superficial_signal];
    % Width of patch
    dipfwhm=[5 5];
elseif strcmp(type,'reversed')
    % White wide - pial narrow
    % Coordinate of each signal
    simpos=[D.inv{1}.mesh.tess_mni.vert(grey_white_vertex,:)
        D.inv{1}.mesh.tess_mni.vert(grey_pial_vertex,:)];
    % Orientation of each dipole
    ormni=[white_unit_norm; pial_unit_norm];
    % Dipole momemnts
    nAmdipmom=[8 6];
    %nAmdipmom=[67.6499 31.7187];
    % Combined signal
    simsignal=[superficial_signal; deep_signal];
    % Width of patch
    dipfwhm=[5 5];
elseif strcmp(type,'pial_only')
    % Pial only
    simpos=[D.inv{1}.mesh.tess_mni.vert(grey_pial_vertex,:)];
    ormni=[pial_unit_norm];
    nAmdipmom=[8];
    %nAmdipmom=[67.6499];
    simsignal=[superficial_signal];
    dipfwhm=[5];
elseif strcmp(type,'white_only')
    % White only
    simpos=[D.inv{1}.mesh.tess_mni.vert(grey_white_vertex,:)];
    ormni=[white_unit_norm];
    nAmdipmom=[6];
    %nAmdipmom=[31.7187];
    simsignal=[deep_signal];
    dipfwhm=[5];
elseif strcmp(type,'reversed_ori')
    % White wide - pial narrow
    % Coordinate of each signal
    simpos=[D.inv{1}.mesh.tess_mni.vert(grey_white_vertex,:)
        D.inv{1}.mesh.tess_mni.vert(grey_pial_vertex,:)];
    % Orientation of each dipole
    ormni=[pial_unit_norm; white_unit_norm];
    % Dipole momemnts
    nAmdipmom=[6 8];
    %nAmdipmom=[31.7187 67.6499];
    % Combined signal
    simsignal=[deep_signal; superficial_signal];
    % Width of patch
    dipfwhm=[5 5];
elseif strcmp(type,'pial_ori')
    % White wide - pial narrow
    % Coordinate of each signal
    simpos=[D.inv{1}.mesh.tess_mni.vert(grey_white_vertex,:)
        D.inv{1}.mesh.tess_mni.vert(grey_pial_vertex,:)];
    % Orientation of each dipole
    ormni=[pial_unit_norm; pial_unit_norm];
    % Dipole momemnts
    nAmdipmom=[6 8];
    %nAmdipmom=[31.7187 67.6499];
    % Combined signal
    simsignal=[deep_signal; superficial_signal];
    % Width of patch
    dipfwhm=[5 5];
elseif strcmp(type,'white_ori')
    % White wide - pial narrow
    % Coordinate of each signal
    simpos=[D.inv{1}.mesh.tess_mni.vert(grey_white_vertex,:)
        D.inv{1}.mesh.tess_mni.vert(grey_pial_vertex,:)];
    % Orientation of each dipole
    ormni=[white_unit_norm; white_unit_norm];
    % Dipole momemnts
    nAmdipmom=[6 8];
    %nAmdipmom=[31.7187 67.6499];
    % Combined signal
    simsignal=[deep_signal; superficial_signal];
    % Width of patch
    dipfwhm=[5 5];
elseif strcmp(type,'none')
    %% White wide - pial narrow
    % Coordinate of each signal
    simpos=[D.inv{1}.mesh.tess_mni.vert(grey_white_vertex,:)
        D.inv{1}.mesh.tess_mni.vert(grey_pial_vertex,:)];
    % Orientation of each dipole
    ormni=[white_unit_norm; pial_unit_norm];
    % Dipole momemnts
    nAmdipmom=[0 0];
    %nAmdipmom=[31.7187 67.6499];
    % Combined signal
    simsignal=[deep_signal; superficial_signal];
    % Width of patch
    dipfwhm=[5 5];
    SNRdB=[];
    noiseFt=100;
end

times=D.time.*1000;
zero_time=times(round(length(times)/2))+.5*(times(round(length(times)/2+1))-times(round(length(times)/2)));
times=times-zero_time;
 
% fig=figure();
% subplot(3,1,1);
% hold all;
% pial_signal=nAmdipmom(2).*squeeze(simsignal(2,:,:))';
% % plot(times, nAmdipmom(2).*squeeze(simsignal(2,:,:))','Color',[240 140 100]./255.0);
% % plot(times, nAmdipmom(2).*squeeze(mean(simsignal(2,:,:),3)),'Color','k');
% % ylabel('Dipole moment (nA)');
% % xlim(times([1 end]));
% % %ylim([-2 13]);
% % 
% % subplot(3,1,2);
% % hold all;
% % cumulative_signal=squeeze(cum_signal(1,:,:))';
% % plot(times, cumulative_signal,'Color',[.75 .75 .75]);
% % plot(times, mean(cumulative_signal),'Color','k');
% % ylabel('Dipole moment (nA)');
% % xlim(times([1 end]));
% % %ylim([-12 7]);
% % 
% % subplot(3,1,3);
% % hold all;
% % plot(times, nAmdipmom(1).*squeeze(simsignal(1,:,:))','Color',[106 175 215]./255.0);
% % plot(times, nAmdipmom(1).*squeeze(mean(simsignal(1,:,:),3)),'Color','k');
% % ylabel('Dipole moment (nA)');
% % xlim(times([1 end]));
% % %ylim([-2 13]);


% Set sim signal to have unit variance
% simsignal=simsignal./repmat(std(simsignal,[],2),1,size(simsignal,2),1);
[D,meshsourceind]=spm_eeg_simulate({D}, sprintf('sim_%s_%d_',type,sim_idx), simpos,...
   simsignal, ormni, woi, noiseFt, SNRdB, [], [], dipfwhm,...
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

plot_model_sim_burst_sensor(subj_info, type, sim_idx, 'base_dir', params.base_dir)

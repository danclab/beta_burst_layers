function preprocess(data_dir)

clear jobs
matlabbatch={};
batch_idx=1;

% Copy data file (before averaging over bursts)
data_file=fullfile(data_dir, 'rcresp_TafdfC.mat');
matlabbatch{batch_idx}.spm.meeg.other.copy.D = {data_file};
matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(data_dir, 'rcresp_TafdfC_shuffled.mat');
batch_idx=batch_idx+1;
spm_jobman('run', matlabbatch); 

% Shuffle times in each trial
D=spm_eeg_load(fullfile(data_dir, 'rcresp_TafdfC_shuffled.mat'));
for trl_idx=1:size(D,3)
    permuted_times=randperm(size(D,2));
    D(:,:,trl_idx)=D(:,permuted_times,trl_idx);
end
save(fullfile(data_dir, 'rcresp_TafdfC_shuffled.mat'),'D');

% Average over bursts
clear jobs
matlabbatch={};
batch_idx=1;

matlabbatch{batch_idx}.spm.meeg.averaging.average.D = {fullfile(data_dir, 'rcresp_TafdfC_shuffled.mat')};
matlabbatch{batch_idx}.spm.meeg.averaging.average.userobust.standard = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.plv = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.prefix = fullfile(data_dir, 'm');
batch_idx=batch_idx+1;

spm_jobman('run', matlabbatch); 
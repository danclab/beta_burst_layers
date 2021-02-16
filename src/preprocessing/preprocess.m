function preprocess(raw_data_dir, output_data_dir)

clear jobs
matlabbatch={};
batch_idx=1;

% List of files for each session
session_files={};
d=dir(raw_data_dir);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..','plots'})) = [];
for idx=1:length(nameFolds)
    session_num=nameFolds{idx};
    session_dir=fullfile(raw_data_dir,session_num);
    session_files{end+1,1}=fullfile(session_dir, sprintf('rcresp_Tafdf%sC.mat', session_num));
end

% If more than one session
if length(nameFolds)>1
    % Merge the files
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.D = session_files;
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.file = '.*';
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelorg = '.*';
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelnew = '#labelorg#';
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.prefix = 'c';
    batch_idx=batch_idx+1;    

    % Copy to subject's output data dir
    [path filename]=fileparts(session_files{end,1});
    matlabbatch{batch_idx}.spm.meeg.other.copy.D = {sprintf('c%s.mat', filename)};
    matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(output_data_dir, 'rcresp_TafdfC.mat');
    batch_idx=batch_idx+1;

    % Delete original merged file
    matlabbatch{batch_idx}.spm.meeg.other.delete.D = {sprintf('c%s.mat', filename)};
    batch_idx=batch_idx+1;
    
% If just one session, copy to subject's output data dir
else
    matlabbatch{batch_idx}.spm.meeg.other.copy.D = {session_files{end,1}};
    matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(output_data_dir, 'rcresp_TafdfC.mat');
    batch_idx=batch_idx+1;
end

% Average over bursts
matlabbatch{batch_idx}.spm.meeg.averaging.average.D = {fullfile(output_data_dir, 'rcresp_TafdfC.mat')};
matlabbatch{batch_idx}.spm.meeg.averaging.average.userobust.standard = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.plv = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.prefix = fullfile(output_data_dir, 'm');
batch_idx=batch_idx+1;

spm_jobman('run', matlabbatch); 
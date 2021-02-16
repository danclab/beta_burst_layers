function plot_subjects_tc(subjects, erd_dir, ers_dir, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

num_spatial_runs=20;

erd_base_dir_parts=strsplit(erd_dir,filesep);
ers_base_dir_parts=strsplit(ers_dir,filesep);

erd_spatial_control_subj_mean_fdiffs=[];
ers_spatial_control_subj_mean_fdiffs=[];

for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);

    g=gifti(fullfile('../../data/surf',sprintf('%s-synth',subj_info.subj_id),'surf','white.ds-pial.ds.link_vector.gii'));
    
    % Get sliding time window results
    data_dir=fullfile('../../output/data',erd_base_dir_parts{end},subj_info.subj_id);
    load(fullfile(data_dir, 'invert_burst_tc_results.mat'));

    sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
    
    for r_idx=1:num_spatial_runs
        % Get sliding time window results
        load(fullfile(data_dir, sprintf('invert_burst_tc_results-random_locations_%d.mat',r_idx)));
        pial_clusters=invert_burst_tc_results.clusters;
        f_diffs=[];
        for c_idx=1:length(pial_clusters)
            f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        end
        erd_spatial_control_subj_mean_fdiffs(r_idx,s_idx,:)=mean(f_diffs,1);
    end
    
    % Get sliding time window results
    data_dir=fullfile('../../output/data',ers_base_dir_parts{end},subj_info.subj_id);
    load(fullfile(data_dir, 'invert_burst_tc_results.mat'));

    sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
    
    % Get source tc
    D=spm_eeg_load(fullfile(data_dir, 'grey_mrcresp_TafdfC.mat'));
    
    for r_idx=1:num_spatial_runs
        % Get sliding time window results
        load(fullfile(data_dir, sprintf('invert_burst_tc_results-random_locations_%d.mat',r_idx)));
        pial_clusters=invert_burst_tc_results.clusters;
        f_diffs=[];
        for c_idx=1:length(pial_clusters)
            f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        end
        ers_spatial_control_subj_mean_fdiffs(r_idx,s_idx,:)=mean(f_diffs,1);
    end

    times=D.time.*1000;
    zero_time=times(round(length(times)/2))+.5*(times(round(length(times)/2+1))-times(round(length(times)/2)));
    times=times-zero_time;
end

for r_idx=1:num_spatial_runs
    figure();
	subplot(2,1,1);
    plot(sliding_tc_times,squeeze(sum(erd_spatial_control_subj_mean_fdiffs(r_idx,:,:),2)));
    subplot(2,1,2);
    plot(sliding_tc_times,squeeze(sum(ers_spatial_control_subj_mean_fdiffs(r_idx,:,:),2)));
end

figure();
subplot(2,1,1);
plot(sliding_tc_times,squeeze(sum(erd_spatial_control_subj_mean_fdiffs,2))');
hold all
plot(sliding_tc_times,squeeze(mean(sum(erd_spatial_control_subj_mean_fdiffs,2)))','k','LineWidth',2);
subplot(2,1,2);
plot(sliding_tc_times,squeeze(sum(ers_spatial_control_subj_mean_fdiffs,2))');
hold all
plot(sliding_tc_times,squeeze(mean(sum(ers_spatial_control_subj_mean_fdiffs,2)))','k','LineWidth',2);
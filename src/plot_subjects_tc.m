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
num_temporal_runs=20;

erd_base_dir_parts=strsplit(erd_dir,filesep);
ers_base_dir_parts=strsplit(ers_dir,filesep);

erd_subj_mean_fdiffs=[];
ers_subj_mean_fdiffs=[];
erd_subj_mean_source_tcs=[];
ers_subj_mean_source_tcs=[];
erd_subj_scaled_mean_source_tcs=[];
ers_subj_scaled_mean_source_tcs=[];
erd_spatial_control_subj_mean_fdiffs=[];
ers_spatial_control_subj_mean_fdiffs=[];
erd_temporal_control_subj_mean_fdiffs=[];
ers_temporal_control_subj_mean_fdiffs=[];

for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);

    g=gifti(fullfile('../data/surf',sprintf('%s-synth',subj_info.subj_id),'surf','white.ds-pial.ds.link_vector.gii'));
    
    % Get sliding time window results
    data_dir=fullfile('../output/data',erd_base_dir_parts{end},subj_info.subj_id);
    load(fullfile(data_dir, 'invert_burst_tc_results.mat'));

    sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
    
    % Get source tc
    D=spm_eeg_load(fullfile(data_dir, 'grey_mrcresp_TafdfC.mat'));
    MU=D.inv{1}.inverse.M*D.inv{1}.inverse.U{1};
    megchans=D.indchantype('meg','good');
    legend_labels={};
    pial_clusters=invert_burst_tc_results.clusters;
    
    source_tcs=[];
    f_diffs=[];
    for c_idx=1:length(pial_clusters)
        verts=size(MU,1)/2+pial_clusters(c_idx).inv_verts;    
        source_tc=mean(MU(verts,:)*squeeze(D(megchans,:,:)),1);
        if mean(source_tc(26:31))>0
            source_tc=source_tc.*-1;
        end
        f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        source_tcs(end+1,:)=source_tc;        
    end
    erd_subj_mean_source_tcs(s_idx,:)=mean(source_tcs,1);
    erd_subj_scaled_mean_source_tcs(s_idx,:)=erd_subj_mean_source_tcs(s_idx,:)./max(abs(erd_subj_mean_source_tcs(s_idx,:)));
    erd_subj_mean_fdiffs(s_idx,:)=mean(f_diffs,1);

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
    for r_idx=1:num_temporal_runs
        load(fullfile(data_dir, sprintf('invert_burst_tc_results-shuffled_times_%d.mat',r_idx)));
        pial_clusters=invert_burst_tc_results.clusters;
        f_diffs=[];
        for c_idx=1:length(pial_clusters)
            f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        end
        erd_temporal_control_subj_mean_fdiffs(r_idx,s_idx,:)=mean(f_diffs,1);
    end
    
    % Get sliding time window results
    data_dir=fullfile('../output/data',ers_base_dir_parts{end},subj_info.subj_id);
    load(fullfile(data_dir, 'invert_burst_tc_results.mat'));

    sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
    
    % Get source tc
    D=spm_eeg_load(fullfile(data_dir, 'grey_mrcresp_TafdfC.mat'));
    MU=D.inv{1}.inverse.M*D.inv{1}.inverse.U{1};
    megchans=D.indchantype('meg','good');
    legend_labels={};
    pial_clusters=invert_burst_tc_results.clusters;
    
    source_tcs=[];
    f_diffs=[];
    for c_idx=1:length(pial_clusters)
        verts=size(MU,1)/2+pial_clusters(c_idx).inv_verts;    
        source_tc=mean(MU(verts,:)*squeeze(D(megchans,:,:)),1);
        if mean(source_tc(26:31))>0
            source_tc=source_tc.*-1;
        end
        f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        source_tcs(end+1,:)=source_tc;   
    end
    ers_subj_mean_source_tcs(s_idx,:)=mean(source_tcs,1);
    ers_subj_scaled_mean_source_tcs(s_idx,:)=ers_subj_mean_source_tcs(s_idx,:)./max(abs(ers_subj_mean_source_tcs(s_idx,:)));
    ers_subj_mean_fdiffs(s_idx,:)=mean(f_diffs,1);

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
    
    % Get sliding time window results
    for r_idx=1:num_temporal_runs
        load(fullfile(data_dir, sprintf('invert_burst_tc_results-shuffled_times_%d.mat',r_idx)));
        pial_clusters=invert_burst_tc_results.clusters;
        f_diffs=[];
        for c_idx=1:length(pial_clusters)
            f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        end
        ers_temporal_control_subj_mean_fdiffs(r_idx,s_idx,:)=mean(f_diffs,1);
    end

    times=D.time.*1000;
    zero_time=times(round(length(times)/2))+.5*(times(round(length(times)/2+1))-times(round(length(times)/2)));
    times=times-zero_time;
end

[wout,lags]=woody(erd_subj_scaled_mean_source_tcs',[],[],'woody','biased');
   
erd_n_aligned=size(erd_subj_scaled_mean_source_tcs,2)-abs(min(lags))-max(lags);
erd_n_inv_aligned=size(erd_subj_mean_fdiffs,2)-abs(min(lags))-max(lags);

erd_aligned_tcs=[];
erd_aligned_mean_f_diffs=[];
erd_aligned_spatial_control_mean_f_diffs=[];
erd_aligned_temporal_control_mean_f_diffs=[];
for subj_idx=1:size(erd_subj_scaled_mean_source_tcs,1)
    start_idx=1+(lags(subj_idx)-min(lags));
    erd_aligned_tcs(subj_idx,:)=erd_subj_mean_source_tcs(subj_idx,start_idx:start_idx+erd_n_aligned-1);
    erd_aligned_mean_f_diffs(subj_idx,:)=erd_subj_mean_fdiffs(subj_idx,start_idx:start_idx+erd_n_inv_aligned-1);
    erd_aligned_spatial_control_mean_f_diffs(:,subj_idx,:)=erd_spatial_control_subj_mean_fdiffs(:,subj_idx,start_idx:start_idx+erd_n_inv_aligned-1);
    erd_aligned_temporal_control_mean_f_diffs(:,subj_idx,:)=erd_temporal_control_subj_mean_fdiffs(:,subj_idx,start_idx:start_idx+erd_n_inv_aligned-1);
    
    erd_aligned_times=times(start_idx:start_idx+erd_n_aligned-1);
    erd_aligned_sliding_tc_times=sliding_tc_times(start_idx:start_idx+erd_n_inv_aligned-1);
end
   

[wout,lags]=woody(ers_subj_scaled_mean_source_tcs',[],[],'woody','biased');
   
ers_n_aligned=size(ers_subj_scaled_mean_source_tcs,2)-abs(min(lags))-max(lags);
ers_n_inv_aligned=size(ers_subj_mean_fdiffs,2)-abs(min(lags))-max(lags);

ers_aligned_tcs=[];
ers_aligned_mean_f_diffs=[];
ers_aligned_spatial_control_mean_f_diffs=[];
ers_aligned_temporal_control_mean_f_diffs=[];
for subj_idx=1:size(ers_subj_scaled_mean_source_tcs,1)
    start_idx=1+(lags(subj_idx)-min(lags));
    ers_aligned_tcs(subj_idx,:)=ers_subj_mean_source_tcs(subj_idx,start_idx:start_idx+ers_n_aligned-1);
    ers_aligned_mean_f_diffs(subj_idx,:)=ers_subj_mean_fdiffs(subj_idx,start_idx:start_idx+ers_n_inv_aligned-1);
    ers_aligned_spatial_control_mean_f_diffs(:,subj_idx,:)=ers_spatial_control_subj_mean_fdiffs(:,subj_idx,start_idx:start_idx+ers_n_inv_aligned-1);
    ers_aligned_temporal_control_mean_f_diffs(:,subj_idx,:)=ers_temporal_control_subj_mean_fdiffs(:,subj_idx,start_idx:start_idx+ers_n_inv_aligned-1);
    
    ers_aligned_times=times(start_idx:start_idx+ers_n_aligned-1);
    ers_aligned_sliding_tc_times=sliding_tc_times(start_idx:start_idx+ers_n_inv_aligned-1);
end

figure();
subplot(2,2,1);
hold on
%plot(times,erd_subj_mean_source_tcs.'*1000);
shadedErrorBar(times,...
    mean(erd_subj_mean_source_tcs.*1000),...
    std(erd_subj_mean_source_tcs.*1000)./sqrt(length(subjects)),...
    'LineProps',{'Color','k','LineWidth',2});
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('pAm/mm^2');

subplot(2,2,2);
hold on
%plot(times,ers_subj_mean_source_tcs'.*1000);
shadedErrorBar(times,...
    mean(ers_subj_mean_source_tcs.*1000),...
    std(ers_subj_mean_source_tcs.*1000)./sqrt(length(subjects)),...
    'LineProps',{'Color','k','LineWidth',2});
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('pAm/mm^2');

subplot(2,2,3);
hold all
plot(sliding_tc_times,sum(erd_subj_mean_fdiffs,1),'k','LineWidth',2);
shadedErrorBar(sliding_tc_times,...
    squeeze(mean(sum(erd_spatial_control_subj_mean_fdiffs,2))),...
    squeeze(std(sum(erd_spatial_control_subj_mean_fdiffs,2)))./sqrt(num_spatial_runs),...
    'LineProps',{'Color','r'});
shadedErrorBar(sliding_tc_times,...
    squeeze(mean(sum(erd_temporal_control_subj_mean_fdiffs,2))),...
    squeeze(std(sum(erd_temporal_control_subj_mean_fdiffs,2)))./sqrt(num_temporal_runs),...
    'LineProps',{'Color','b'});
plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylim([-75 75]);
ylabel('\Delta F');
xlabel('Time (ms)');

subplot(2,2,4);
hold all
plot(sliding_tc_times,sum(ers_subj_mean_fdiffs,1),'k','LineWidth',2);
shadedErrorBar(sliding_tc_times,...
    squeeze(mean(sum(ers_spatial_control_subj_mean_fdiffs,2))),...
    squeeze(std(sum(ers_spatial_control_subj_mean_fdiffs,2)))./sqrt(num_spatial_runs),...
    'LineProps',{'Color','r'});
shadedErrorBar(sliding_tc_times,...
    squeeze(mean(sum(ers_temporal_control_subj_mean_fdiffs,2))),...
    squeeze(std(sum(ers_temporal_control_subj_mean_fdiffs,2)))./sqrt(num_temporal_runs),...
    'LineProps',{'Color','b'});
plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylim([-75 75]);
legend('data','spatial control','temporal control');
ylabel('\Delta F');
xlabel('Time (ms)');

figure();
subplot(2,2,1);
hold all
%plot(erd_aligned_times,erd_aligned_tcs'.*1000);
shadedErrorBar(erd_aligned_times,...
    mean(erd_aligned_tcs.*1000),...
    std(erd_aligned_tcs.*1000)./sqrt(length(subjects)),...
    'LineProps',{'Color','k','LineWidth',2});
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('pAm/mm^2');

subplot(2,2,2);
hold all
%plot(ers_aligned_times,ers_aligned_tcs'.*1000);
shadedErrorBar(ers_aligned_times,...
    mean(ers_aligned_tcs.*1000),...
    std(ers_aligned_tcs.*1000)./sqrt(length(subjects)),...
    'LineProps',{'Color','k','LineWidth',2});
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('pAm/mm^2');

subplot(2,2,3);
hold all
plot(erd_aligned_sliding_tc_times,sum(erd_aligned_mean_f_diffs,1),'k','LineWidth',2);
shadedErrorBar(erd_aligned_sliding_tc_times,...
    squeeze(mean(sum(erd_aligned_spatial_control_mean_f_diffs,2))),...
    squeeze(std(sum(erd_aligned_spatial_control_mean_f_diffs,2)))./sqrt(num_spatial_runs),...
    'LineProps',{'Color','r'});
shadedErrorBar(erd_aligned_sliding_tc_times,...
    squeeze(mean(sum(erd_aligned_temporal_control_mean_f_diffs,2))),...
    squeeze(std(sum(erd_aligned_temporal_control_mean_f_diffs,2)))./sqrt(num_spatial_runs),...
    'LineProps',{'Color','b'});
plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylim([-75 75]);
ylabel('\Delta F');
xlabel('Time (ms)');

subplot(2,2,4);
hold all
plot(ers_aligned_sliding_tc_times,sum(ers_aligned_mean_f_diffs,1),'k','LineWidth',2);
shadedErrorBar(ers_aligned_sliding_tc_times,...
    squeeze(mean(sum(ers_aligned_spatial_control_mean_f_diffs,2))),...
    squeeze(std(sum(ers_aligned_spatial_control_mean_f_diffs,2)))./sqrt(num_spatial_runs),...
    'LineProps',{'Color','r'});
shadedErrorBar(ers_aligned_sliding_tc_times,...
    squeeze(mean(sum(ers_aligned_temporal_control_mean_f_diffs,2))),...
    squeeze(std(sum(ers_aligned_temporal_control_mean_f_diffs,2)))./sqrt(num_spatial_runs),...
    'LineProps',{'Color','b'});
plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylim([-75 75]);
legend('data','spatial control','temporal control');
ylabel('\Delta F');
xlabel('Time (ms)');

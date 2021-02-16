function plot_subjects_tc(subjects, erd_dir, ers_dir, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

window_lengths=[2 5 10 15];
labels={};
for i=1:length(window_lengths)
    labels{i}=sprintf('%dms', window_lengths(i)*4);
end

erd_base_dir_parts=strsplit(erd_dir,filesep);
ers_base_dir_parts=strsplit(ers_dir,filesep);

erd_subj_mean_fdiffs={};
ers_subj_mean_fdiffs={};
sliding_tc_times={};
erd_subj_scaled_mean_source_tcs=[];
ers_subj_scaled_mean_source_tcs=[];

for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);

    % Get sliding time window results
    data_dir=fullfile('../../output/data',erd_base_dir_parts{end},subj_info.subj_id);
    for i=1:length(window_lengths)
        if window_lengths(i)==10
            fname='invert_burst_tc_results.mat';
        else
            fname=sprintf('invert_burst_tc_results_%d_window_size.mat', window_lengths(i));
        end
        load(fullfile(data_dir, fname));
        pial_clusters=invert_burst_tc_results.clusters;
        f_diffs=[];
        for c_idx=1:length(pial_clusters)
            f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        end
        erd_subj_mean_fdiffs{i}(s_idx,:)=mean(f_diffs,1);
        sliding_tc_times{i}=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
    end
        
    % Get source tc
    D=spm_eeg_load(fullfile(data_dir, 'grey_mrcresp_TafdfC.mat'));
    MU=D.inv{1}.inverse.M*D.inv{1}.inverse.U{1};
    megchans=D.indchantype('meg','good');
    legend_labels={};
    pial_clusters=invert_burst_tc_results.clusters;
    
    source_tcs=[];
    for c_idx=1:length(pial_clusters)
        verts=size(MU,1)/2+pial_clusters(c_idx).inv_verts;    
        source_tc=mean(MU(verts,:)*squeeze(D(megchans,:,:)),1);
        if mean(source_tc(26:31))>0
            source_tc=source_tc.*-1;
        end
        source_tcs(end+1,:)=source_tc;        
    end
    mean_source_tcs=mean(source_tcs,1);
    erd_subj_scaled_mean_source_tcs(s_idx,:)=mean_source_tcs./max(abs(mean_source_tcs));
    
    % Get sliding time window results
    data_dir=fullfile('../../output/data',ers_base_dir_parts{end},subj_info.subj_id);
    for i=1:length(window_lengths)
        if window_lengths(i)==10
            fname='invert_burst_tc_results.mat';
        else
            fname=sprintf('invert_burst_tc_results_%d_window_size.mat', window_lengths(i));
        end
        load(fullfile(data_dir, fname));
        pial_clusters=invert_burst_tc_results.clusters;
        f_diffs=[];
        for c_idx=1:length(pial_clusters)
            f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        end
        ers_subj_mean_fdiffs{i}(s_idx,:)=mean(f_diffs,1);
        sliding_tc_times{i}=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
    end
    
    % Get source tc
    D=spm_eeg_load(fullfile(data_dir, 'grey_mrcresp_TafdfC.mat'));
    MU=D.inv{1}.inverse.M*D.inv{1}.inverse.U{1};
    megchans=D.indchantype('meg','good');
    legend_labels={};
    pial_clusters=invert_burst_tc_results.clusters;
    
    source_tcs=[];
    for c_idx=1:length(pial_clusters)
        verts=size(MU,1)/2+pial_clusters(c_idx).inv_verts;    
        source_tc=mean(MU(verts,:)*squeeze(D(megchans,:,:)),1);
        if mean(source_tc(26:31))>0
            source_tc=source_tc.*-1;
        end
        source_tcs(end+1,:)=source_tc;        
    end
    mean_source_tcs=mean(source_tcs,1);
    ers_subj_scaled_mean_source_tcs(s_idx,:)=mean_source_tcs./max(abs(mean_source_tcs));
    
    times=D.time.*1000;
    zero_time=times(round(length(times)/2))+.5*(times(round(length(times)/2+1))-times(round(length(times)/2)));
    times=times-zero_time;
end

addpath('..');
[wout,lags]=woody(erd_subj_scaled_mean_source_tcs',[],[],'woody','biased');
   
erd_n_aligned=size(erd_subj_scaled_mean_source_tcs,2)-abs(min(lags))-max(lags);

erd_aligned_tcs=[];
erd_aligned_mean_f_diffs={};
erd_aligned_sliding_tc_times={};
for subj_idx=1:size(erd_subj_scaled_mean_source_tcs,1)
    start_idx=1+(lags(subj_idx)-min(lags));
    erd_aligned_tcs(subj_idx,:)=erd_subj_scaled_mean_source_tcs(subj_idx,start_idx:start_idx+erd_n_aligned-1);
    for i=1:length(window_lengths)
        erd_n_inv_aligned=size(erd_subj_mean_fdiffs{i},2)-abs(min(lags))-max(lags);
        erd_aligned_mean_f_diffs{i}(subj_idx,:)=erd_subj_mean_fdiffs{i}(subj_idx,start_idx:start_idx+erd_n_inv_aligned-1);
        erd_aligned_sliding_tc_times{i}=sliding_tc_times{i}(start_idx:start_idx+erd_n_inv_aligned-1);
    end
    
    erd_aligned_times=times(start_idx:start_idx+erd_n_aligned-1);
    
end
   

[wout,lags]=woody(ers_subj_scaled_mean_source_tcs',[],[],'woody','biased');
   
ers_n_aligned=size(ers_subj_scaled_mean_source_tcs,2)-abs(min(lags))-max(lags);

ers_aligned_tcs=[];
ers_aligned_mean_f_diffs={};
ers_aligned_sliding_tc_times={};
for subj_idx=1:size(ers_subj_scaled_mean_source_tcs,1)
    start_idx=1+(lags(subj_idx)-min(lags));
    ers_aligned_tcs(subj_idx,:)=ers_subj_scaled_mean_source_tcs(subj_idx,start_idx:start_idx+ers_n_aligned-1);
    for i=1:length(window_lengths)
        ers_n_inv_aligned=size(ers_subj_mean_fdiffs{i},2)-abs(min(lags))-max(lags);
        ers_aligned_mean_f_diffs{i}(subj_idx,:)=ers_subj_mean_fdiffs{i}(subj_idx,start_idx:start_idx+ers_n_inv_aligned-1);
        ers_aligned_sliding_tc_times{i}=sliding_tc_times{i}(start_idx:start_idx+ers_n_inv_aligned-1);
    end    
    ers_aligned_times=times(start_idx:start_idx+ers_n_aligned-1);
    
end
rmpath('..');

colors=cbrewer('seq','Blues',length(window_lengths)+1);

figure();
subplot(2,2,1);
hold all
for i=1:length(window_lengths)
    color=colors(i+1,:);
    if window_lengths(i)==10
        color='k';
    end
    plot(sliding_tc_times{i},squeeze(sum(erd_subj_mean_fdiffs{i}(:,:),1)),'Color',color,'LineWidth',2);
end
plot([sliding_tc_times{1}(1) sliding_tc_times{1}(end)],[3 3],'k--');
plot([sliding_tc_times{1}(1) sliding_tc_times{1}(end)],[-3 -3],'k--');
legend(labels);
xlim([sliding_tc_times{1}(1) sliding_tc_times{1}(end)]);
%ylim([-75 75]);
ylabel('Sum mean \Delta F');
xlabel('Time (ms)');

subplot(2,2,2);
hold all
for i=1:length(window_lengths)
    color=colors(i+1,:);
    if window_lengths(i)==10
        color='k';
    end
    plot(sliding_tc_times{i},squeeze(sum(ers_subj_mean_fdiffs{i}(:,:),1)),'Color',color,'LineWidth',2);
end
plot([sliding_tc_times{1}(1) sliding_tc_times{1}(end)],[3 3],'k--');
plot([sliding_tc_times{1}(1) sliding_tc_times{1}(end)],[-3 -3],'k--');
xlim([sliding_tc_times{1}(1) sliding_tc_times{1}(end)]);
%ylim([-75 75]);
ylabel('Sum mean \Delta F');
xlabel('Time (ms)');

subplot(2,2,3);
hold all
for i=1:length(window_lengths)
    color=colors(i+1,:);
    if window_lengths(i)==10
        color='k';
    end
    plot(erd_aligned_sliding_tc_times{i},squeeze(sum(erd_aligned_mean_f_diffs{i}(:,:),1)),'Color',color,'LineWidth',2);
end
plot([sliding_tc_times{1}(1) sliding_tc_times{1}(end)],[3 3],'k--');
plot([sliding_tc_times{1}(1) sliding_tc_times{1}(end)],[-3 -3],'k--');
xlim([sliding_tc_times{1}(1) sliding_tc_times{1}(end)]);
%ylim([-75 75]);
ylabel('Sum aligned mean \Delta F');
xlabel('Time (ms)');

subplot(2,2,4);
hold all
for i=1:length(window_lengths)
    color=colors(i+1,:);
    if window_lengths(i)==10
        color='k';
    end
    plot(ers_aligned_sliding_tc_times{i},squeeze(sum(ers_aligned_mean_f_diffs{i}(:,:),1)),'Color',color,'LineWidth',2);
end
plot([sliding_tc_times{1}(1) sliding_tc_times{1}(end)],[3 3],'k--');
plot([sliding_tc_times{1}(1) sliding_tc_times{1}(end)],[-3 -3],'k--');
xlim([sliding_tc_times{1}(1) sliding_tc_times{1}(end)]);
%ylim([-75 75]);
ylabel('Sum aligned mean \Delta F');
xlabel('Time (ms)');

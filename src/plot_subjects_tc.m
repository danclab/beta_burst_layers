function plot_subjects_tc(subjects, varargin)

defaults = struct('base_dir', '../../data/JB_BUTTON_LOCKED_d3_ers');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

base_dir_parts=strsplit(params.base_dir,filesep);
subj_mean_bc_fdiffs=[];
subj_scaled_mean_source_tcs=[];

for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);

    % Get sliding time window results
    data_dir=fullfile('../output/data',base_dir_parts{end},subj_info.subj_id);
    load(fullfile(data_dir, 'invert_burst_tc_results.mat'));

    sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
    
    % Get source tc
    D=spm_eeg_load(fullfile(data_dir, 'grey_mrcresp_TafdfC.mat'));
    MU=D.inv{1}.inverse.M*D.inv{1}.inverse.U{1};
    megchans=D.indchantype('meg','good');
    legend_labels={};
    pial_clusters=invert_burst_tc_results.clusters;
    figure();
    hold all;
    for c_idx=1:length(pial_clusters)
        plot(mean(pial_clusters(c_idx).bc_f_diff,1));
        legend_labels{end+1}=num2str(c_idx);
    end
    legend(legend_labels);
    title(sprintf('subject %d', s_idx));
    
    source_tcs=[];
    bc_f_diffs=[];
    for c_idx=1:length(pial_clusters)
        verts=size(MU,1)/2+pial_clusters(c_idx).inv_verts;    
        source_tc=mean(MU(verts,:)*squeeze(D(megchans,:,:)),1);
        if mean(source_tc(26:31))>0
            source_tc=source_tc.*-1;
        end
        bc_f_diffs(end+1,:)=mean(pial_clusters(c_idx).bc_f_diff,1);
        source_tcs(end+1,:)=source_tc;        
    end
    mean_source_tcs=mean(source_tcs,1);
    subj_scaled_mean_source_tcs(s_idx,:)=mean_source_tcs./max(abs(mean_source_tcs));
    subj_mean_bc_fdiffs(s_idx,:)=mean(bc_f_diffs,1);
    
    figure();
    hold all
    plot(mean(bc_f_diffs,1),'k');
    plot(xlim(),[3 3],'k--');
    plot(xlim(),[-3 -3],'k--');

    title(sprintf('subject %d', s_idx));
        
    times=D.time.*1000;
    zero_time=times(round(length(times)/2))+.5*(times(round(length(times)/2+1))-times(round(length(times)/2)));
    times=times-zero_time;
end

[wout,lags]=woody(subj_scaled_mean_source_tcs',[],[],'woody','biased');
   
n_aligned=size(subj_scaled_mean_source_tcs,2)-abs(min(lags))-max(lags);
n_inv_aligned=size(subj_mean_bc_fdiffs,2)-abs(min(lags))-max(lags);
cn=1;
aligned_tcs=[];
aligned_mean_bc_f_diffs=[];
for subj_idx=1:size(subj_scaled_mean_source_tcs,1)
    start_idx=1+(lags(subj_idx)-min(lags));
    aligned_tcs(cn,:)=subj_scaled_mean_source_tcs(subj_idx,start_idx:start_idx+n_aligned-1);
    aligned_mean_bc_f_diffs(cn,:)=subj_mean_bc_fdiffs(subj_idx,start_idx:start_idx+n_inv_aligned-1);
    aligned_times=times(start_idx:start_idx+n_aligned-1);
    aligned_sliding_tc_times=sliding_tc_times(start_idx:start_idx+n_inv_aligned-1);
    cn=cn+1;
end
    

figure();
subplot(2,2,1);
hold on
plot(times,subj_scaled_mean_source_tcs');
plot(times,mean(subj_scaled_mean_source_tcs),'k','LineWidth',2);
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('Scaled source TCs');

subplot(2,2,2);
hold all
plot(aligned_times,aligned_tcs');
plot(aligned_times,mean(aligned_tcs),'k','LineWidth',2);
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('Aligned scaled source TCs');

subplot(2,2,3);
hold all
plot(sliding_tc_times,sum(subj_mean_bc_fdiffs,1),'LineWidth',2);
plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('Sum mean BC \Delta F');
xlabel('Time (ms)');

subplot(2,2,4);
hold all
plot(aligned_sliding_tc_times,sum(aligned_mean_bc_f_diffs,1),'LineWidth',2);
plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('Sum aligned mean BC \Delta F');
xlabel('Time (ms)');
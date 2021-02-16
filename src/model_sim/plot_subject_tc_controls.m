function plot_subject_tc_controls(subj_info, varargin)

defaults = struct('base_dir', 'JB_BUTTON_LOCKED_d3_ers');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

base_dir_parts=strsplit(params.base_dir,filesep);

all_dist=[];
all_prox=[];
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
end

ds_dist=ft_preproc_resample(all_dist,4.0000e+04,250,'resample');
ds_prox=ft_preproc_resample(all_prox,4.0000e+04,250,'resample');

n_sims=10;
n_sims_reversed=10;
n_sims_pial_only=10;
n_sims_white_only=10;

% Get sliding time window results
data_dir=fullfile('../../output/data',base_dir_parts{end},subj_info.subj_id,'model_sim');
% Get sliding time window results
for sim_idx=1:n_sims
    load(fullfile(data_dir, sprintf('invert_burst_tc_results_standard_%d.mat',sim_idx)));

    sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);

    % Get source tc
    D=spm_eeg_load(fullfile(data_dir, sprintf('grey_msim_standard_%d_grey_rcresp_TafdfC.mat',sim_idx)));
    MU=D.inv{1}.inverse.M*D.inv{1}.inverse.U{1};
    megchans=D.indchantype('meg','good');
    legend_labels={};
    pial_clusters=invert_burst_tc_results.clusters;

    source_tcs=[];
    f_diffs=[];
    for c_idx=1:length(pial_clusters)
        verts=size(MU,1)/2+pial_clusters(c_idx).inv_verts;    
        source_tc=mean(MU(verts,:)*squeeze(D(megchans,:,:)),1);
        if mean(source_tc(23:27))>0
            source_tc=source_tc.*-1;
        end
        source_tc=source_tc./max(abs(source_tc));
        f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        source_tcs(end+1,:)=source_tc;        
    end
    subj_mean_source_tcs(sim_idx,:)=mean(source_tcs,1);
    subj_mean_fdiffs(sim_idx,:)=mean(f_diffs,1);
end

% Get sliding time window results
for sim_idx=1:n_sims_reversed
    load(fullfile(data_dir, sprintf('invert_burst_tc_results_reversed_%d.mat',sim_idx)));

    sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);

    % Get source tc
    D=spm_eeg_load(fullfile(data_dir, sprintf('grey_msim_reversed_%d_grey_rcresp_TafdfC.mat',sim_idx)));
    MU=D.inv{1}.inverse.M*D.inv{1}.inverse.U{1};
    megchans=D.indchantype('meg','good');
    legend_labels={};
    pial_clusters=invert_burst_tc_results.clusters;

    source_tcs=[];
    f_diffs=[];
    for c_idx=1:length(pial_clusters)
        verts=size(MU,1)/2+pial_clusters(c_idx).inv_verts;    
        source_tc=mean(MU(verts,:)*squeeze(D(megchans,:,:)),1);
        if mean(source_tc(23:27))<0
            source_tc=source_tc.*-1;
        end
        source_tc=source_tc./max(abs(source_tc));
        f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        source_tcs(end+1,:)=source_tc;        
    end
    subj_mean_source_tcs_reversed(sim_idx,:)=mean(source_tcs,1);
    subj_mean_fdiffs_reversed(sim_idx,:)=mean(f_diffs,1);
end

% Get sliding time window results
for sim_idx=1:n_sims_pial_only
    load(fullfile(data_dir, sprintf('invert_burst_tc_results_pial_only_%d.mat',sim_idx)));

    sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);

    % Get source tc
    D=spm_eeg_load(fullfile(data_dir, sprintf('grey_msim_pial_only_%d_grey_rcresp_TafdfC.mat',sim_idx)));
    MU=D.inv{1}.inverse.M*D.inv{1}.inverse.U{1};
    megchans=D.indchantype('meg','good');
    legend_labels={};
    pial_clusters=invert_burst_tc_results.clusters;

    source_tcs=[];
    f_diffs=[];
    for c_idx=1:length(pial_clusters)
        verts=size(MU,1)/2+pial_clusters(c_idx).inv_verts;    
        source_tc=mean(MU(verts,:)*squeeze(D(megchans,:,:)),1);
        if mean(source_tc(23:27))>0
            source_tc=source_tc.*-1;
        end
        source_tc=source_tc./max(abs(source_tc));
        f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
        source_tcs(end+1,:)=source_tc;        
    end
    subj_mean_source_tcs_pial_only(sim_idx,:)=mean(source_tcs,1);
    subj_mean_fdiffs_pial_only(sim_idx,:)=mean(f_diffs,1);
end
% 
% % Get sliding time window results
% for sim_idx=1:n_sims_white_only
%     load(fullfile(data_dir, sprintf('invert_burst_tc_results_white_only_%d.mat',sim_idx)));
% 
%     sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
% 
%     % Get source tc
%     D=spm_eeg_load(fullfile(data_dir, sprintf('grey_msim_white_only_%d_grey_rcresp_TafdfC.mat',sim_idx)));
%     MU=D.inv{1}.inverse.M*D.inv{1}.inverse.U{1};
%     megchans=D.indchantype('meg','good');
%     legend_labels={};
%     pial_clusters=invert_burst_tc_results.clusters;
% 
%     source_tcs=[];
%     f_diffs=[];
%     for c_idx=1:length(pial_clusters)
%         verts=size(MU,1)/2+pial_clusters(c_idx).inv_verts;    
%         source_tc=mean(MU(verts,:)*squeeze(D(megchans,:,:)),1);
%         if mean(source_tc(23:27))<0
%             source_tc=source_tc.*-1;
%         end
%         source_tc=source_tc./max(abs(source_tc));
%         f_diffs(end+1,:)=mean(pial_clusters(c_idx).f_diff,1);
%         source_tcs(end+1,:)=source_tc;        
%     end
%     subj_mean_source_tcs_white_only(sim_idx,:)=mean(source_tcs,1);
%     subj_mean_fdiffs_white_only(sim_idx,:)=mean(f_diffs,1);
% end

deep_signal = zeros(length(D.time),size(ds_prox,1));
deep_signal(round(length(D.time)/2-size(ds_prox,2)/2)+1:round(length(D.time)/2+size(ds_prox,2)/2),:)=ds_prox';
superficial_signal = zeros(length(D.time),size(ds_dist,1));
superficial_signal(round(length(D.time)/2-size(ds_dist,2)/2)+1:round(length(D.time)/2+size(ds_dist,2)/2),:)=ds_dist';
deep_mag=max(abs(squeeze(mean(deep_signal,2))));
deep_signal=deep_signal./deep_mag;
superficial_mag=max(abs(squeeze(mean(superficial_signal,2))));
superficial_signal=superficial_signal./superficial_mag;
deep_moment=6;
superficial_moment=8;

times=D.time.*1000;
zero_time=times(round(length(times)/2))+.5*(times(round(length(times)/2+1))-times(round(length(times)/2)));
times=times-zero_time;

fig=figure();
set(fig,'PaperUnits','centimeters');
set(fig,'PaperPosition',[0 0 40.64 15.24]);
set(fig,'Position',[0 0 1120 420]);

subplot(3,2,1);
hold all;
plot(times, -deep_moment.*deep_signal','Color',[106 175 215]./255.0);
plot(times, -deep_moment.*mean(deep_signal,2),'Color','k');
ylabel('Dipole moment (nA)');
xlim(times([1 end]));
ylim([-10 1]);

subplot(3,2,3);
hold all
plot(times, -superficial_moment.*superficial_signal'-deep_moment.*deep_signal','Color',[.75 .75 .75]);
plot(times, -superficial_moment.*mean(superficial_signal,2)-deep_moment.*mean(deep_signal,2),'Color','k');
xlim(times([1 end]));

subplot(3,2,5);
hold all
plot(times, -superficial_moment.*superficial_signal','Color',[240 140 100]./255.0);
plot(times, -superficial_moment.*mean(superficial_signal,2),'Color','k');
xlim(times([1 end]));
ylim([-1 8]);
ylabel('Dipole moment (nA)');
xlabel('Time (ms)');

subplot(3,2,[2 4 6]);
hold all
plot(sliding_tc_times,sum(subj_mean_fdiffs_reversed),'k','LineWidth',2);
plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
%ylim([-75 75]);
ylabel('\Delta F');
xlabel('Time (ms)');


fig=figure();
hold on
plot(times,subj_mean_source_tcs_reversed');
shadedErrorBar(times,...
    mean(subj_mean_source_tcs_reversed,1),...
    std(subj_mean_source_tcs_reversed,[],1)./sqrt(n_sims_reversed),...
    'LineProps',{'Color','k','LineWidth',2});
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('Dipole moment (normalized)');


fig=figure();
set(fig,'PaperUnits','centimeters');
set(fig,'PaperPosition',[0 0 40.64 15.24]);
set(fig,'Position',[0 0 1120 420]);

subplot(3,2,1);
hold all;
plot(times, superficial_moment.*superficial_signal','Color',[240 140 100]./255.0);
plot(times, superficial_moment.*mean(superficial_signal,2),'Color','k');
ylabel('Dipole moment (nA)');
xlim(times([1 end]));
ylim([-8 1]);

subplot(3,2,3);
hold all
plot(times, superficial_moment.*superficial_signal','Color',[.75 .75 .75]);
plot(times, superficial_moment.*mean(superficial_signal,2),'Color','k');
xlim(times([1 end]));

subplot(3,2,5);
xlim(times([1 end]));
ylim([-8 1]);
ylabel('Dipole moment (nA)');
xlabel('Time (ms)');

subplot(3,2,[2 4 6]);
hold all
plot(sliding_tc_times,sum(subj_mean_fdiffs_pial_only),'k','LineWidth',2);
plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('\Delta F');
xlabel('Time (ms)');


fig=figure();
hold on
plot(times,subj_mean_source_tcs_pial_only');
shadedErrorBar(times,...
    mean(subj_mean_source_tcs_pial_only,1),...
    std(subj_mean_source_tcs_pial_only,[],1)./sqrt(n_sims_pial_only),...
    'LineProps',{'Color','k','LineWidth',2});
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('Dipole moment (normalized)');

% 
% fig=figure();
% set(fig,'PaperUnits','centimeters');
% set(fig,'PaperPosition',[0 0 40.64 15.24]);
% set(fig,'Position',[0 0 1120 420]);
% 
% subplot(3,2,1);
% hold all;
% ylabel('Dipole moment (nA)');
% xlim(times([1 end]));
% ylim([-8 1]);
% 
% subplot(3,2,3);
% hold all
% plot(times, deep_moment.*deep_signal','Color',[.75 .75 .75]);
% plot(times, deep_moment.*mean(deep_signal,2),'Color','k');
% xlim(times([1 end]));
% 
% subplot(3,2,5);
% hold all;
% plot(times, deep_moment.*deep_signal','Color',[106 175 215]./255.0);
% plot(times, deep_moment.*mean(deep_signal,2),'Color','k');
% xlim(times([1 end]));
% ylim([-1 10]);
% ylabel('Dipole moment (nA)');
% xlabel('Time (ms)');
% 
% subplot(3,2,[2 4 6]);
% hold all
% plot(sliding_tc_times,sum(subj_mean_fdiffs_white_only),'k','LineWidth',2);
% plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
% plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
% xlim([sliding_tc_times(1) sliding_tc_times(end)]);
% ylabel('\Delta F');
% xlabel('Time (ms)');
% 
% 
% fig=figure();
% hold on
% plot(times,subj_mean_source_tcs_white_only');
% shadedErrorBar(times,...
%     mean(subj_mean_source_tcs_white_only,1),...
%     std(subj_mean_source_tcs_white_only,[],1)./sqrt(n_sims_white_only),...
%     'LineProps',{'Color','k','LineWidth',2});
% xlim([sliding_tc_times(1) sliding_tc_times(end)]);
% ylabel('Dipole moment (normalized)');

function plot_subject_tc(subj_info, varargin)

defaults = struct('base_dir', 'JB_BUTTON_LOCKED_d3_ers');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

base_dir_parts=strsplit(params.base_dir,filesep);

n_sims=10;

output_data_dir=fullfile('../../output/data',params.base_dir,subj_info.subj_id,'simple_sim');

% Copy data file    
reg_file=fullfile(output_data_dir, 'grey_rcresp_TafdfC.mat');

% Load inverted data file
D=spm_eeg_load(reg_file);   

% time zero is midpoint of WOI
zero_time=D.time((length(D.time)-1)/2+1);

% Signal in deep layers
deep_signal_width=.025; % 25ms
deep_signal=exp(-((D.time-zero_time).^2)/(2*deep_signal_width^2));
deep_moment=4.5;

% Signal in superficial layers
superficial_signal_width=.01; % 10ms
superficial_signal=exp(-((D.time-zero_time).^2)/(2*superficial_signal_width^2));
superficial_moment=6;

% Get sliding time window results
data_dir=fullfile('../../output/data',base_dir_parts{end},subj_info.subj_id,'simple_sim');
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

times=D.time.*1000;
zero_time=times(round(length(times)/2))+.5*(times(round(length(times)/2+1))-times(round(length(times)/2)));
times=times-zero_time;

fig=figure();
set(fig,'PaperUnits','centimeters');
set(fig,'PaperPosition',[0 0 60.96 15.24]);
set(fig,'Position',[0 0 1680 420]);

subplot(3,3,1);
hold all;
plot(times, -superficial_moment.*superficial_signal','Color',[240 140 100]./255.0);
ylabel('Dipole moment (nA)');
xlim(times([1 end]));
%ylim([-8 1]);

subplot(3,3,4);
hold all;
plot(times, -superficial_moment.*superficial_signal'+deep_moment.*deep_signal','Color','k');
ylabel('Dipole moment (nA)');
xlim(times([1 end]));
%ylim([-6 8]);

subplot(3,3,7);
hold all
plot(times, deep_moment.*deep_signal','Color',[106 175 215]./255.0);
xlim(times([1 end]));
%ylim([-1 10]);
ylabel('Dipole moment (nA)');
xlabel('Time (ms)');

subplot(3,3,[2 5 8]);
plot(times,subj_mean_source_tcs');
shadedErrorBar(times,...
    mean(subj_mean_source_tcs,1),...
    std(subj_mean_source_tcs,[],1)./sqrt(n_sims),...
    'LineProps',{'Color','k','LineWidth',2});
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
ylabel('Dipole moment (normalized)');

subplot(3,3,[3 6 9]);
hold all
plot(sliding_tc_times,sum(subj_mean_fdiffs),'k','LineWidth',2);
plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
xlim([sliding_tc_times(1) sliding_tc_times(end)]);
%ylim([-75 75]);
ylabel('\Delta F');
xlabel('Time (ms)');





% fig=figure();
% set(fig,'PaperUnits','centimeters');
% set(fig,'PaperPosition',[0 0 40.64 15.24]);
% set(fig,'Position',[0 0 1120 420]);
% 
% subplot(3,2,1);
% plot(times, superficial_signal.*-superficial_moment,'Color',[240 140 100]./255.0, 'LineWidth',2);
% ylabel('Dipole moment (nA)');
% xlim(times([1 end]));
% 
% subplot(3,2,3);
% hold on
% plot(times,subj_mean_source_tcs');
% shadedErrorBar(times,...
%     mean(subj_mean_source_tcs,1),...
%     std(subj_mean_source_tcs,[],1)./sqrt(n_sims),...
%     'LineProps',{'Color','k','LineWidth',2});
% xlim([sliding_tc_times(1) sliding_tc_times(end)]);
% ylabel('Dipole moment (normalized)');
% 
% subplot(3,2,5);
% plot(times, deep_signal.*deep_moment,'Color',[106 175 215]./255.0, 'LineWidth',2);
% xlim(times([1 end]));
% ylabel('Dipole moment (nA)');
% xlabel('Time (ms)');
% 
% subplot(3,2,[2 4 6]);
% hold all
% plot(sliding_tc_times,sum(subj_mean_fdiffs),'k','LineWidth',2);
% plot([sliding_tc_times(1) sliding_tc_times(end)],[3 3],'k--');
% plot([sliding_tc_times(1) sliding_tc_times(end)],[-3 -3],'k--');
% xlim([sliding_tc_times(1) sliding_tc_times(end)]);
% ylim([-500 500]);
% ylabel('\Delta F');
% xlabel('Time (ms)');

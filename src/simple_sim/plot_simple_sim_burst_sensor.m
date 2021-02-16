function plot_simple_sim_burst_sensor(subj_info, type, idx, varargin)

defaults = struct('base_dir','../../output/data/JB_BUTTON_LOCKED_d3_ers',...
    'channel', 'MLP34', 'n_temp_modes', 4, 'data_type', 'mean_evoked', 'flipped', false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

base_dir_parts=strsplit(params.base_dir,'/');

% Where to save plots
plot_dir=fullfile('../../output/figures',base_dir_parts{end},...
    subj_info.subj_id, 'simple_sim');
mkdir(plot_dir);

data_dir=fullfile(params.base_dir,subj_info.subj_id,'simple_sim');
reg_file=fullfile(data_dir, sprintf('sim_%s_%d_grey_rcresp_TafdfC.mat',type,idx));
D=spm_eeg_load(reg_file);

times=D.time;
zero_time=times(1)+(times(end)-times(1))/2;

ntrials=size(D,3);
times=(times-zero_time)*1000;

meg_ch_idx=D.indchantype('MEG');

% Position of each meg channel
ch_pos=D.coor2D(meg_ch_idx);
% Label for each meg channel
ch_labels=D.chanlabels(meg_ch_idx);
trials=setdiff([1:ntrials],D.badtrials);

data_vals=D(meg_ch_idx,:,trials);

if length(params.channel)==0
    mean_trial=mean(data_vals,3);
    mag_trial=max(mean_trial,[],2)-min(mean_trial,[],2);
    chan_idx=find(mag_trial==max(mag_trial));
else
    chan_idx=find(strcmp(D.chanlabels(meg_ch_idx),params.channel));
end
disp(sprintf('Using channel %s', D.chanlabels{meg_ch_idx(chan_idx)}));

mean_tc=squeeze(mean(data_vals(chan_idx,:,:),3));
stderr_tc=squeeze(std(data_vals(chan_idx,:,:),[],3))./sqrt(size(data_vals,3));

if params.flipped
    data_vals=-1.*data_vals;
    mean_tc=squeeze(mean(data_vals(chan_idx,:,:),3));
    stderr_tc=squeeze(std(data_vals(chan_idx,:,:),[],3))./sqrt(size(data_vals,3));
end
peak_time=find(mean_tc==min(mean_tc),1);
left_min=find(mean_tc(1:peak_time-1)==max(mean_tc(1:peak_time-1)));
right_min=peak_time+find(mean_tc(peak_time+1:end)==max(mean_tc(peak_time+1:end)));

mean_scalp_vals=squeeze(mean(data_vals(:,[left_min peak_time right_min],:),3));
    
fig=figure('Position',[1 1 1600 400],'PaperUnits','points',...
    'PaperPosition',[1 1 1600 400],'PaperPositionMode','manual');
for i=1:size(mean_scalp_vals,2)
    ax=subplot(2,4,[i i+4]);
    in.f=fig;
    in.ParentAxes=ax;
    in.noButtons=true;
    in.type='MEG';
    in.min=min(mean_scalp_vals(:));
    %in.min=-300;
    in.max=max(mean_scalp_vals(:));
    %in.max=200;
    [ZI,f]=spm_eeg_plotScalpData(mean_scalp_vals(:,i),ch_pos,ch_labels,in);
    children=get(f,'Children');
    d=get(children(2),'UserData');
    set(d.ht(chan_idx),'visible','on');
    xdata=get(d.hp,'XData');
    ydata=get(d.hp,'YData');
    set(d.hp,'XData',xdata(chan_idx));
    set(d.hp,'YData',ydata(chan_idx));
end

 
subplot(2,4,4);
hold all;
for i=1:length(trials)
    trial_idx=trials(i);
    plot(times, squeeze(data_vals(chan_idx,:,trial_idx)));
end
%ylim([-1250 1250]);
yl=ylim();
plot([times(left_min) times(left_min)],yl,'r','LineWidth',2);
plot([times(peak_time) times(peak_time)],yl,'r','LineWidth',2);
plot([times(right_min) times(right_min)],yl,'r','LineWidth',2);
hold off;
xlim([times(1) times(end)]);
ylabel('Field Intensitiy (fT)');

subplot(2,4,8);
hold on;
shadedErrorBar(times, mean_tc, stderr_tc, 'lineprops','b');
%ylim([-250 250]);
yl=ylim();
plot([times(left_min) times(left_min)],yl,'r','LineWidth',2);
plot([times(peak_time) times(peak_time)],yl,'r','LineWidth',2);
plot([times(right_min) times(right_min)],yl,'r','LineWidth',2);
hold off;
xlim([times(1) times(end)]);
xlabel('Time (ms)');
ylabel('Field Intensitiy (fT)');
% 
% saveas(fig, fullfile(plot_dir, 'sensor_data.png'), 'png');
% saveas(fig, fullfile(plot_dir, 'sensor_data.eps'), 'epsc');
% saveas(fig, fullfile(plot_dir, 'sensor_data.fig'), 'fig');
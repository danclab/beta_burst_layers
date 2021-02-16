function plot_subjects_burst_sensor(subjects, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',  
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

erd_dir='../output/data/JB_BUTTON_LOCKED_d3_erd2_2_28_20';
ers_dir='../output/data/JB_BUTTON_LOCKED_d3_ers';

spm('defaults','eeg');

erd_subj_tcs=[];
ers_subj_tcs=[];
erd_flipped=[0 0 0 0 0 0 1 1];
ers_flipped=[0 1 0 0 0 1 1 0];

subj_ids={};
for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);
    subj_dir=fullfile(erd_dir, subj_info.subj_id);
    D=spm_eeg_load(fullfile(subj_dir, 'rcresp_TafdfC.mat'));

    subj_ids{end+1}=num2str(s_idx);
    
    times=D.time;
    zero_time=times(1)+(times(end)-times(1))/2;

    ntrials=size(D,3);
    times=(times-zero_time)*1000;

    meg_ch_idx=D.indchantype('MEG');

    trials=setdiff([1:ntrials],D.badtrials);

    data_vals=D(meg_ch_idx,:,trials);

    %mean_trial=mean(data_vals,3);
    %mag_trial=max(mean_trial,[],2)-min(mean_trial,[],2);
    %chan_idx=find(mag_trial==max(mag_trial));
    chan_idx=find(strcmp(D.chanlabels(meg_ch_idx),'MLP34'));
    disp(sprintf('Using channel %s', D.chanlabels{meg_ch_idx(chan_idx)}));

    mean_tc=squeeze(mean(data_vals(chan_idx,:,:),3));
    
    if erd_flipped(s_idx)
        data_vals=-1.*data_vals;
        mean_tc=squeeze(mean(data_vals(chan_idx,:,:),3));
    end
    
    erd_subj_tcs(s_idx,:)=mean_tc;
end

[wout,lags]=woody(erd_subj_tcs',[],[],'woody','biased');   
n_aligned=size(erd_subj_tcs,2)-abs(min(lags))-max(lags);

aligned_erd_tcs=[];
for s_idx=1:size(erd_subj_tcs,1)
    start_idx=1+(lags(s_idx)-min(lags));
    aligned_erd_tcs(s_idx,:)=erd_subj_tcs(s_idx,start_idx:start_idx+n_aligned-1);    
    aligned_erd_times=times(start_idx:start_idx+n_aligned-1);
end
zero_time=aligned_erd_times(round(length(aligned_erd_times)/2))+.5*(aligned_erd_times(round(length(aligned_erd_times)/2+1))-aligned_erd_times(round(length(aligned_erd_times)/2)));
aligned_erd_times=aligned_erd_times-zero_time;

for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);
    subj_dir=fullfile(ers_dir, subj_info.subj_id);
    D=spm_eeg_load(fullfile(subj_dir, 'rcresp_TafdfC.mat'));

    times=D.time;
    zero_time=times(1)+(times(end)-times(1))/2;

    ntrials=size(D,3);
    times=(times-zero_time)*1000;

    meg_ch_idx=D.indchantype('MEG');

    trials=setdiff([1:ntrials],D.badtrials);

    data_vals=D(meg_ch_idx,:,trials);

%     mean_trial=mean(data_vals,3);
%     mag_trial=max(mean_trial,[],2)-min(mean_trial,[],2);
%     chan_idx=find(mag_trial==max(mag_trial));
    chan_idx=find(strcmp(D.chanlabels(meg_ch_idx),'MLP34'));
    disp(sprintf('Using channel %s', D.chanlabels{meg_ch_idx(chan_idx)}));

    mean_tc=squeeze(mean(data_vals(chan_idx,:,:),3));
    
    if ers_flipped(s_idx)
        data_vals=-1.*data_vals;
        mean_tc=squeeze(mean(data_vals(chan_idx,:,:),3));
    end
    
    ers_subj_tcs(s_idx,:)=mean_tc;
end

[wout,lags]=woody(ers_subj_tcs',[],[],'woody','biased');   
n_aligned=size(ers_subj_tcs,2)-abs(min(lags))-max(lags);

aligned_ers_tcs=[];
for s_idx=1:size(ers_subj_tcs,1)
    start_idx=1+(lags(s_idx)-min(lags));
    aligned_ers_tcs(s_idx,:)=ers_subj_tcs(s_idx,start_idx:start_idx+n_aligned-1);    
    aligned_ers_times=times(start_idx:start_idx+n_aligned-1);
end
zero_time=aligned_ers_times(round(length(aligned_ers_times)/2))+.5*(aligned_ers_times(round(length(aligned_ers_times)/2+1))-aligned_ers_times(round(length(aligned_ers_times)/2)));
aligned_ers_times=aligned_ers_times-zero_time;


fig=figure();
subplot(2,2,1);
plot(aligned_erd_times, aligned_erd_tcs');
legend(subj_ids);
xlim(aligned_erd_times([1 end]));
xlabel('Time (ms)');
ylabel('Field Intensitiy (fT)');
subplot(2,2,2);
plot(aligned_ers_times, aligned_ers_tcs');
xlim(aligned_ers_times([1 end]));
xlabel('Time (ms)');
subplot(2,2,3);
shadedErrorBar(aligned_erd_times, mean(aligned_erd_tcs), std(aligned_erd_tcs)./sqrt(length(subjects)), 'lineprops','b');
xlim(aligned_erd_times([1 end]));
xlabel('Time (ms)');
ylabel('Field Intensitiy (fT)');
subplot(2,2,4);
shadedErrorBar(aligned_ers_times, mean(aligned_ers_tcs), std(aligned_ers_tcs)./sqrt(length(subjects)), 'lineprops','b');
xlim(aligned_ers_times([1 end]));
xlabel('Time (ms)');

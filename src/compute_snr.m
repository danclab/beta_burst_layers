function compute_snr(subjects, data_dir)

spm('defaults','eeg');

subj_snrs=[];
for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);
    subj_dir=fullfile('../output/data',data_dir, subj_info.subj_id);
    D=spm_eeg_load(fullfile(subj_dir, 'rcresp_TafdfC.mat'));

    % Get good trials
    ntrials=size(D,3);
    trials=setdiff([1:ntrials],D.badtrials);

    % Position of each meg channel
    meg_ch_idx=D.indchantype('MEG');
    
    % Data from MEG channels and good trials
    data_vals=D(meg_ch_idx,:,trials);

    % Average over trials (bursts)
    mean_data=squeeze(mean(data_vals,3));
    
    % Subtract average from data to get noise
    noise_data=data_vals-repmat(mean_data,1,1,length(trials));
    
    % Standard deviation of signal (over time)
    allchanstd=std(mean_data,[],2);
    % Mean standard deviation of signal over channels
    meanrmssignal=mean(allchanstd);
    
    % Average noise over trials (bursts)
    mean_noise_data=squeeze(mean(noise_data,3));
    
    % Standard deviation of noise (over time)
    allchanstd=std(mean_noise_data,[],2);
    % Mean standard deviation of noise over channels
    meanrmsnoise=mean(allchanstd);
    
    snr=-20*log10(meanrmsnoise/meanrmssignal);
    subj_snrs(s_idx)=snr;
end

mean_snr=mean(subj_snrs);
se_snr=std(subj_snrs)/sqrt(length(subjects));
disp(sprintf('SNR: M=%.3f, SE=%.3f dB', mean_snr, se_snr));
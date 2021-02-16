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

srate=40000;
time=linspace(-90,90,srate*.18+1);
deep_signal = zeros(length(time),size(all_prox,1));
deep_signal(round(length(time)/2-size(all_prox,2)/2)+1:round(length(time)/2+size(all_prox,2)/2),:)=all_prox';
superficial_signal = zeros(length(time),size(all_dist,1));
superficial_signal(round(length(time)/2-size(all_dist,2)/2)+1:round(length(time)/2+size(all_dist,2)/2),:)=all_dist';

cum_dip=mean(superficial_signal+deep_signal,2);

min_freq = 0;
max_freq = 60;
num_frex = 100;

% frequencies vector
frex = linspace(min_freq, max_freq, num_frex);

%% wavelet cycles - variable : min 4 max 10
range_cycles = [3 10];
cylvec = logspace(log10(range_cycles(1)), log10(range_cycles(end)), num_frex)./ (2*pi*frex);

%% wavelet parameters
wavtime = -1:1/srate:1; % length of wavelet
half_wave = (length(wavtime)-1)/2;
        
%% FFT parameters
nWave = length(wavtime);
nData = length(time);
nConv = nWave + nData - 1;
        
tf_data=zeros(length(frex), nData);
            
%% Run wavelet convolution
for fi=1:length(frex) % loop through all frequencies
                    
    %% Create wavelate
    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*cylvec(fi)^2));
    waveletX = fft(wavelet, nConv); % fft of wavelet
    waveletX = waveletX ./ max(waveletX); % normalize fft of wavelet
                    
    data = fft(cum_dip', nConv);
                        
    %% run convolution
    data_conv = ifft(waveletX .* data);
    data_conv = data_conv(half_wave+1:end-half_wave);
                        
    %% compute power
    tf_data(fi,:) = abs(data_conv).^2;
                        
end

figure();
subplot(2,1,1);
plot(time,cum_dip,'k');
xlim([time(1) time(end)]);
ylabel('Dipole (nAm');
subplot(2,1,2);
imagesc(time,frex,tf_data);
set(gca,'ydir','normal');
curr_pos=get(gca,'position');
colorbar();
set(gca,'position',curr_pos);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
function plot_results_f(subj_info, data_dir, type)

spm('defaults','eeg');

g=gifti(fullfile('../../../output/data', data_dir, subj_info.subj_id,'pial_mrcresp_TafdfC_1_t-2000_-1800_f_1.gii'));
[max_val,sim_vertex]=max(g.cdata(:));
% Load meshes
subj_surf_dir=fullfile('../../../data/surf', sprintf('%s-synth',...
    subj_info.subj_id),'surf');
pial_fname=fullfile(subj_surf_dir,'pial.ds.link_vector.gii');
pial=gifti(pial_fname);
sim_coord=pial.vertices(sim_vertex,:);

snrs=[-50 -45 -40 -35 -30 -25 -20];
n_sims=50;

base_dir_parts=strsplit(data_dir,filesep);
base_data_dir=fullfile('../../../output/data',base_dir_parts{end},subj_info.subj_id,'simple_sim','snr_-100');
reg_file=fullfile(base_data_dir, 'grey_rcresp_TafdfC.mat');

% Load inverted data file
D=spm_eeg_load(reg_file);   

all_f_diffs=[];
tail=[];
peak=[];

for s_idx=1:length(snrs)
    snr=snrs(s_idx);
    base_data_dir=fullfile('../../../output/data',base_dir_parts{end},subj_info.subj_id,'simple_sim',sprintf('snr_%d',snr));
    for idx=1:n_sims
        load(fullfile(base_data_dir, sprintf('invert_burst_tc_results_%s_%d.mat',type,idx)));
        sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
        f_diffs=[];
        for j=1:length(invert_burst_tc_results.clusters)
            f_diffs(j,:)=mean(invert_burst_tc_results.clusters(j).f_diff,1);
        end
        f_diff=mean(f_diffs,1);
        %figure();
        %plot(f_diff);
        all_f_diffs(s_idx,idx,:)=f_diff;

        peak_idx=find((sliding_tc_times>-10) & (sliding_tc_times<10));
        peak(s_idx,idx)=max(f_diff(peak_idx));
        
        tail_idx=find(((sliding_tc_times>-40) & (sliding_tc_times<-10)) | ((sliding_tc_times>10) & (sliding_tc_times<40)));
        tail(s_idx,idx)=min(f_diff(tail_idx));
        
    end
    
    % Normalized f_diff
%     sum_f_diff=squeeze(sum(all_f_diffs(s_idx,:,:),2));
%     norm_f_diff=sum_f_diff./max(abs(sum_f_diff));
%     n_aligned=36;%length(norm_f_diff)-abs(min(lags))-max(lags);
%     start_idx=6;%1+(lags(1)-min(lags));
%     aligned_prediction=prediction(start_idx:start_idx+n_aligned-1);
%     start_idx=1;%1+(lags(2)-min(lags));
%     aligned_f_diff=norm_f_diff(start_idx:start_idx+n_aligned-1);
%     dist=sqrt(sum((aligned_prediction'-aligned_f_diff).^2));
%     rsqr=rsquare(aligned_prediction',aligned_f_diff);
%     
%     all_dists(s_idx)=dist;
%     all_rsqs(s_idx)=rsqr;
end      

figure();
for s_idx=1:length(snrs)
    subplot(2,ceil(length(snrs)/2),s_idx);
    hold all
    shadedErrorBar(sliding_tc_times,squeeze(mean(all_f_diffs(s_idx,:,:),2)),squeeze(std(all_f_diffs(s_idx,:,:),[],2))./sqrt(n_sims),'LineProps',{'Color','b'});
    plot(sliding_tc_times([1 end]), [-3 -3], 'k--');
    plot(sliding_tc_times([1 end]), [3 3], 'k--');
    xlim(sliding_tc_times([1 end]));
    ylim([-60 40]);
    %ylim([-450 300]);
    ylabel('\Delta F');
    title(sprintf('%d dB', snrs(s_idx)));
end

figure();
hold all;
shadedErrorBar(snrs,squeeze(mean(tail,2)),squeeze(std(tail,[],2))/sqrt(n_sims),'LineProps',{'Color','b'});
shadedErrorBar(snrs,squeeze(mean(peak,2)),squeeze(std(peak,[],2))/sqrt(n_sims),'LineProps',{'Color','b'});
plot(xlim(), [-3 -3], 'k--');
plot(xlim(), [3 3], 'k--');
xlabel('SNR (dB)');
ylabel('\Delta F');
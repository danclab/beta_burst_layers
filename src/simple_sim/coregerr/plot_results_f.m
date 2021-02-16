function plot_results_f(subj_info, data_dir, type)

spm('defaults','eeg');

coregerrs=[0 .25 .5 .75 1 1.5 2 3 4 5];

n_sims=50;

base_dir_parts=strsplit(data_dir,filesep);
base_data_dir=fullfile('../../../output/data',base_dir_parts{end},subj_info.subj_id,'simple_sim',sprintf('coregerr_%d',coregerrs(1)));
reg_file=fullfile(base_data_dir, 'grey_rcresp_TafdfC.mat');

% Load inverted data file
D=spm_eeg_load(reg_file);   

all_f_diffs=[];
peak=[];
tail=[];

for c_idx=1:length(coregerrs)
    coregerr=coregerrs(c_idx);
    base_data_dir=fullfile('../../../output/data',base_dir_parts{end},subj_info.subj_id,'simple_sim',sprintf('coregerr_%d',coregerr));
    for idx=1:n_sims
        load(fullfile(base_data_dir, sprintf('invert_burst_tc_results_%s_%d.mat',type,idx)));
        sliding_tc_times=invert_burst_tc_results.times(invert_burst_tc_results.left_idx:invert_burst_tc_results.right_idx);
        f_diffs=[];
        for j=1:length(invert_burst_tc_results.clusters)
            f_diffs(j,:)=mean(invert_burst_tc_results.clusters(j).f_diff,1);
        end
        f_diff=mean(f_diffs,1);
%         figure();
%         plot(sliding_tc_times,f_diff);
        all_f_diffs(c_idx,idx,:)=f_diff;
        
        peak_idx=find((sliding_tc_times>-10) & (sliding_tc_times<10));
        peak(c_idx,idx)=max(f_diff(peak_idx));
        
        tail_idx=find(((sliding_tc_times>-40) & (sliding_tc_times<-10)) | ((sliding_tc_times>10) & (sliding_tc_times<40)));
        tail(c_idx,idx)=min(f_diff(tail_idx));
    end
end      

figure();
for c_idx=1:length(coregerrs)
    subplot(2,ceil(length(coregerrs)/2),c_idx);
    hold all
    shadedErrorBar(sliding_tc_times,squeeze(mean(all_f_diffs(c_idx,:,:),2)),squeeze(std(all_f_diffs(c_idx,:,:),[],2))./sqrt(n_sims),'LineProps',{'Color','b'});
    plot(sliding_tc_times([1 end]), [-3 -3], 'k--');
    plot(sliding_tc_times([1 end]), [3 3], 'k--');
    xlim(sliding_tc_times([1 end]));
    ylim([-60 40]);
    ylabel('\Delta F');
    title(sprintf('%d mm/deg', coregerrs(c_idx)));
end

figure();
hold all
shadedErrorBar(coregerrs,squeeze(mean(tail,2)),squeeze(std(tail,[],2))/sqrt(n_sims),'LineProps',{'Color','b'});
shadedErrorBar(coregerrs,squeeze(mean(peak,2)),squeeze(std(peak,[],2))/sqrt(n_sims),'LineProps',{'Color','b'});
plot(xlim(), [-3 -3], 'k--');
plot(xlim(), [3 3], 'k--');
xlabel('Coreg error (mm/deg)');
ylabel('\Delta F');
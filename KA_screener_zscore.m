function [out] = KA_screener_zscore(cell_to_process)

if ~exist('Zone_plots', 'dir')
    mkdir('Zone_plots')
end
c_ord = linspecer(5); % nice colours.
%% load data

%%% load the events. 
evt = LoadEvents([]);
% add in check for multiple recording periods.  Some seem to have a pre and post recoding.
s_rec_idx = find(contains(evt.label, 'Starting Record')); 
e_rec_idx = find(contains(evt.label, 'Stopping Record')); 

nRec = length(evt.t{s_rec_idx}); 
if nRec >1
    for iR = nRec:-1:1
        rec_dur(iR) = evt.t{e_rec_idx}(iR) - evt.t{s_rec_idx}(iR);
    end
    [~, task_rec_idx] = max(rec_dur); 
    
else
   task_rec_idx = 1;  
   rec_dur = evt.t{e_rec_idx}(task_rec_idx) - evt.t{s_rec_idx}(task_rec_idx);
end

% check for sessions that are too short. If less than 15mins skip. 
if max(rec_dur)/60 < 10
    fprintf('<strong>Minumum recording duration (15mins) not met: %2.1fmins</strong>\n', max(rec_dur)/60);
    out = 'too short';  
    return
end


% relabel the events file for arm/feeders
evt.label{3} = 'NorthPellets';
evt.label{4} = 'EastPellets';
evt.label{5} = 'SouthPellets';
evt.label{6} = 'WestPellets';
evt.label{7} = 'TotalPellets';

% change the order to fit the ZoneIn
evt_og = evt; % temp version
evt.t{4} = evt_og.t{6};  %swap east and west
evt.t{6} = evt_og.t{4};  %swap east and west

evt.label{4} = evt_og.label{6};  %swap east and west
evt.label{6} = evt_og.label{4};  %swap east and west



%%% load the position %%%
cfg_pos.convFact = [560/142 480/142]; 
out.pos = LoadPos(cfg_pos);
out.pos = restrict(out.pos, evt.t{s_rec_idx}(task_rec_idx), evt.t{e_rec_idx}(task_rec_idx)); % restrict position to time on track. 

out.velo = getLinSpd([], out.pos); 
out.velo.data = interp1(out.velo.tvec,out.velo.data(1,:),out.pos.tvec,'linear');

% smooth speed over 0.5 seconds
out.velo_smooth = out.velo; 
% out.velo_smooth.data = smooth(out.velo.data, round(1/mode(diff(out.pos.tvec)))*.5)'; % smooth with gaussian 


%%% load the spikes %%%
cfg = [];
cfg.getTTnumbers = 0;
cfg.fc = {cell_to_process};
out.S = LoadSpikes(cfg);

out.S = restrict(out.S, evt.t{s_rec_idx}(task_rec_idx), evt.t{e_rec_idx}(task_rec_idx)); % restrict the spikes recording periods.avoids odd thing where spike trains contains zeros.  MClust issue?

if length(out.S.t{1}) / (out.pos.tvec(end) - out.pos.tvec(1)) < .2
    fprintf('<strong>Minumum mean spike rate (0.5Hz) not met: %2.1fhz</strong>\n', length(out.S.t{1}) / (out.pos.tvec(end) - out.pos.tvec(1)));
    out =  length(out.S.t{1}) / (out.pos.tvec(end) - out.pos.tvec(1));  
    return
end

% interp to match position and spikes. 

spk_x = interp1(out.pos.tvec,out.pos.data(1,:),out.S.t{1},'linear');
spk_y = interp1(out.pos.tvec,out.pos.data(2,:),out.S.t{1},'linear');



%%% load PM control script vars %%%
PM_dir = dir('PM*.mat');
load(PM_dir.name)

Zone_names = {'North', 'West', 'South', 'East', 'All'};

for ii = length(EnteringZoneTime):-1:1
    Zone_cord(ii,:) = c_ord(ZoneIn(ii),:);
end

zone_idx = nearest_idx3(EnteringZoneTime/1000000, out.pos.tvec);
%% Summary of spiking in time and space.
figure(101)
subplot(2,3,1)
hold on
plot(out.pos.data(1,:), out.pos.data(2,:), 'color', [.7 .7 .7])

plot(spk_x,spk_y, '.r')
set(gca,'YDir','reverse')
axis off

subplot(2,3,2)
hold on
plot(out.pos.data(1,:), out.pos.data(2,:), 'color', [.7 .7 .7])
gscatter(nan(1,4), nan(1,4),1:4, c_ord(1:4,:)); % space filler to get the colors in the legend to work. 
for ii = 1:length(zone_idx) % hack to get the colours to work properly. 
    gscatter(out.pos.data(1,zone_idx(ii)), out.pos.data(2,zone_idx(ii)), ZoneIn(ii), Zone_cord(ii,:), 'o', 3)
% gscatter(out.pos.data(1,zone_idx), out.pos.data(2,zone_idx), ZoneIn, Zone_cord, 'o+x.', 3)
end
axis off
set(gca,'YDir','reverse')
legend(['pos', Zone_names(1:4)])
% add in waveforms if they exist.  WIP. 
% for ii = 1:size(mWV,2)
%     plot(500:(160/32):655,(mWV(:,ii)*.010)+400, 'color',c_ord(ii,:), 'linewidth', 2)
% end

text(max(out.pos.data(1,:))*.1, max(out.pos.data(2,:))*.9, {['Firing rate = ' num2str(length(out.S.t{1})/(out.pos.tvec(end) - out.pos.tvec(1)),3) 'Hz'] ;['Mode ISI = ' num2str((mode(diff(out.S.t{1})))*10000, 5) 'ms'] ; ['Median ISI = ' num2str((median(diff(out.S.t{1})))*10000, 5) 'ms']})

subplot(2,3,3)
histogram(diff(out.S.t{1})*10000, 0:20:2000);
xlabel('ISI ms')
ylabel('Spike count')


subplot(2,3, 4:6)
cfg_rast = [];

hold on
for ii = unique(ZoneIn)
    F_idx = find(ZoneIn == ii);
    e_idx = nearest_idx3(EnteringZoneTime(F_idx)/1000000, out.pos.tvec);
    for jj = 1:length(e_idx)
        rectangle('position', [out.pos.tvec(e_idx(jj)) - 2, .5, 4, 1],'edgecolor', c_ord(ii,:), 'facecolor',[c_ord(ii,:) 0.2])
        text(out.pos.tvec(e_idx(jj)), 1.55, Zone_names{ii}(1), 'color', c_ord(ii,:), 'HorizontalAlignment', 'center')
    end
end

PlotSpikeRaster2(cfg_rast, out.S);
xlim([out.pos.tvec(nearest_idx3(EnteringZoneTime(1)/1000000, out.pos.tvec))-2 out.pos.tvec(nearest_idx3(EnteringZoneTime(end)/1000000, out.pos.tvec))+2])


%% get the mean and std for the data to zscore the rate

cfg_z.binsize = 0.1;
cfg_z.gauss_window = 1;
cfg_z.gauss_sd = 0.02;

% set up gau kernal
gauss_window = cfg_z.gauss_window./cfg_z.binsize; % 1 second window
gauss_SD = cfg_z.gauss_sd./cfg_z.binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./cfg_z.binsize; % normalize by binsize

tbin_edges = evt.t{1}(1):cfg_z.binsize:evt.t{2}(end);




spk_count = histc(out.S.t{1},tbin_edges);
S_gau_sdf = conv2(spk_count(1:end-1),gk,'same'); % convolve with gaussian window
gau_mean = mean(S_gau_sdf);
gau_std = std(S_gau_sdf);




%% spike Peth
% prepare the PETHS
cfg_peth = [];
cfg_peth.window = [-2.8 2.8];
%     cfg_peth.dt = 0.0025; % fine resulition for regular PETH.
cfg_peth.dt = cfg_z.binsize; % wider bins for 'Lap' plot later.
% cfg_peth.plot_type = 'zscore';
cfg_peth.plot_type = 'raw';


%% generate a shuffle distribution 
%     figure(9999)
    clip_idx = ceil(cfg_peth.window(2) * 1/mode(diff(out.pos.tvec))); % get clip index to prevent sampling too close to the statt 
    rand_tidx = datasample(out.pos.tvec(clip_idx:end-clip_idx),10000); % get nShuff timestamps a pseudo feeder times. 
    
    
%     rand_tidx = sort(rand_tidx);
    
    [shuff.outputS, shuff.outputIT, shuff.outputGau, shuff.mean_S_gau, shuff.pre_stim_means, shuff.post_stim_means,~, ~, shuff.Z{ii}] = SpikePETH(cfg_peth, out.S, rand_tidx);
    
        std_S_gau = std(shuff.outputGau,[],2); 

    
%         hold on
%         yyaxis left
%     plot(shuff.outputIT, shuff.mean_S_gau, '-', 'color', [.2 .2 .2])
%     plot(shuff.outputIT, shuff.mean_S_gau +std_S_gau*1.96, '-.', 'color', [.2 .2 .2])
%     plot(shuff.outputIT, shuff.mean_S_gau -std_S_gau*1.96, '-.', 'color', [.2 .2 .2])
%     
%     cfg_peth.z_mean = mean(shuff.mean_S_gau); 
%     cfg_peth.z_std = mean(std_S_gau);

% 
    % re-define the 'pre' and 'post' means
    pre_idx = nearest_idx3([-2.8 -1.8], shuff.outputIT);
    act_idx = nearest_idx3([-0.5 0.5], shuff.outputIT);  
    rew_idx = nearest_idx3([1 2], shuff.outputIT);  
    
    
    shuff.pre_stim_means = nanmean(shuff.outputGau(pre_idx(1):pre_idx(2),:),1); 
    shuff.act_stim_means = nanmean(shuff.outputGau(act_idx(1):act_idx(2),:),1); 
    shuff.rew_stim_means = nanmean(shuff.outputGau(rew_idx(1):rew_idx(2),:),1);
    
    z_mean_act_pre = mean(shuff.act_stim_means - shuff.pre_stim_means);
    z_std_act_pre = std(shuff.act_stim_means - shuff.pre_stim_means);
    
    z_mean_rew_pre = mean(shuff.rew_stim_means - shuff.pre_stim_means);
    z_std_rew_pre = std(shuff.rew_stim_means - shuff.pre_stim_means);   
    
    z_mean_rew_act = mean(shuff.rew_stim_means - shuff.act_stim_means);
    z_std_rew_act = std(shuff.rew_stim_means - shuff.act_stim_means);   

    
%     
%     figure(919)
%     subplot(3,1,1)
%     hist(shuff.pre_stim_means - shuff.act_stim_means, 50);
%     ylabel('pre - act');
%     vline([-std([shuff.pre_stim_means - shuff.act_stim_means])*1.96, std([shuff.pre_stim_means - shuff.act_stim_means])*1.96])
%     xlim([-20 20]);
%     
%     subplot(3,1,2)
%     hist(shuff.rew_stim_means - shuff.act_stim_means, 50);
%     ylabel('act- rew');
%     vline([-std([shuff.rew_stim_means - shuff.act_stim_means])*1.96, std([shuff.rew_stim_means - shuff.act_stim_means])*1.96])
%     xlim([-20 20]);
%     
%     subplot(3,1,3)
%     hist(shuff.rew_stim_means - shuff.pre_stim_means, 50);
%     ylabel('rew - pre');
%     vline([-std([shuff.rew_stim_means - shuff.pre_stim_means])*1.96, std([shuff.rew_stim_means - shuff.pre_stim_means])*1.96])
%     xlim([-20 20]);
%     

%% get the PETHs for the actual events. 
for ii = unique(ZoneIn)
    F_idx = find(ZoneIn == ii);
    figure(ii)
    [All_trial.outputS{ii}, All_trial.outputIT{ii}, All_trial.outputGau{ii}, All_trial.mean_S_gau{ii}, All_trial.pre_stim_means{ii}, All_trial.post_stim_means{ii},~, ~, All_trial.Z{ii}] = SpikePETH(cfg_peth, out.S, EnteringZoneTime(F_idx)/1000000);
    
    
    All_trial.mean_S_gau_z{ii} = (All_trial.mean_S_gau{ii} - shuff.mean_S_gau)./std_S_gau;
    
%     % re-define the 'pre' and 'post' means
%     pre_idx = nearest_idx3([-2.8 -1.8], All_trial.outputIT{ii});
%     act_idx = nearest_idx3([-0.5 0.5], All_trial.outputIT{ii});  
%     rew_idx = nearest_idx3([1.5 2.5], All_trial.outputIT{ii});  
    
    
    All_trial.pre_stim_diff{ii} = nanmean(All_trial.outputGau{ii}(pre_idx(1):pre_idx(2),:),1); 
    All_trial.act_stim_diff{ii} = nanmean(All_trial.outputGau{ii}(act_idx(1):act_idx(2),:),1); 
    All_trial.rew_stim_diff{ii} = nanmean(All_trial.outputGau{ii}(rew_idx(1):rew_idx(2),:),1);
    
    All_trial.pre_stim_means{ii} = nanmean(All_trial.mean_S_gau{ii}(pre_idx(1):pre_idx(2),:),1); 
    All_trial.act_stim_means{ii} = nanmean(All_trial.mean_S_gau{ii}(act_idx(1):act_idx(2),:),1); 
    All_trial.rew_stim_means{ii} = nanmean(All_trial.mean_S_gau{ii}(rew_idx(1):rew_idx(2),:),1);
    
    
    
    % zscore pre/post blocks realtive to shuffle

    All_trial.act_pre_z{ii} = ([All_trial.act_stim_means{ii} - All_trial.pre_stim_means{ii}]- z_mean_act_pre)./z_std_act_pre;
    All_trial.rew_pre_z{ii} = ([All_trial.rew_stim_means{ii} - All_trial.pre_stim_means{ii}]- z_mean_rew_pre)./z_std_rew_pre;
    All_trial.rew_act_z{ii} = ([All_trial.rew_stim_means{ii} - All_trial.act_stim_means{ii}]- z_mean_rew_act)./z_std_rew_act;

    
    
    % test for baseline vs action
    [All_trial.H_act_mod{ii}, All_trial.p_act_mod{ii},~] = ttest2(All_trial.act_stim_diff{ii} - All_trial.pre_stim_diff{ii},0,'Tail','both', 'alpha', 0.025); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
    [All_trial.H_pre_a{ii}, All_trial.p_pre_a{ii},~] = ttest2(All_trial.act_stim_diff{ii} - All_trial.pre_stim_diff{ii},0,'Tail','right', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
    [All_trial.H_act{ii}, All_trial.p_act{ii},~] = ttest2(All_trial.act_stim_diff{ii} - All_trial.pre_stim_diff{ii},0,'Tail','left', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

    if All_trial.H_act_mod{ii} == 1
        fprintf('<strong>%s</strong> Cell:  %s has significantly activity change at action (vs baseline) <strong>(p = %0.3f) at %s arm </strong>\n', mfilename, cell_to_process, All_trial.p_act_mod{ii}, Zone_names{ii});
    else
        fprintf('<strong>%s</strong> Cell:  %s is no significantly activity change at action (vs baseline) (p = %0.3f) at %s arm \n', mfilename, cell_to_process, All_trial.p_act_mod{ii}, Zone_names{ii});
    end
    
        % test for baseline vs action
    [All_trial.H_rew_mod{ii}, All_trial.p_rew_mod{ii},~] = ttest(All_trial.rew_stim_diff{ii} - All_trial.pre_stim_diff{ii},0,'Tail','both', 'alpha', 0.025); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
    [All_trial.H_pre_r{ii}, All_trial.p_pre_r{ii},~] = ttest(All_trial.rew_stim_diff{ii} - All_trial.pre_stim_diff{ii},0,'Tail','right', 'alpha', 0.05); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
    [All_trial.H_rew{ii}, All_trial.p_rew{ii},~] = ttest(All_trial.rew_stim_diff{ii} - All_trial.pre_stim_diff{ii},0,'Tail','left', 'alpha', 0.05); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

    if All_trial.H_rew_mod{ii} == 1
        fprintf('<strong>%s</strong> Cell:  %s has significantly activity change at reward (vs baseline) <strong>(p = %0.3f) at %s arm </strong>\n', mfilename, cell_to_process, All_trial.p_rew_mod{ii}, Zone_names{ii});
    else
        fprintf('<strong>%s</strong> Cell:  %s is no significantly activity change at reward (vs baseline) (p = %0.3f) at %s arm \n', mfilename, cell_to_process, All_trial.p_rew_mod{ii}, Zone_names{ii});
    end
    
    
    title(['PETH for Entry: ' Zone_names{ii} ' | Trials: ' num2str(length(F_idx))]);

    % test for sig modulation using t-test. 
%     [All_trial.H{ii}, All_trial.p{ii},~] = ttest(All_trial.post_stim_means{ii} - All_trial.pre_stim_means{ii},0,'Tail','both', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
%     [All_trial.H_pre{ii}, All_trial.p_pre{ii},~] = ttest(All_trial.act_stim_means{ii} - All_trial.pre_stim_means{ii},0,'Tail','left', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
%     [All_trial.H_act{ii}, All_trial.p_act{ii},~] = ttest(All_trial.act_stim_means{ii} - All_trial.pre_stim_means{ii},0,'Tail','right', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
%     [All_trial.H_post{ii}, All_trial.p_post{ii},~] = ttest(All_trial.rew_stim_means{ii} - All_trial.pre_stim_means{ii},0,'Tail','right', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
% 
%     

    
    saveas(gcf, ['Zone_plots' filesep out.S.label{1}(1:end-2) '_' Zone_names{ii} '_all_zcore.png'])
    saveas(gcf, ['Zone_plots' filesep out.S.label{1}(1:end-2)  '_' Zone_names{ii} '_all_zcore.fig'])
    
    print(gcf,['Zone_plots' filesep out.S.label{1}(1:end-2)  '_' Zone_names{ii} '_all_zcore'],'-depsc')
%     saveas_eps([out.S.label{1}(1:end-2)  '_' Feeder_names{ii} '_' cfg_peth.plot_type], cd)


end


% plot for all rewards this time using zscore. 
figure(1001)
cfg_all = cfg_peth;
cfg_all.evt_color_mat = Zone_cord;
cfg_all.plot_type = 'raw'; %'zscore';
cfg_all.markersize = 10; 
[All_trial.outputS{5}, All_trial.outputIT{5}, All_trial.outputGau{5}, All_trial.mean_S_gau{5}, All_trial.pre_stim_means{5}, All_trial.post_stim_means{5},~, ~, All_trial.Z{5}] = SpikePETH(cfg_all, out.S, EnteringZoneTime/1000000);
title(['PETH for Entry: ' Zone_names{5} ' | Trials: ' num2str(length(EnteringZoneTime))])


    All_trial.pre_stim_diff{5} = nanmean(All_trial.outputGau{5}(pre_idx(1):pre_idx(2),:),1); 
    All_trial.act_stim_diff{5} = nanmean(All_trial.outputGau{5}(act_idx(1):act_idx(2),:),1); 
    All_trial.rew_stim_diff{5} = nanmean(All_trial.outputGau{5}(rew_idx(1):rew_idx(2),:),1); 
    
    All_trial.pre_stim_means{5} = nanmean(All_trial.mean_S_gau{5}(pre_idx(1):pre_idx(2),:),1); 
    All_trial.act_stim_means{5} = nanmean(All_trial.mean_S_gau{5}(act_idx(1):act_idx(2),:),1); 
    All_trial.rew_stim_means{5} = nanmean(All_trial.mean_S_gau{5}(rew_idx(1):rew_idx(2),:),1);
    
    
    
    % zscore pre/post blocks realtive to shuffle

    All_trial.act_pre_z{5} = ([All_trial.act_stim_means{5} - All_trial.pre_stim_means{ii}]- z_mean_act_pre)./z_std_act_pre;
    All_trial.rew_pre_z{5} = ([All_trial.rew_stim_means{5} - All_trial.pre_stim_means{ii}]- z_mean_rew_pre)./z_std_rew_pre;
    All_trial.rew_act_z{5} = ([All_trial.rew_stim_means{5} - All_trial.act_stim_means{ii}]- z_mean_rew_act)./z_std_rew_act;


[All_trial.H_act_mod{5}, All_trial.p_act_mod{5},~] = ttest(All_trial.act_stim_diff{5} - All_trial.pre_stim_diff{5},0,'Tail','both', 'alpha', 0.025); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
[All_trial.H_pre_a{5}, All_trial.p_pre_a{5},~] = ttest(All_trial.act_stim_diff{5} - All_trial.pre_stim_diff{5},0,'Tail','right', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
[All_trial.H_act{5}, All_trial.p_act{5},~] = ttest(All_trial.act_stim_diff{5} - All_trial.pre_stim_diff{5},0,'Tail','left', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

if All_trial.H_act_mod{5} == 1
    fprintf('<strong>%s</strong> Cell:  %s has significantly activity change at action (vs baseline) <strong>(p = %0.3f) at %sarm </strong>\n', mfilename, cell_to_process, All_trial.p_act_mod{5}, Zone_names{5});
else
    fprintf('<strong>%s</strong> Cell:  %s is no significantly activity change at action (vs baseline) (p = %0.3f) at %s arm \n', mfilename, cell_to_process, All_trial.p_act_mod{5}, Zone_names{5});
end

[All_trial.H_rew_mod{5}, All_trial.p_rew_mod{5},~] = ttest(All_trial.rew_stim_diff{5} - All_trial.pre_stim_diff{5},0,'Tail','both', 'alpha', 0.025); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
[All_trial.H_pre_r{5}, All_trial.p_pre_r{5},~] = ttest(All_trial.rew_stim_diff{5} - All_trial.pre_stim_diff{5},0,'Tail','right', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
[All_trial.H_rew{5}, All_trial.p_rew{5},~] = ttest(All_trial.rew_stim_diff{5} - All_trial.pre_stim_diff{5},0,'Tail','left', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

if All_trial.H_rew_mod{5} == 1
    fprintf('<strong>%s</strong> Cell:  %s has significantly activity change at reward (vs baseline) <strong>(p = %0.3f) at %s arm </strong>\n', mfilename, cell_to_process, All_trial.p_rew_mod{5}, Zone_names{5});
elseif All_trial.H_rew_mod{5} == 1
        chil = get(gca, 'Children');
    chil(4).LineStyle = '-.'; 
else
    fprintf('<strong>%s</strong> Cell:  %s is no significantly activity change at reward (vs baseline) (p = %0.3f) at %s arm \n', mfilename, cell_to_process, All_trial.p_rew_mod{5}, Zone_names{5});
    chil = get(gca, 'Children');
    chil(4).LineStyle = '--'; 
end

% % test for sig modulation using t-test. 
% [All_trial.H{5}, All_trial.p{5},~] = ttest2(All_trial.post_stim_means{5} - All_trial.pre_stim_means{5},0,'Tail','both', 'alpha', 0.025); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
% [All_trial.H_pre{5}, All_trial.p_pre{5},~] = ttest2(All_trial.post_stim_means{5} - All_trial.pre_stim_means{5},0,'Tail','left', 'alpha', 0.05); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.
% [All_trial.H_post{5}, All_trial.p_post{5},~] = ttest2(All_trial.post_stim_means{5} - All_trial.pre_stim_means{5},0,'Tail','right', 'alpha', 0.05); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

% update plot based sig of all trial types.  '-' for sig, '--' for not. 

% if All_trial.H{5} ==1
%     fprintf('<strong>%s</strong> Cell:  %s is significantly modulated by reward  <strong> all reward (p = %0.3d)</strong>\n', mfilename, cell_to_process, All_trial.p{5});
% else
%         fprintf('<strong>%s</strong> Cell:  %s is not significantly modulated by reward  (p = %0.3f) at %s arm \n', mfilename, cell_to_process, All_trial.p{5});
% end

%% add the individual feeders to the overall PETH
figure(1001)
subplot(212)
hold on
for ii =unique(ZoneIn)
    if All_trial.H_rew_mod{ii} == 1
        plot(All_trial.outputIT{1}, All_trial.mean_S_gau{ii},'-',  'color', [c_ord(ii,:), .8])
    elseif All_trial.H_act_mod{ii} == 1
        plot(All_trial.outputIT{1}, All_trial.mean_S_gau{ii},'-.',  'color', [c_ord(ii,:), .8])
    else
        plot(All_trial.outputIT{1}, All_trial.mean_S_gau{ii},'--',  'color', [c_ord(ii,:), .8])
    end
end

y_lim = ([min([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}]) max([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}])]);
ylim(y_lim); 
% move the pre post means up
chil = get(gca, 'Children');
chil(5).Position = [chil(5).Position(1) max([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}])*.3 chil(5).Position(3)]; 
chil(6).Position = [chil(6).Position(1) max([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}])*.5 chil(6).Position(3)]; 
% adjust the vertical line at 0 which is a rectanlge. 
chil(7).Position = [chil(7).Position(1) y_lim(1) chil(7).Position(3) y_lim(2) - y_lim(1)]; 
leg = legend(['All', Zone_names(1:4)], 'Orientation', 'horizontal', 'Location', 'northwest','FontSize',14);
set(leg, 'box', 'off')

%% add velocity
velo_window = [cfg_peth.window(1)*floor(1/mode(diff(out.velo_smooth.tvec))), cfg_peth.window(2)*floor(1/mode(diff(out.velo_smooth.tvec)))]; 
all_velo = NaN(length(EnteringZoneTime), (abs(velo_window(1)) + abs(velo_window(2)) +1));
for ii = length(EnteringZoneTime):-1:1
    this_idx = nearest_idx3(EnteringZoneTime(ii)/1000000, out.velo_smooth.tvec);
    
        if this_idx < abs(velo_window(1)) || velo_window(2)+this_idx > length(out.velo_smooth.data)
            continue
        end

        
    all_velo(ii,:) = out.velo_smooth.data((velo_window(1)+this_idx):(velo_window(2)+this_idx));
    
end

velo_mean = nanmedian(all_velo, 1);
velo_tvec = cfg_peth.window(1) : 1/floor(1/mode(diff(out.velo_smooth.tvec))):cfg_peth.window(2);

yyaxis right
plot(velo_tvec, velo_mean);
ylabel('speed (cm/s)')
%%

SetFigure([], gcf);

    saveas(gcf, ['Zone_plots' filesep out.S.label{1}(1:end-2) '_all_zcore.png'])
    saveas(gcf, ['Zone_plots' filesep out.S.label{1}(1:end-2) '_all_zcore.fig'])
    print(gcf,['Zone_plots' filesep out.S.label{1}(1:end-2) '_all_zcore'],'-depsc')


    close all
    
    %% create a rate map for each zone and the all zone
    figure(918)
         title('Rate Map')
% 
%     cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = 5; % speed limit in cm/sec
%     iv_fast = TSDtoIV(cfg,out.velo_smooth); % only keep intervals with speed above thresh
%     
%     pos_r = restrict(out.pos,iv_fast);
%     S_r = restrict(out.S,iv_fast);
%     
%     spk_x = interp1(pos_r.tvec,pos_r.data(1,:),S_r.t{1},'linear');
%     spk_y = interp1(pos_r.tvec,pos_r.data(2,:),S_r.t{1},'linear');
%     
spk_x = interp1(out.pos.tvec,out.pos.data(1,:),out.S.t{1},'linear');
spk_y = interp1(out.pos.tvec,out.pos.data(2,:),out.S.t{1},'linear');

    % set up bins
    SET_xmin = 0; SET_ymin = 0; % set up bins
    SET_xmax = 180; SET_ymax = 180;
    SET_xBinSz = (SET_xmax - SET_xmin)/60; SET_yBinSz = (SET_xmax - SET_xmin)/60;
    
    
    x_edges = SET_xmin:SET_xBinSz:SET_xmax;
    y_edges = SET_ymin:SET_yBinSz:SET_ymax;
    
    % set up gaussian
    kernel = gausskernel([4 4],2); % 2d gaussian in bins
    
      % compute occupancy
    occ_hist = hist3(out.pos.data(1:2,:)', 'edges', {y_edges x_edges});
%     occ_hist = histcn(out.pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
    low_occ_idx = find(occ_hist <4 ); 
    occ_hist(low_occ_idx) = 0; 

%     occ_hist = conv2(occ_hist,kernel,'same');
    
    no_occ_idx = find(occ_hist == 0 ); % NaN out bins never visited
    occ_hist(no_occ_idx) = NaN;
%     
    occ_hist = occ_hist .* mode(diff(out.pos.tvec)); % convert samples to seconds using video frame rate (30 Hz)

    subplot(1,3,1)
    pcolor(occ_hist'); shading flat; axis off; c =colorbar;
    c.Label.String = 'nFrames'; c.Label.FontSize = 12; 
    c.Ticks = [min(c.Ticks) max(c.Ticks)];
    a = get(c, 'position');
    set(c,'Position',[a(1) a(2) 0.025 0.2]);
    set(gca, 'yDir', 'reverse')
    title('occupancy');
    

    % get the spike map
    spk_hist = hist3([spk_x, spk_y], 'edges', {y_edges x_edges});
%         spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    
%     spk_hist = conv2(spk_hist,kernel,'same');  
    spk_hist(no_occ_idx) = NaN;
    
    subplot(1,3,2)
    pcolor(spk_hist'); shading flat; axis off; c=colorbar;
    c.Label.String = 'nSpikes'; c.Label.FontSize = 12; 
    c.Ticks = [min(c.Ticks) max(c.Ticks)];
    a = get(c, 'position');
    set(c,'Position',[a(1) a(2) 0.025 0.2]);
    set(gca, 'yDir', 'reverse')

    title('spikes');
    
    % rate map
    tc = spk_hist./occ_hist;
    nan_idx = isnan(tc); 
    tc(nan_idx) = 0; 
    tc = conv2(tc,kernel,'same');
    tc(no_occ_idx) = NaN;

    
    subplot(1,3,3)
    pcolor(tc'); shading flat; axis off; c= colorbar; 
    c.Label.String = 'Hz'; c.Label.FontSize = 12; 
    c.Ticks = [min(c.Ticks) max(c.Ticks)];
    a = get(c, 'position');
    set(c,'Position',[a(1) a(2) 0.025 0.2]);
    set(gca, 'yDir', 'reverse')
    title('rate map');
        
        SetFigure([], gcf);
    set(gcf, 'position', [178,418,1578,468]); 

    saveas(gcf, ['Zone_plots' filesep out.S.label{1}(1:end-2) '_rate.png'])
    saveas(gcf, ['Zone_plots' filesep out.S.label{1}(1:end-2) '_rate.fig'])
    print(gcf,['Zone_plots' filesep out.S.label{1}(1:end-2) '_rate'],'-depsc')


    %% 
     figure(919)
     title('Entry Rate Map')
    
    % restrict to times around entries
    pos_r = restrict(out.pos,(EnteringZoneTime/1000000)-5, (EnteringZoneTime/1000000)+5);
    S_r = restrict(out.S,(EnteringZoneTime/1000000)-5, (EnteringZoneTime/1000000)+5);
    
    cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = 5; % speed limit in cm/sec
    iv_fast = TSDtoIV(cfg,out.velo_smooth); % only keep intervals with speed above thresh
    
    pos_r = restrict(pos_r,iv_fast);
    S_r = restrict(S_r,iv_fast);
    
    spk_x = interp1(pos_r.tvec,pos_r.data(1,:),S_r.t{1},'linear');
    spk_y = interp1(pos_r.tvec,pos_r.data(2,:),S_r.t{1},'linear');
%     
% spk_x = interp1(out.pos.tvec,out.pos.data(1,:),out.S.t{1},'linear');
% spk_y = interp1(out.pos.tvec,out.pos.data(2,:),out.S.t{1},'linear');

    % set up bins
    SET_xmin = 0; SET_ymin = 0; % set up bins
    SET_xmax = 180; SET_ymax = 180;
    SET_xBinSz = (SET_xmax - SET_xmin)/60; SET_yBinSz = (SET_xmax - SET_xmin)/60;
    
    
    x_edges = SET_xmin:SET_xBinSz:SET_xmax;
    y_edges = SET_ymin:SET_yBinSz:SET_ymax;
    
    % set up gaussian
    kernel = gausskernel([4 4],2); % 2d gaussian in bins
    
      % compute occupancy
    occ_hist = hist3(out.pos.data(1:2,:)', 'edges', {y_edges x_edges});
%     occ_hist = histcn(out.pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
    low_occ_idx = find(occ_hist <4 ); 
    occ_hist(low_occ_idx) = 0; 

%     occ_hist = conv2(occ_hist,kernel,'same');
    
    no_occ_idx = find(occ_hist == 0 ); % NaN out bins never visited
    occ_hist(no_occ_idx) = NaN;
%     
    occ_hist = occ_hist .* mode(diff(out.pos.tvec)); % convert samples to seconds using video frame rate (30 Hz)

    subplot(1,3,1)
    pcolor(occ_hist'); shading flat; axis off; c =colorbar;
    c.Label.String = 'nFrames'; c.Label.FontSize = 12; 
    c.Ticks = [min(c.Ticks) max(c.Ticks)];
    a = get(c, 'position');
    set(c,'Position',[a(1) a(2) 0.025 0.2]);
    set(gca, 'yDir', 'reverse')
    title('occupancy');
    

    % get the spike map
    spk_hist = hist3([spk_x, spk_y], 'edges', {y_edges x_edges});
%         spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    
%     spk_hist = conv2(spk_hist,kernel,'same');  
    spk_hist(no_occ_idx) = NaN;
    
    subplot(1,3,2)
    pcolor(spk_hist'); shading flat; axis off; c=colorbar;
    c.Label.String = 'nSpikes'; c.Label.FontSize = 12; 
    c.Ticks = [min(c.Ticks) max(c.Ticks)];
    a = get(c, 'position');
    set(c,'Position',[a(1) a(2) 0.025 0.2]);
    set(gca, 'yDir', 'reverse')

    title('spikes');
    
    % rate map
    tc = spk_hist./occ_hist;
    nan_idx = isnan(tc); 
    tc(nan_idx) = 0; 
    tc = conv2(tc,kernel,'same');
    tc(no_occ_idx) = NaN;

    
    subplot(1,3,3)
    pcolor(tc'); shading flat; axis off; c= colorbar; 
    c.Label.String = 'Hz'; c.Label.FontSize = 12; 
    c.Ticks = [min(c.Ticks) max(c.Ticks)];
    a = get(c, 'position');
    set(c,'Position',[a(1) a(2) 0.025 0.2]);
    set(gca, 'yDir', 'reverse')
    title('entry rate map');
        
        SetFigure([], gcf);
    set(gcf, 'position', [178,418,1578,468]); 

    saveas(gcf, ['Zone_plots' filesep out.S.label{1}(1:end-2) '_entry_rate.png'])
    saveas(gcf, ['Zone_plots' filesep out.S.label{1}(1:end-2) '_entry_rate.fig'])
    print(gcf,['Zone_plots' filesep out.S.label{1}(1:end-2) '_entry_rate'],'-depsc')
    
    
    %% add the individual feeders to the overall PETH [within trial zscore]
% figure(1001)
% subplot(212)
% hold on
% for ii =unique(FeedersFired)
%     if All_trial.H{ii} == 1
%         plot(All_trial.outputIT{1}, All_trial.Z{ii},'-',  'color', [c_ord(ii,:), .8])
%     else
%         plot(All_trial.outputIT{1}, All_trial.Z{ii},'--',  'color', [c_ord(ii,:), .8])
%     end
% end
% 
% ylim([min([All_trial.Z{1}; All_trial.Z{2}; All_trial.Z{3}; All_trial.Z{4}; All_trial.Z{5}]) max([All_trial.Z{1}; All_trial.Z{2}; All_trial.Z{3}; All_trial.Z{4}; All_trial.Z{5}])])
% 
% % move the pre post means up
% chil = get(gca, 'Children');
% chil(5).Position = [chil(5).Position(1) max([All_trial.Z{1}; All_trial.Z{2}; All_trial.Z{3}; All_trial.Z{4}; All_trial.Z{5}])*.6 chil(5).Position(3)]; 
% chil(6).Position = [chil(6).Position(1) max([All_trial.Z{1}; All_trial.Z{2}; All_trial.Z{3}; All_trial.Z{4}; All_trial.Z{5}])*.7 chil(6).Position(3)]; 
% 
% 
% leg = legend(['All', Feeder_names(1:4)], 'Orientation', 'horizontal', 'Location', 'northwest','FontSize',14);
% set(leg, 'box', 'off')
% 
% SetFigure([], gcf)
% 
%     saveas(gcf, [out.S.label{1}(1:end-2) '_all_' cfg_peth.plot_type '.png'])
%     saveas(gcf, [out.S.label{1}(1:end-2)  '_all_' cfg_peth.plot_type '.fig'])
%     print(gcf,[out.S.label{1}(1:end-2)  '_all_' cfg_peth.plot_type],'-depsc')
%% get the shuffles
% nShuff = 1000;
% cfg_shuff.plot ='off';
% cfg_shuff.window = cfg_peth.window;
% cfg_shuff.dt = cfg_peth.dt;
% 
% 
% tic
% fprintf('<strong>%s</strong> starting shuffle ...', mfilename);
% for iS = nShuff:-1:1
%     this_S = out.S;
%     for iT = 1:length(out.S.t)
%         this_S.t{iT} = (this_S.t{iT}(end)-this_S.t{iT}(1)).*rand(size(this_S.t{iT})) + this_S.t{iT}(1);
%     end
%     [~, ~, ~, Shuff.mean_gau(:,iS), Shuff.pre_means(:,iS), Shuff.post_means(:,iS)] = SpikePETH(cfg_shuff, this_S, FeederTimes./1000000);
%     
%     
% end
% toc
% 
% % get a shuffle of the spike train as a stand in for the baseline recording
% % period.  
% % tic
% % for iS = 100:-1:1
% %         shuff_FR(:,iS) = MS_get_gau_sdf([], pos.tvec, this_S.t{1}); 
% % end
% % toc
% % 
% % shuff_FR_mean  = std(shuff_FR, [], 1);
% 
% % get distribution of Reward index values for shuffle (post mean FR / pre mean FR)
% Shuff.Ridx = mean(Shuff.mean_gau(z_idx:end,:),1) - mean(Shuff.mean_gau(1:z_idx,:),1);


%% basic stats and distribution of cell vs shuffle.
% [Shuff.H, Shuff.p, Shuff.CI] = ttest(Shuff.Ridx,0, 'Alpha', 0.01); % test that shuff does not differ from 0.
% 
% if Shuff.H == 1
%     fprintf('<strong>%s</strong> Shuffle distribution was different from zero (p = %0.3d)\n', mfilename,  Shuff.p);
% end
% 
% figure(202)
% histogram(Shuff.Ridx, 'DisplayStyle', 'stairs')
% vline([prctile(Shuff.Ridx, 95), prctile(Shuff.Ridx, 99)], {'--r', '--r'})
% text(prctile(Shuff.Ridx, 95), max(ylim)*.95, '95%tile','color', 'r', 'fontweight', 'bold', 'fontsize', 8, 'fontname', 'helvetica')
% text(prctile(Shuff.Ridx, 99), max(ylim)*.95, '99%tile','color', 'r', 'fontweight', 'bold', 'fontsize', 8, 'fontname', 'helvetica')
% 
% vline([Shuff.CI(1), Shuff.CI(2)], {'--k', '--k'})
% text(Shuff.CI(1), max(ylim)*.95, 'L99CI', 'fontweight', 'bold', 'fontsize', 8, 'fontname', 'helvetica')
% text(Shuff.CI(2), max(ylim)*.95, 'U99CI', 'fontweight', 'bold', 'fontsize', 8, 'fontname', 'helvetica')
% 
% 
% h_vl = vline(mean(All_trial.post_stim_means{5} - All_trial.pre_stim_means{5}), 'k');
% h_vl.LineWidth = 3;
% text(mean(All_trial.post_stim_means{5} - All_trial.pre_stim_means{5}), max(ylim)*.95, 'This cell', 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'helvetica')
% 
% % add in the mean and std for the shuffle to the all trial plot
% % figure(1001)
% % subplot(212)
% % hold on
% % plot(outputIT{5}, mean(Shuff.mean_gau,2),'color', c_ord(2,:))
% % plot(outputIT{5}, mean(Shuff.mean_gau,2)+std(Shuff.mean_gau,[],2),'--', 'color', c_ord(2,:));
% % plot(outputIT{5}, mean(Shuff.mean_gau,2)-std(Shuff.mean_gau,[],2),'--', 'color', c_ord(2,:));

%% collect everything
out_temp = out; 
out = All_trial;
out.S = out_temp.S;
out.pos = out_temp.pos;
out.velo = out_temp.velo; 
% out.Shuff = Shuff; 
out.Zone_names = Zone_names;
out.Zone_cord = Zone_cord;
out.Zone_times = EnteringZoneTime;



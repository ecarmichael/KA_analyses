function [out] = KA_screener_feeder(cell_to_process)


c_ord = linspecer(5); % nice colours.

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
%% load data
cfg = [];
cfg.getTTnumbers = 0;
cfg.fc = {cell_to_process};
out.S = LoadSpikes(cfg);

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


cfg_pos.convFact = [560/142 480/142]; 
out.pos = LoadPos(cfg_pos);

% % circ smooth
% figure(919)
% plot(out.pos.tvec, unwrap(out.pos.data(3,:)*pi/180-pi), out.pos.tvec,smooth(unwrap(out.pos.data(3,:)*pi/180-pi), floor((1/mode(diff(out.pos.tvec))/4))));
% 
% figure(929)
% plot(out.pos.tvec,out.pos.data(3,:), out.pos.tvec,mod((180/pi)*smooth(unwrap(out.pos.data(3,:)*pi/180-pi), floor((1/mode(diff(out.pos.tvec))/4)))+180, 360))
% 
% 
% out.hd_smooth = rad2deg((unwrap(deg2rad(out.pos.data(3,:)))));%, floor((1/mode(diff(out.pos.tvec))/4))); % smooth hd over .25s


out.S = restrict(out.S, evt.t{1}(1), evt.t{2}(end)); % restrict the spikes recording periods.avoids odd thing where spike trains contains zeros.  MClust issue?
if length(out.S.t{1}) / (out.pos.tvec(end) - out.pos.tvec(1)) < .25
    fprintf('<strong>Minumum mean spike rate (0.5Hz) not met: %2.1fhz</strong>\n', length(out.S.t{1}) / (out.pos.tvec(end) - out.pos.tvec(1)));
    out =  length(out.S.t{1}) / (out.pos.tvec(end) - out.pos.tvec(1));  
    return
end



out.velo = getLinSpd([], out.pos); 
out.velo.data = interp1(out.velo.tvec,out.velo.data(1,:),out.pos.tvec,'linear');



% smooth speed over 0.5 seconds
out.velo_smooth = out.velo; 
out.velo_smooth.data = smooth(out.velo.data, round(1/mode(diff(out.pos.tvec)))*.5)'; % smooth with gaussian 

spk_x = interp1(out.pos.tvec,out.pos.data(1,:),out.S.t{1},'linear');
spk_y = interp1(out.pos.tvec,out.pos.data(2,:),out.S.t{1},'linear');


%  load PM control script vars
PM_dir = dir('PM*.mat');
load(PM_dir.name)

Feeder_names = {'North', 'West', 'South', 'East', 'All'};
Feeder_mag = [3 3 1 1];
Feeder_type = {'Banana', 'Grain', 'Banana', 'Grain'};
% 
% % check for trials beyond the recording time
beyond_idx = (FeederTimes/1000000) > out.pos.tvec(end);

FeedersFired = FeedersFired(~beyond_idx); 
FeederTimes = FeederTimes(~beyond_idx); 


for ii = length(FeederTimes):-1:1
    Feeder_cord(ii,:) = c_ord(FeedersFired(ii),:);
end

%% make a movie of the tracking for debugging
% figure(909)
% hold on
% plot(out.pos.data(1,:), out.pos.data(2,:), 'color', [.8 .8 .8])
% xlim([min(out.pos.data(1,:)) max(out.pos.data(1,:))]);
% ylim([min(out.pos.data(2,:)) max(out.pos.data(2,:))]);
% set(gca, 'YDir','reverse', 'XDir','reverse')
% Fs = mode(diff(out.pos.tvec)); 
% F1 = nearest_idx(FeederTimes(5)/1000000, out.pos.tvec); 
% 
% for ii = 1:2:length(out.pos.tvec)
%     if ii > F1 - 30 && ii < F1+30
%             h = plot(out.pos.data(1,ii), out.pos.data(2,ii), 'ob');
%     else
%     h = plot(out.pos.data(1,ii), out.pos.data(2,ii), 'or');
%     end
%     title(['time: ' num2str(out.pos.tvec(ii))])% - out.pos.tvec(1))]) 
%    drawnow 
%    pause(0.005)
%    delete(h)
% end

%% infer nosepoke from FeedersFired and tracking data. 
% S_pix_x = [85 98]; S_pix_y = 130;      
% W_pix_x = 25;     W_pix_y = [68, 82];
% N_pix_x = [70 85]; N_pix_y = 10;
% E_pix_x = 145;      E_pix_y = [50, 68]; 

S_pix_x = [70 100]; S_pix_y = 130;      
W_pix_x = 25;     W_pix_y = [55, 85];
N_pix_x = [70 100]; N_pix_y = 10;
E_pix_x = 145;      E_pix_y = [55, 85]; 

trial_vec= NaN(1,length(out.pos.tvec)); 

    
figure(109)
clf
hold on
plot(out.pos.data(1,:), out.pos.data(2,:), 'color', [.8 .8 .8])
xlim([min(out.pos.data(1,:)) max(out.pos.data(1,:))]);
ylim([min(out.pos.data(2,:)) max(out.pos.data(2,:))]);
set(gca, 'YDir','reverse', 'XDir','reverse')

h = hline(N_pix_y); h.Color = c_ord(1,:); h.LineStyle = '--'; 
h = vline(N_pix_x); h(1).Color = c_ord(1,:); h(1).LineStyle = '--'; h(2).Color = c_ord(1,:); h(2).LineStyle = '--'; 
h = vline(E_pix_x); h.Color = c_ord(2,:); h.LineStyle = '--'; 
h = hline(E_pix_y); h(1).Color = c_ord(2,:); h(1).LineStyle = '--'; h(2).Color = c_ord(2,:); h(2).LineStyle = '--'; 
h = hline(S_pix_y); h.Color = c_ord(3,:); h.LineStyle = '--'; 
h = vline(S_pix_x); h(1).Color = c_ord(3,:); h(1).LineStyle = '--'; h(2).Color = c_ord(3,:); h(2).LineStyle = '--'; 
h = vline(W_pix_x); h.Color = c_ord(4,:); h.LineStyle = '--'; 
h = hline(W_pix_y); h(1).Color = c_ord(4,:); h(1).LineStyle = '--'; h(2).Color = c_ord(4,:); h(2).LineStyle = '--'; 

plot(median([min(out.pos.data(1,:)) max(out.pos.data(1,:))]),  median([min(out.pos.data(2,:)) max(out.pos.data(2,:))]), 'x', 'MarkerSize', 44)

for iF = 1:length(FeedersFired)
    if iF == length(FeedersFired)
        this_trial.pos = restrict(out.pos, FeederTimes(iF)/1000000, out.pos.tvec(end));
        this_trial.velo_smooth = restrict(out.velo_smooth, FeederTimes(iF)/1000000,out.pos.tvec(end));
        trial_vec(nearest_idx(FeederTimes(iF)/1000000, out.pos.tvec): end)= FeedersFired(iF);% make an array of the trial type.
    else
        this_trial.pos = restrict(out.pos, FeederTimes(iF)/1000000, FeederTimes(iF+1)/1000000);
        this_trial.velo_smooth = restrict(out.velo_smooth, FeederTimes(iF)/1000000, FeederTimes(iF+1)/1000000);
        trial_vec(nearest_idx(FeederTimes(iF)/1000000, out.pos.tvec): nearest_idx(FeederTimes(iF+1)/1000000, out.pos.tvec)) = FeedersFired(iF); % make an array of the trial type.
    end
    
    %get the points within the corresponding pixels
    if FeedersFired(iF) == 1 % zone 1 'North'
        trial_idx = this_trial.pos.data(2,:) < N_pix_y; 
        trial_idx = trial_idx & (this_trial.pos.data(1,:) > N_pix_x(1))  & (this_trial.pos.data(1,:) < N_pix_x(2));
%         trial_idx = trial_idx & this_trial.pos.data(3,:) ;
        
    elseif FeedersFired(iF) == 2 % zone 2 'West'
        trial_idx = this_trial.pos.data(1,:) < W_pix_x;
        trial_idx = trial_idx & this_trial.pos.data(2,:) > W_pix_y(1)  & this_trial.pos.data(2,:) < W_pix_y(2);

    elseif FeedersFired(iF) == 3 % zone 3 'South'
        trial_idx = this_trial.pos.data(2,:) > S_pix_y; 
        trial_idx = trial_idx & this_trial.pos.data(1,:) > S_pix_x(1)  & this_trial.pos.data(1,:) < S_pix_x(2);

    else % zone 4 'East'
        trial_idx = this_trial.pos.data(1,:) > E_pix_x;
        trial_idx = trial_idx & this_trial.pos.data(2,:) > E_pix_y(1)  & this_trial.pos.data(2,:) < E_pix_y(2);
    end
    
    % get periods of no movement. 
    trial_idx = trial_idx & (this_trial.velo_smooth.data <= 2);
    % get the longest 

    [~, p_idx, p_w] = findpeaks(double(trial_idx), 'MinPeakWidth',floor(.2/mode(diff(this_trial.pos.tvec))));
    if isempty(p_idx)
       error_trial(iF)  = 1; enter_t(iF) = NaN; exit_t(iF) = NaN; 
%                plot(this_trial.pos.data(1,:), this_trial.pos.data(2,:),'--', 'color', [c_ord(FeedersFired(iF),:) .2], 'linewidth', 2)
%         pause
       continue
    else
        error_trial(iF) = 0; 
    end
%     [~, longest_pause] = max(p_idx); 

        enter_idx = p_idx(1);
        exit_idx = enter_idx +p_w(1); 

        % plot only the entry to the feeder
        
%         plot(this_trial.pos.data(1,:), this_trial.pos.data(2,:), 'color', [c_ord(FeedersFired(iF),:) .1], 'linewidth', 2)
        plot(this_trial.pos.data(1,1:exit_idx), this_trial.pos.data(2,1:exit_idx), 'color', [c_ord(FeedersFired(iF),:) .8], 'linewidth', 2)

% for ii = enter_idx:exit_idx
    plot(this_trial.pos.data(1,1), this_trial.pos.data(2,1),'d', 'color', c_ord(FeedersFired(iF),:)*.8, 'markersize', 14);
    plot(this_trial.pos.data(1,enter_idx), this_trial.pos.data(2,enter_idx),'o', 'color', c_ord(FeedersFired(iF),:)*.8, 'markersize', 14);
    plot(this_trial.pos.data(1,exit_idx), this_trial.pos.data(2,exit_idx),'x', 'color', c_ord(FeedersFired(iF),:)*.8, 'markersize', 14);

    %     plot(this_trial.pos.data(1,trial_idx(ii)), this_trial.pos.data(2,trial_idx(ii)),'o', 'color', c_ord(FeedersFired(iF),:), 'markersize', 10)
%   h = text(20, 120, num2str(iF), 'fontsize', 24);
%     pause(1)
%     delete(h); delete(hp)
% end
        drawnow
        pause(.05)
        disp(iF)
    enter_t(iF) = this_trial.pos.tvec(enter_idx); 
    exit_t(iF) = this_trial.pos.tvec(exit_idx); 
end

    fprintf('Trial hit rate: %0.2f%%\n', (1 -  sum(error_trial)/length(error_trial))*100)
legend({'path','Center',  'N trial','N enter', 'N poke', 'N exit', 'E trial',' E enter', 'E poke', 'E exit', 'S trial', 'S enter', 'S poke', 'S exit', 'W trial', 'W enter','W poke', 'W exit'})

%% test position and heading in time

figure(1010)
cla
hold on
subplot(5,1,1)
plot(this_trial.pos.tvec, this_trial.pos.data(1,:), 'b'); xlim([this_trial.pos.tvec(1) this_trial.pos.tvec(end)]); ylabel('x pos'); 
subplot(5,1,2)
plot(this_trial.pos.tvec, this_trial.pos.data(2,:), 'g'); xlim([this_trial.pos.tvec(1) this_trial.pos.tvec(end)]); ylabel('y pos'); 
subplot(5,1,3)
plot(this_trial.pos.tvec, this_trial.pos.data(3,:), 'm'); xlim([this_trial.pos.tvec(1) this_trial.pos.tvec(end)]); ylabel('HD'); 
subplot(5,1,4)
plot(this_trial.pos.tvec, this_trial.velo_smooth.data, 'r'); xlim([this_trial.pos.tvec(1) this_trial.pos.tvec(end)]); ylabel('Speed'); 
subplot(5,1,5)
plot(this_trial.pos.tvec, trial_idx, 'k'); ylim([-.5 1.5]); xlim([this_trial.pos.tvec(1) this_trial.pos.tvec(end)]); ylabel('Trial Idx'); 


legend({'X', 'Y', 'HD', 'Speed', 'Trial idx'})
figure(109)
legend({'path', 'N trial', 'N poke', 'E trial', 'E poke', 'S trial', 'S poke', 'W trial', 'W poke'})


%% make a replay of the tracking

% figure(99)
% plot(out.pos.data(1,:), out.pos.data(2,:), 'color', [.8 .8 .8])
%     hold on
% 
% 
% for ii = 301:1:floor(length(out.pos.tvec)/10)
%     
%     xlim([min(out.pos.data(1,:)) max(out.pos.data(1,:))]);
%     ylim([min(out.pos.data(2,:)) max(out.pos.data(2,:))]);
%     h1 = plot(out.pos.data(1,floor(ii - 300) : ii), out.pos.data(2,floor(ii - 300) : ii), 'color', [c_ord(trial_vec(ii),:) .2]);
%     h2 = plot(out.pos.data(1,ii), out.pos.data(2,ii), 'o', 'color', c_ord(trial_vec(ii),:), 'markersize', 10);
% 
%     drawnow 
% %     pause(.5)
%         delete(h1);
%     delete(h2); 
% end


% drawnow



%% Summary of spiking in time and space.
figure(1091)
subplot(2,3,1:2)
hold on
hold on
plot(out.pos.data(1,:), out.pos.data(2,:), 'color', [.7 .7 .7])

plot(spk_x,spk_y, '.r')
axis off

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
for ii = unique(FeedersFired)
    F_idx = find(FeedersFired == ii);
%     e_idx = nearest_idx3(FeederTimes(F_idx)/1000000, out.pos.tvec);
    F_idx(isnan(enter_t(F_idx))) = []; 
    e_idx = nearest_idx3(enter_t(F_idx), out.pos.tvec);
    e_diff = exit_t(F_idx) - enter_t(F_idx); 
    for jj = 1:length(e_idx)
        rectangle('position', [out.pos.tvec(e_idx(jj)),.5, e_diff(jj) , 1],'edgecolor', c_ord(ii,:), 'facecolor',[c_ord(ii,:) 0.2])
        text(out.pos.tvec(e_idx(jj)), 1.55, Feeder_names{ii}(1), 'color', c_ord(ii,:), 'HorizontalAlignment', 'center')
    end
end

PlotSpikeRaster2(cfg_rast, out.S);
xlim([out.pos.tvec(nearest_idx3(enter_t(F_idx(1)), out.pos.tvec))-5 out.pos.tvec(nearest_idx3(enter_t(F_idx(end)), out.pos.tvec))+5])


%% get the mean and std for the data to zscore the rate

cfg_z.binsize = 0.2;
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
cfg_peth.window = [-2.5 2.5];
%     cfg_peth.dt = 0.0025; % fine resulition for regular PETH.
cfg_peth.dt = cfg_z.binsize; % wider bins for 'Lap' plot later.
% cfg_peth.plot_type = 'zscore';
cfg_peth.plot_type = 'raw';
cfg_peth.z_mean = gau_mean; 
cfg_peth.z_std = gau_std;


for ii = unique(FeedersFired)
    F_idx = find(FeedersFired == ii);
    F_idx(isnan(enter_t(F_idx))) = []; 
    if isempty(F_idx)
        All_trial.H{ii} = NaN; 
        All_trial.mean_S_gau{ii} = NaN;
        continue
    end
    figure(ii)
    [All_trial.outputS{ii}, All_trial.outputIT{ii}, All_trial.outputGau{ii}, All_trial.mean_S_gau{ii}, All_trial.pre_stim_means{ii}, All_trial.post_stim_means{ii},~, ~, All_trial.Z{ii}] = SpikePETH(cfg_peth, out.S, enter_t(F_idx));
    
    title(['PETH for Feeder: ' Feeder_names{ii} ' | Reward mag: ' num2str(Feeder_mag(ii)) ' | ' Feeder_type{ii} ' | Trials: ' num2str(length(F_idx))]);
%     z_idx = find(All_trial.outputIT{ii} == 0);
%     All_trial.Z{ii} = (All_trial.mean_S_gau{ii} - mean(All_trial.mean_S_gau{ii}(1:z_idx)))./ std(All_trial.mean_S_gau{ii}(1:z_idx)); 
%     All_trial.Z{ii} = (All_trial.mean_S_gau{ii} - cfg_peth.z_mean)./ cfg_peth.z_std; 

    % test for sig modulation using t-test. 
    [All_trial.H{ii}, All_trial.p{ii},~] = ttest(All_trial.post_stim_means{ii} - All_trial.pre_stim_means{ii},0,'Tail','both', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

    [All_trial.H_pre{ii}, All_trial.p_pre{ii},~] = ttest(All_trial.post_stim_means{ii} - All_trial.pre_stim_means{ii},0,'Tail','left', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

    [All_trial.H_post{ii}, All_trial.p_post{ii},~] = ttest(All_trial.post_stim_means{ii} - All_trial.pre_stim_means{ii},0,'Tail','right', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

    
    if All_trial.H{ii} == 1
        fprintf('<strong>%s</strong> Cell:  %s has significantly increased activity following reward <strong>(p = %0.3f) at %s arm </strong>\n', mfilename, cell_to_process, All_trial.p{ii}, Feeder_names{ii});
    else
        fprintf('<strong>%s</strong> Cell:  %s is not significantly increased activity following reward (p = %0.3f) at %s arm \n', mfilename, cell_to_process, All_trial.p{ii}, Feeder_names{ii});
    end
    
    saveas(gcf, [out.S.label{1}(1:end-2) '_' Feeder_names{ii} '_all_zcore_enter.png'])
    saveas(gcf, [out.S.label{1}(1:end-2)  '_' Feeder_names{ii} '_all_zcore_enter.fig'])
    
    print(gcf,[out.S.label{1}(1:end-2)  '_' Feeder_names{ii} '_all_zcore_enter'],'-depsc')
%     saveas_eps([out.S.label{1}(1:end-2)  '_' Feeder_names{ii} '_' cfg_peth.plot_type], cd)


end


% plot for all rewards this time using zscore. 
figure(1001)
cfg_all = cfg_peth;
cfg_all.evt_color_mat = Feeder_cord;
cfg_all.plot_type = 'raw'; %'zscore';
cfg_all.markersize = 10; 
[All_trial.outputS{5}, All_trial.outputIT{5}, All_trial.outputGau{5}, All_trial.mean_S_gau{5}, All_trial.pre_stim_means{5}, All_trial.post_stim_means{5},~, ~, All_trial.Z{5}] = SpikePETH(cfg_all, out.S, enter_t);
title(['PETH for Feeder: ' Feeder_names{5} ' | Trials: ' num2str(length(enter_t))])


% % conver the PETH into a zscore using the pre-event period as the baseline.
% z_idx = find(All_trial.outputIT{5} == 0);
% All_trial.Z{5} = (All_trial.mean_S_gau{5} - cfg_peth.z_mean)./ cfg_peth.z_std; 
% 

% test for sig modulation using t-test. 
[All_trial.H{5}, All_trial.p{5},~] = ttest(All_trial.post_stim_means{5} - All_trial.pre_stim_means{5},0,'Tail','both', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

[All_trial.H_pre{5}, All_trial.p_pre{5},~] = ttest(All_trial.post_stim_means{5} - All_trial.pre_stim_means{5},0,'Tail','left', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

[All_trial.H_post{5}, All_trial.p_post{5},~] = ttest(All_trial.post_stim_means{5} - All_trial.pre_stim_means{5},0,'Tail','right', 'alpha', 0.01); % use a t-test to see if the distribution of post FR - pre FR across trials is greater than 0. Implying that the cell has significantly higher firing rate following reward delivery.

% update plot based sig of all trial types.  '-' for sig, '--' for not. 

if All_trial.H{5} ==1
    fprintf('<strong>%s</strong> Cell:  %s is significantly modulated by reward  <strong> all reward (p = %0.3d)</strong>\n', mfilename, cell_to_process, All_trial.p{5});
else
        fprintf('<strong>%s</strong> Cell:  %s is not significantly modulated by reward  (p = %0.3f) at %s arm \n', mfilename, cell_to_process, All_trial.p{5});
    chil = get(gca, 'Children');
    chil(4).LineStyle = '--'; 
end


%% add the individual feeders to the overall PETH
figure(1001)
subplot(212)
hold on
for ii =unique(FeedersFired)
    if ~isempty(All_trial.mean_S_gau{ii})
    if All_trial.H{ii} == 1
        plot(All_trial.outputIT{1}, All_trial.mean_S_gau{ii},'-',  'color', [c_ord(ii,:), .8])
    elseif isnan(All_trial.H{ii})
        continue
    else
        plot(All_trial.outputIT{1}, All_trial.mean_S_gau{ii},'--',  'color', [c_ord(ii,:), .8])
    end
    else
        All_trial.mean_S_gau{ii} = nan(size(All_trial.mean_S_gau{5})); 
    end
end

ylim([min([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}]) max([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}])])

% move the pre post means up
chil = get(gca, 'Children');
chil(end-3).Position = [chil(end-3).Position(1) max([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}])*.6 chil(end-3).Position(3)]; 
chil(end-2).Position = [chil(end-2).Position(1) max([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}])*.7 chil(end-2).Position(3)]; 

%% get the speed for each trial and superimpose it. 
velo_window = [cfg_peth.window(1)*floor(1/mode(diff(out.velo_smooth.tvec))), cfg_peth.window(2)*floor(1/mode(diff(out.velo_smooth.tvec)))]; 
all_velo = NaN(length(FeederTimes), (abs(velo_window(1)) + abs(velo_window(2)) +1));

enter_t_correct = enter_t(~error_trial); 
for ii = length(enter_t_correct):-1:1
    this_idx = nearest_idx(enter_t_correct(ii), out.velo_smooth.tvec);
    
        if this_idx < abs(velo_window(1)) || (this_idx + abs(velo_window(1)) > length(out.velo_smooth.data))
            continue
        end

    all_velo(ii,:) = out.velo_smooth.data((velo_window(1)+this_idx):(velo_window(2)+this_idx));
end

velo_mean = nanmedian(all_velo, 1);
velo_tvec = cfg_peth.window(1) : 1/floor(1/mode(diff(out.velo_smooth.tvec))):cfg_peth.window(2);

yyaxis right
% plot(velo_tvec, velo_mean, '--', );
hv = shadedErrorBar(velo_tvec,velo_mean,nanstd(all_velo) / sqrt(size(all_velo,1)));
hv.mainLine.Color =  [.5 .5 .5 .3];
hv.edge(1).Color =  [.5 .5 .5 .3];
hv.edge(2).Color =  [.5 .5 .5 .3];

hv.patch.FaceAlpha =  .2;

ylabel('speed (cm/s)')


%%
% y_lim = ylim
% rectangle('position', [0 min(mean_S_gau) 0.001  abs(max(mean_S_gau))*10], 'facecolor', [[4,172,218]./2550.5], 'edgecolor', [[4,172,218]./255 0.5])


subplot(212)
y_lim = ylim; 
rectangle('position', [0 y_lim(1) 0.001  y_lim(2) - y_lim(1)], 'facecolor', [[4,172,218]./255 0.5], 'edgecolor', [[4,172,218]./255 0.5])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [.5 .5 .5]; 

leg = legend(['All', Feeder_names(1:4) , 'speed'], 'Orientation', 'horizontal', 'Location', 'northwest','FontSize',14);
set(leg, 'box', 'off')

SetFigure([], gcf);

    saveas(gcf, [out.S.label{1}(1:end-2) '_all_zcore_enter.png'])
    saveas(gcf, [out.S.label{1}(1:end-2)  '_all_zcore_enter.fig'])
    print(gcf,[out.S.label{1}(1:end-2)  '_all_zcore_enter'],'-depsc')


    
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
out.error_trials = error_trial;
out.enter_t = enter_t; 
out.exit_t = exit_t;
% out.Shuff = Shuff; 
out.Feeder_names = Feeder_names;
out.Feeder_mag = Feeder_mag;
out.Feeder_type = Feeder_type;
out.Feeder_cord = Feeder_cord;



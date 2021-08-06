function [out] = KA_screener(cell_to_process)


c_ord = linspecer(5); % nice colours.
%% load data
cfg = [];
cfg.getTTnumbers = 0;
cfg.fc = {cell_to_process};
out.S = LoadSpikes(cfg);

if length(out.S.t{1}) < 1200
    out = []; 
    return
end
%load('TT1.ntt_01-wv.mat');

evt = LoadEvents([]);
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

out.S = restrict(out.S, evt.t{1}(1), evt.t{2}(end)); % restrict the spikes recording periods.avoids odd thing where spike trains contains zeros.  MClust issue?

out.velo = getLinSpd([], out.pos); 
out.velo.data = interp1(out.velo.tvec,out.velo.data(1,:),out.pos.tvec,'linear');

% smooth speed over 0.5 seconds
out.velo_smooth = out.velo; 
% out.velo_smooth.data = smooth(out.velo.data, round(1/mode(diff(out.pos.tvec)))*.5)'; % smooth with gaussian 

spk_x = interp1(out.pos.tvec,out.pos.data(1,:),out.S.t{1},'linear');
spk_y = interp1(out.pos.tvec,out.pos.data(2,:),out.S.t{1},'linear');


%  load PM control script vars
PM_dir = dir('PM*.mat');
load(PM_dir.name)

Feeder_names = {'North', 'West', 'South', 'East', 'All'};
Feeder_mag = [3 3 1 1];
Feeder_type = {'Banana', 'Grain', 'Banana', 'Grain'};

for ii = length(FeederTimes):-1:1
    Feeder_cord(ii,:) = c_ord(FeedersFired(ii),:);
end



%% Summary of spiking in time and space.
figure(101)
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
    e_idx = nearest_idx3(FeederTimes(F_idx)/1000000, out.pos.tvec);
    for jj = 1:length(e_idx)
        rectangle('position', [out.pos.tvec(e_idx(jj)) - 5, .5, 10, 1],'edgecolor', c_ord(ii,:), 'facecolor',[c_ord(ii,:) 0.2])
        text(out.pos.tvec(e_idx(jj)), 1.55, Feeder_names{ii}(1), 'color', c_ord(ii,:), 'HorizontalAlignment', 'center')
    end
end

PlotSpikeRaster2(cfg_rast, out.S);
xlim([out.pos.tvec(nearest_idx3(FeederTimes(1)/1000000, out.pos.tvec))-5 out.pos.tvec(nearest_idx3(FeederTimes(end)/1000000, out.pos.tvec))+5])


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
cfg_peth.window = [-5 5];
%     cfg_peth.dt = 0.0025; % fine resulition for regular PETH.
cfg_peth.dt = cfg_z.binsize; % wider bins for 'Lap' plot later.
% cfg_peth.plot_type = 'zscore';
cfg_peth.plot_type = 'raw';
cfg_peth.z_mean = gau_mean; 
cfg_peth.z_std = gau_std;


for ii = unique(FeedersFired)
    F_idx = find(FeedersFired == ii);
    figure(ii)
    [All_trial.outputS{ii}, All_trial.outputIT{ii}, All_trial.outputGau{ii}, All_trial.mean_S_gau{ii}, All_trial.pre_stim_means{ii}, All_trial.post_stim_means{ii},~, ~, All_trial.Z{ii}] = SpikePETH(cfg_peth, out.S, FeederTimes(F_idx)/1000000);
    
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
    
    saveas(gcf, [out.S.label{1}(1:end-2) '_' Feeder_names{ii} '_all_zcore.png'])
    saveas(gcf, [out.S.label{1}(1:end-2)  '_' Feeder_names{ii} '_all_zcore.fig'])
    
    print(gcf,[out.S.label{1}(1:end-2)  '_' Feeder_names{ii} '_all_zcore'],'-depsc')
%     saveas_eps([out.S.label{1}(1:end-2)  '_' Feeder_names{ii} '_' cfg_peth.plot_type], cd)


end


% plot for all rewards this time using zscore. 
figure(1001)
cfg_all = cfg_peth;
cfg_all.evt_color_mat = Feeder_cord;
cfg_all.plot_type = 'raw'; %'zscore';
cfg_all.markersize = 10; 
[All_trial.outputS{5}, All_trial.outputIT{5}, All_trial.outputGau{5}, All_trial.mean_S_gau{5}, All_trial.pre_stim_means{5}, All_trial.post_stim_means{5},~, ~, All_trial.Z{5}] = SpikePETH(cfg_all, out.S, FeederTimes/1000000);
title(['PETH for Feeder: ' Feeder_names{5} ' | Trials: ' num2str(length(FeederTimes))])


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
    if All_trial.H{ii} == 1
        plot(All_trial.outputIT{1}, All_trial.mean_S_gau{ii},'-',  'color', [c_ord(ii,:), .8])
    else
        plot(All_trial.outputIT{1}, All_trial.mean_S_gau{ii},'--',  'color', [c_ord(ii,:), .8])
    end
end

ylim([min([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}]) max([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}])])

% move the pre post means up
chil = get(gca, 'Children');
chil(5).Position = [chil(5).Position(1) max([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}])*.6 chil(5).Position(3)]; 
chil(6).Position = [chil(6).Position(1) max([All_trial.mean_S_gau{1}; All_trial.mean_S_gau{2}; All_trial.mean_S_gau{3}; All_trial.mean_S_gau{4}; All_trial.mean_S_gau{5}])*.7 chil(6).Position(3)]; 

%% get the speed for each trial and superimpose it. 
velo_window = [cfg_peth.window(1)*floor(1/mode(diff(out.velo_smooth.tvec))), cfg_peth.window(2)*floor(1/mode(diff(out.velo_smooth.tvec)))]; 
all_velo = NaN(length(FeederTimes), (abs(velo_window(1)) + abs(velo_window(2)) +1));
for ii = length(FeedersFired):-1:1
    this_idx = nearest_idx(FeederTimes(ii)/1000000, out.velo_smooth.tvec);
    
        if this_idx < abs(velo_window(1))
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

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [.5 .5 .5]; 

leg = legend(['All', Feeder_names(1:4) , 'speed'], 'Orientation', 'horizontal', 'Location', 'northwest','FontSize',14);
set(leg, 'box', 'off')

SetFigure([], gcf);

    saveas(gcf, [out.S.label{1}(1:end-2) '_all_zcore.png'])
    saveas(gcf, [out.S.label{1}(1:end-2)  '_all_zcore.fig'])
    print(gcf,[out.S.label{1}(1:end-2)  '_all_zcore'],'-depsc')


    
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
% out.Shuff = Shuff; 
out.Feeder_names = Feeder_names;
out.Feeder_mag = Feeder_mag;
out.Feeder_type = Feeder_type;
out.Feeder_cord = Feeder_cord;



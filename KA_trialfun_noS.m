function [data] = KA_trialfun_noS(plot_dir)

if nargin <1
    plot_flag = 0;
elseif nargin == 1
    plot_flag = 1;
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
task_rec_idx_s = task_rec_idx;
task_rec_idx_e = task_rec_idx;

if exist('PM-2021-04-29-09_41_50.mat', 'file') % odd session with a momentary stop in recording.
    task_rec_idx_s = 1;
    task_rec_idx_e = 2;
    rec_dur = evt.t{e_rec_idx}(task_rec_idx_e) - evt.t{s_rec_idx}(task_rec_idx_s);
end

% check for sessions that are too short. If less than 15mins skip.
if max(rec_dur)/60 < 10
    fprintf('<strong>Minumum recording duration (15mins) not met: %2.1fmins</strong>\n', max(rec_dur)/60);
    data = [];
    return
end


% relabel the events file for arm/feeders
evt.label{3} = 'NorthPellets';
evt.label{4} = 'EastPellets';
evt.label{5} = 'SdatahPellets';
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
data.pos = LoadPos(cfg_pos);
data.pos = restrict(data.pos, evt.t{s_rec_idx}(task_rec_idx_s), evt.t{e_rec_idx}(task_rec_idx_e)); % restrict position to time on track.

data.velo = getLinSpd([], data.pos);
data.velo.data = interp1(data.velo.tvec,data.velo.data(1,:),data.pos.tvec,'linear');

% smooth speed over 0.5 seconds
data.velo_smooth = data.velo;
data.velo_smooth.data = smooth(data.velo.data, round(1/mode(diff(data.pos.tvec)))*.5)'; % smooth with gaussian


%%% load the spikes %%%
% cfg = [];
% cfg.getTTnumbers = 0;
% t_files =  dir('*.t64');
% cfg.fc = {t_files.name};
% data.S = LoadSpikes(cfg);
% for iS = length(data.S.t):-1:1
%     this_fname = strrep(data.S.label{iS}, '.t64', ''); 
%     if exist([this_fname '-wv.mat'], 'file')
%         load([this_fname '-wv.mat'])
%         load([this_fname '-CluQual.mat'])
%         data.S.waves{iS}.mWV = mWV;
%         data.S.waves{iS}.xrange = xrange;
%         data.S.Q{iS} = CluSep; 
% 
%         clear mWV xrange CluSep
%     end
% end


% data.S = restrict(data.S, evt.t{s_rec_idx}(task_rec_idx_s), evt.t{e_rec_idx}(task_rec_idx_e)); % restrict the spikes recording periods.avoids odd thing where spike trains contains zeros.  MClust issue?

% keep_idx = ones(length(data.S.t));
% for iS = length(data.S.t):-1:1
%     if length(data.S.t{iS}) / (data.pos.tvec(end) - data.pos.tvec(1)) < min_fr
%         fprintf('<strong>Minumum mean spike rate (%0.2fHz) not met: %2.1fhz</strong>\n',min_fr, length(data.S.t{iS}) / (data.pos.tvec(end) - data.pos.tvec(1)));
% 
%         keep_idx(iS) = 0;
%     end
% 
% end


%%% load PM control script vars %%%
PM_dir = dir('PM*.mat');
load(PM_dir.name)

Zone_names = {'North', 'West', 'South', 'East', 'All'};
Feeder_mag = [3 3 1 1];
Feeder_type = {'Banana', 'Grain', 'Banana', 'Grain'};

keep_idx = (FeederTimes/1000000 > data.pos.tvec(1))  & (FeederTimes/1000000 < data.pos.tvec(end));
FeedersTimes = FeederTimes(keep_idx);  % remove any events after the end of the recording.
FeedersFired = FeedersFired(keep_idx);

% if rotation == 90
%     F_temp = NaN(size(FeedersFired));
%     F_temp(FeedersFired == 4) = 1;
%     F_temp(FeedersFired == 1) = 2;
%     F_temp(FeedersFired == 2) = 3;
%     F_temp(FeedersFired == 3) = 4;
%     FeedersFired = F_temp;
%     clear F_temp
% end
%
% zone_idx = nearest_idx3(EnteringZoneTime/1000000, data.pos.tvec);


%% infer nosepoke from FeedersFired and tracking data.
% S_pix_x = [85 98]; S_pix_y = 130;
% W_pix_x = 25;     W_pix_y = [68, 82];
% N_pix_x = [70 85]; N_pix_y = 10;
% E_pix_x = 145;      E_pix_y = [50, 68];

S_pix_x = [70 100]; S_pix_y = 130;
W_pix_x = 25;     W_pix_y = [55, 85];
N_pix_x = [70 100]; N_pix_y = 10;
E_pix_x = 145;      E_pix_y = [55, 85];

trial_vec= NaN(1,length(data.pos.tvec));

if plot_flag
    figure(109)
    clf
    hold on
    plot(data.pos.data(1,:), data.pos.data(2,:), 'color', [.8 .8 .8])
    xlim([min(data.pos.data(1,:)) max(data.pos.data(1,:))]);
    ylim([min(data.pos.data(2,:)) max(data.pos.data(2,:))]);
    set(gca, 'YDir','reverse', 'XDir','reverse')
    
    h = hline(N_pix_y); h.Color = c_ord(1,:); h.LineStyle = '--';
    h = vline(N_pix_x); h(1).Color = c_ord(1,:); h(1).LineStyle = '--'; h(2).Color = c_ord(1,:); h(2).LineStyle = '--';
    h = vline(E_pix_x); h.Color = c_ord(2,:); h.LineStyle = '--';
    h = hline(E_pix_y); h(1).Color = c_ord(2,:); h(1).LineStyle = '--'; h(2).Color = c_ord(2,:); h(2).LineStyle = '--';
    h = hline(S_pix_y); h.Color = c_ord(3,:); h.LineStyle = '--';
    h = vline(S_pix_x); h(1).Color = c_ord(3,:); h(1).LineStyle = '--'; h(2).Color = c_ord(3,:); h(2).LineStyle = '--';
    h = vline(W_pix_x); h.Color = c_ord(4,:); h.LineStyle = '--';
    h = hline(W_pix_y); h(1).Color = c_ord(4,:); h(1).LineStyle = '--'; h(2).Color = c_ord(4,:); h(2).LineStyle = '--';
    
    plot(median([min(data.pos.data(1,:)) max(data.pos.data(1,:))]),  median([min(data.pos.data(2,:)) max(data.pos.data(2,:))]), 'x', 'MarkerSize', 44)
end

for iF = 1:length(FeedersFired)
    if iF == length(FeedersFired)
        this_trial.pos = restrict(data.pos, FeederTimes(iF)/1000000, data.pos.tvec(end));
        this_trial.velo_smooth = restrict(data.velo_smooth, FeederTimes(iF)/1000000,data.pos.tvec(end));
        trial_vec(nearest_idx3(FeederTimes(iF)/1000000, data.pos.tvec): end)= FeedersFired(iF);% make an array of the trial type.
    else
        this_trial.pos = restrict(data.pos, FeederTimes(iF)/1000000, FeederTimes(iF+1)/1000000);
        this_trial.velo_smooth = restrict(data.velo_smooth, FeederTimes(iF)/1000000, FeederTimes(iF+1)/1000000);
        trial_vec(nearest_idx3(FeederTimes(iF)/1000000, data.pos.tvec): nearest_idx3(FeederTimes(iF+1)/1000000, data.pos.tvec)) = FeedersFired(iF); % make an array of the trial type.
    end
    
    %get the points within the corresponding pixels
    if FeedersFired(iF) == 1 % zone 1 'North'
        trial_idx = this_trial.pos.data(2,:) < N_pix_y;
        trial_idx = trial_idx & (this_trial.pos.data(1,:) > N_pix_x(1))  & (this_trial.pos.data(1,:) < N_pix_x(2));
        %         trial_idx = trial_idx & this_trial.pos.data(3,:) ;
        
    elseif FeedersFired(iF) == 2 % zone 2 'West'
        trial_idx = this_trial.pos.data(1,:) < W_pix_x;
        trial_idx = trial_idx & this_trial.pos.data(2,:) > W_pix_y(1)  & this_trial.pos.data(2,:) < W_pix_y(2);
        
    elseif FeedersFired(iF) == 3 % zone 3 'Sdatah'
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
    if plot_flag
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
    end
    %        pause(.05)
    %disp(iF)
    enter_t(iF) = this_trial.pos.tvec(enter_idx);
    exit_t(iF) = this_trial.pos.tvec(exit_idx);
end

fprintf('Trial hit rate: %0.2f%%\n', (1 -  sum(error_trial)/length(error_trial))*100)



for ii = length(FeedersFired):-1:1
    Zone_cord(ii,:) = c_ord(FeedersFired(ii),:);
end

zone_idx = nearest_idx3(enter_t(~logical(error_trial)), data.pos.tvec);
error_idx = nearest_idx3(enter_t(logical(error_trial)), data.pos.tvec);
rew_in = FeedersFired(~logical(error_trial));
rew_in_err = FeedersFired(logical(error_trial));

rew_t = enter_t(~logical(error_trial));
error_t = enter_t(logical(error_trial));

rew_cord = Zone_cord(~logical(error_trial),:);
rew_err_cord = Zone_cord(logical(error_trial),:);



app_idx = nearest_idx3(FeedersTimes(~logical(error_trial)), data.pos.tvec);
error_idx = nearest_idx3(enter_t(logical(error_trial)), data.pos.tvec);

app_in = FeedersFired(~logical(error_trial));
app_in_err = FeedersFired(logical(error_trial));

rew_t = enter_t(~logical(error_trial));
error_t = enter_t(logical(error_trial));

rew_cord = Zone_cord(~logical(error_trial),:);
rew_err_cord = Zone_cord(logical(error_trial),:);

if plot_flag
    if ~exist(plot_dir, 'dir')
        mkdir(plot_dir)
    end
    legend({'path','Center',  'N trial','N enter', 'N poke', 'N exit', 'E trial',' E enter', 'E poke', 'E exit', 'S trial', 'S enter', 'S poke', 'S exit', 'W trial', 'W enter','W poke', 'W exit'})

    title([data.pos.cfg.SessionID ' | ' num2str(sum(rew_in==1)) 'N' ' | ' num2str(sum(rew_in==2)) 'W' ' | ' num2str(sum(rew_in==3)) 'S' ' | ' num2str(sum(rew_in==4)) 'E'])
    saveas(gcf, [plot_dir filesep data.pos.cfg.SessionID '_Behav_trace.png'])
    close(109)
end
%% Summary of spiking in time and space.
% 
% if plot_flag
%     
%     spk_x = interp1(data.pos.tvec,data.pos.data(1,:),data.S.t{1},'linear');
%     spk_y = interp1(data.pos.tvec,data.pos.data(2,:),data.S.t{1},'linear');
%     
%     
%     
%     figure(101)
%     subplot(2,3,1)
%     hold on
%     plot(data.pos.data(1,:), data.pos.data(2,:), 'color', [.7 .7 .7])
%     
%     plot(spk_x,spk_y, '.r')
%     set(gca,'YDir','reverse')
%     axis off
%     
%     subplot(2,3,2)
%     hold on
%     plot(data.pos.data(1,:), data.pos.data(2,:), 'color', [.7 .7 .7])
%     gscatter(nan(1,4), nan(1,4),1:4, c_ord(1:4,:)); % space filler to get the colors in the legend to work.
%     for ii = 1:length(zone_idx) % hack to get the colours to work properly.
%         gscatter(data.pos.data(1,zone_idx(ii)), data.pos.data(2,zone_idx(ii)), rew_in(ii), rew_cord(ii,:), 'o', 3)
%         % gscatter(data.pos.data(1,zone_idx), data.pos.data(2,zone_idx), ZoneIn, Zone_cord, 'o+x.', 3)
%     end
%     axis off
%     set(gca,'YDir','reverse')
%     legend(['pos', Zone_names(1:4)])
%     % add in waveforms if they exist.  WIP.
%     % for ii = 1:size(mWV,2)
%     %     plot(500:(160/32):655,(mWV(:,ii)*.010)+400, 'color',c_ord(ii,:), 'linewidth', 2)
%     % end
%     
%     text(max(data.pos.data(1,:))*.1, max(data.pos.data(2,:))*.9, {['Firing rate = ' num2str(length(data.S.t{1})/(data.pos.tvec(end) - data.pos.tvec(1)),3) 'Hz'] ;['Mode ISI = ' num2str((mode(diff(data.S.t{1})))*10000, 5) 'ms'] ; ['Median ISI = ' num2str((median(diff(data.S.t{1})))*10000, 5) 'ms']})
%     
%     subplot(2,3,3)
%     histogram(diff(data.S.t{1})*10000, 0:1:500);
%     xlabel('ISI ms')
%     ylabel('Spike count')
%     
%     
%     subplot(2,3, 4:6)
%     cfg_rast = [];
%     
%     hold on
%     % for ii = unique(ZoneIn)
%     %     F_idx = find(ZoneIn == ii);
%     %     e_idx = nearest_idx3(EnteringZoneTime(F_idx)/1000000, data.pos.tvec);
%     %     for jj = 1:length(e_idx)
%     %         rectangle('position', [data.pos.tvec(e_idx(jj)) - 2, .5, 4, 1],'edgecolor', c_ord(ii,:), 'facecolor',[c_ord(ii,:) 0.2])
%     %         text(data.pos.tvec(e_idx(jj)), 1.55, Zone_names{ii}(1), 'color', c_ord(ii,:), 'HorizontalAlignment', 'center')
%     %     end
%     % end
%     
%     for ii = unique(FeedersFired)
%         F_idx = find(FeedersFired == ii);
%         %     e_idx = nearest_idx3(FeederTimes(F_idx)/1000000, data.pos.tvec);
%         F_idx(isnan(enter_t(F_idx))) = [];
%         e_idx = nearest_idx3(enter_t(F_idx), data.pos.tvec);
%         e_diff = exit_t(F_idx) - enter_t(F_idx);
%         for jj = 1:length(e_idx)
%             rectangle('position', [data.pos.tvec(e_idx(jj)),.5, e_diff(jj) , 1],'edgecolor', c_ord(ii,:), 'facecolor',[c_ord(ii,:) 0.2])
%             text(data.pos.tvec(e_idx(jj)), 1.55, Zone_names{ii}(1), 'color', c_ord(ii,:), 'HorizontalAlignment', 'center')
%         end
%     end
%     
%     PlotSpikeRaster2(cfg_rast, data.S);
%     xlim([data.pos.tvec(nearest_idx3(EnteringZoneTime(1)/1000000, data.pos.tvec))-2 data.pos.tvec(nearest_idx3(EnteringZoneTime(end)/1000000, data.pos.tvec))+2])
%     
% end

%% dataput

data.evt = evt;
data.rew.in = rew_in;
data.rew.t = rew_t;
data.rew.err_in = rew_in_err;
data.rew.err_t = error_t;
data.rew.cord = rew_cord;
data.rew.err_cord = rew_err_cord;

data.app.in = rew_in;
data.app.t = FeedersTimes(~logical(error_trial))/1000000;
data.app.err_in = rew_in_err;
data.app.err_t = FeedersTimes(logical(error_trial))/1000000;
data.app.cord = rew_cord;
data.app.err_cord = rew_err_cord;

data.PM.FeedersTimes = FeedersTimes/1000000; 
data.PM.FeedersFired = FeedersFired; 
data.PM.Feeder_mag = Feeder_mag;
data.PM.Feeder_type = Feeder_type; 

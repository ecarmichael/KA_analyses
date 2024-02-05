
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

%%
for iF = 1:length(FeedersFired)
    if iF == length(FeedersFired)
        this_trial.pos = restrict(out.pos, FeederTimes(iF)/1000000, out.pos.tvec(end));
        this_trial.velo_smooth = restrict(out.velo_smooth, FeederTimes(iF)/1000000,out.pos.tvec(end));
        trial_vec(nearest_idx3(FeederTimes(iF)/1000000, out.pos.tvec): end)= FeedersFired(iF);% make an array of the trial type.
    else
        this_trial.pos = restrict(out.pos, FeederTimes(iF)/1000000, FeederTimes(iF+1)/1000000);
        this_trial.velo_smooth = restrict(out.velo_smooth, FeederTimes(iF)/1000000, FeederTimes(iF+1)/1000000);
        trial_vec(nearest_idx3(FeederTimes(iF)/1000000, out.pos.tvec): nearest_idx3(FeederTimes(iF+1)/1000000, out.pos.tvec)) = FeedersFired(iF); % make an array of the trial type.
    end
    
    %f FeedersFired(iF) == 2
    %     hs = scatter(this_trial.pos.data(1,1), this_trial.pos.data(2,1), 55, 'k', 's')
    text(this_trial.pos.data(1,1)+(iF/3), this_trial.pos.data(2,1)+(iF/3),num2str(FeedersFired(iF)), 'color', c_ord(FeedersFired(iF),:), 'fontsize', 14)
    %     h = scatter(this_trial.pos.data(1,:), this_trial.pos.data(2,:), 5, winter(length(this_trial.pos.data(2,:))))
    h = plot(this_trial.pos.data(1,:), this_trial.pos.data(2,:),'color', [c_ord(FeedersFired(iF),:) .3], 'linewidth', 2);
    %     disp(iF)
    %     pause(.5)
    drawnow
    %     delete(h);
    %end
end




%%  check for sessions without aligned NLX and maze events


cd(data_dir)
% get all the sessions
this_dir = dir('*DONE');
sess_list = [];
for ii = 1:length(this_dir)
    if strcmp(this_dir(ii).name(1), '.') % check for hidden dirs
        continue
    else
        sess_list{ii} = this_dir(ii).name;
    end
end
sess_list =   sess_list(~cellfun('isempty',sess_list));

success = []; FR = [];
% loop over sessions in the data dir.
for iS =1:length(sess_list)
    
    if ismember(sess_list{iS}, omit_list)
        success(iS) = 99;
        continue
    end
    
    cd([data_dir filesep sess_list{iS}])
    
    % example
    
    cells_to_process = FindFiles('*.t64');
    
    % check if there are any .t files.  if not continue.
    if isempty(cells_to_process)
        success(iS) = 404;
        
        continue
    end
    
    if ~isempty(dir('*VT*.zip')) && isempty(dir('*.nvt'))
        unzip('VT1.zip')
    end
    
    
    %% load the NLX position and events
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
    
    cfg_pos.verbose = 0;
    cfg_pos.convFact = [560/142 480/142];
    out.pos = LoadPos(cfg_pos);
    out.pos = restrict(out.pos, evt.t{s_rec_idx}(task_rec_idx_s), evt.t{e_rec_idx}(task_rec_idx_e)); % restrict position to time on track.
    
    
    %% load the maze events
    
    PM_dir = dir('PM*.mat');
    load(PM_dir.name)
    
    Zone_names = {'North', 'West', 'South', 'East', 'All'};
    Feeder_mag = [3 3 1 1];
    Feeder_type = {'Banana', 'Grain', 'Banana', 'Grain'};
    
    
   keep_idx = (FeederTimes/1000000 > out.pos.tvec(1))  & (FeederTimes/1000000 < out.pos.tvec(end));
   
   if sum(~keep_idx) > 0
    fprintf('<strong>#%.0f %s %.0f/%.0f events fell outside of the NLX recroding (%0.2f%%)</strong>\n',iS, sess_list{iS}, sum(~keep_idx),length(keep_idx), sum(~keep_idx)/length(keep_idx))
   else    
   
    fprintf('#%.0f %s %.0f/%.0f events fell outside of the NLX recroding (%0.2f%%)\n',iS,  sess_list{iS}, sum(~keep_idx),length(keep_idx), sum(~keep_idx)/length(keep_idx))
   end
end
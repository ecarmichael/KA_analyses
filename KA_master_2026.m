%% KA screener master script.
% load data
usr_name = char(java.lang.System.getProperty('user.name'));

if ismac
    addpath(genpath(['/Users/' usr_name '/Documents/Github/vandermeerlab/code-matlab/shared']))
    addpath(genpath(['/Users/' usr_name '/Documents/Github/EC_State']));
    addpath(genpath(['/Users/' usr_name '/Documents/Github/KA_analyses']));
    % addpath(genpath(['/Users/' usr_name '/Documents/Github/NeuroLearning/M1_compiled_loaders_functions/releaseDec2015/binaries'])); % Apple silicon compiled Nlx Loaders

    data_dir = ['/Users/' usr_name '/Williams Lab Dropbox/Eric Carmichael/for_eric_only']; % where all the NLX data is.
    % inter_dir = '/Users/jericcarmichael/Dropbox/KA_Data/inter_';  % where to save the outputs.
    inter_dir = ['/Users/' usr_name '/Williams Lab Dropbox/Eric Carmichael/KA_Data/inter_reward_23'];
    plot_dir = ['/Users/' usr_name '/Williams Lab Dropbox/Eric Carmichael/KA_Data/Behav_plots'];

elseif ispc
    % load data
    addpath(genpath(['C:\Users\' usr_name '\Documents\GitHub\vandermeerlab\code-matlab\shared']))
    addpath(genpath(['C:\Users\' usr_name '\Documents\GitHub\EC_State']));
    addpath(genpath(['C:\Users\' usr_name '\Documents\GitHub\KA_analyses']));
    data_dir = ['C:\Users\' usr_name '\Williams Lab Dropbox\Eric Carmichael\for_eric_only']; % where all the NLX data is.
    inter_dir = ['C:\Users\' usr_name '\Williams Lab Dropbox\Eric Carmichael\KA_Data\inter_reward_23'];
    inter_dir_app = ['C:\Users\' usr_name '\Williams Lab Dropbox\Eric Carmichael\KA_Data\inter_reward_23_approach'];
    plot_dir = ['C:\Users\' usr_name '\Williams Lab Dropbox\Eric Carmichael\KA_Data\Behav_plots'];

else
    addpath(genpath('/home/ecar/Github/vandermeerlab/code-matlab/shared'))
    addpath(genpath('/home/ecar/Github/EC_State'));
    addpath(genpath('/home/ecar/Github/KA_analyses'));
    data_dir = '/lustre06/project/6064766/ecar/for_eric_only'; % where all the NLX data is.
    % inter_dir = '/Users/jericcarmichael/Dropbox/KA_Data/inter_';  % where to save the outputs.
    inter_dir = '/lustre06/project/6064766/ecar/KA_Data/inter_approach_2p5_new_sig';

end

cd(data_dir); % move to the data dir.

% make an intermediate directory if it doesn't exist.
if ~exist(inter_dir,'dir')
    mkdir(inter_dir)
end

% flagged sessions which contain some oddity like events outside of the
% recording

omit_list = {'C4_3_C3_2021-02-25_DONE',...
    'C5_2_O7_2021-04-30_DONE',... % Feeders start way before the recording.
    'C6_3_O1_2021-09-24_DONE',...
    'C6_3_O4_2021-09-27_DONE',...
    'C6_3_O5_2021-09-29_DONE',...
    'C6_4_O6_2021-09-27_DONE',... % HS likly fell out.
    'C6_4_E1_2021-10-07_DONE',... % 3 recordings.
    'C3_2_O7_2020-09-03_DONE',... % spike is all noise.
    'C6_3_E1_2021-10-13_DONE',... % short and split up. likly HS issue. No cells; 
    };

omit_cells =  {'C1_1 O6_2020-07-12_DONE_maze_data_TT1_01.t64',...
    'C3_3_O2_2020-09-03_DONE_maze_data_TT4_02.t64',...
    'C3_3_O4_2020-09-05_DONE_maze_data_TT4_01.t64',...
    'C3_3_O6_2020-09-07_DONE_maze_data_TT2_01.t64',...
    'C5_3_C3_2021-04-24_DONE_maze_data_TT2_01.t64',...
    'C5_3_O1_2021-04-25_DONE_maze_data_TT3_01.t64',...
    'C5_3_O2_2021-04-26_DONE_maze_data_TT3_01.t64',...
    'C5_3_O4_2021-04-28_DONE_maze_data_TT2_01.t64',...
    'C5_3_O4_2021-04-28_DONE_maze_data_TT4_01.t64',...
    'C6_3_O6_2021-09-30_DONE_maze_data_TT2_02.t64'};

min_fr = 0.1;

c_ord = linspecer(5);
c_map = [.05 .05 0.05 ; c_ord(2,:)]; %[6, 25, 34; 249, 160, 27]/255; % SUNS colour b/c KA.

[parent_path] = fileparts(inter_dir);
[~, parent_dir] = fileparts(parent_path);

f_id = {'North', 'West', 'South',  'East', 'Overall'};
flav_id = {'Grape x3','Orange x3', 'Grape x1', 'Orange x1', '   '};

rng(1111, 'twister') % set the rng seed

%% loop over sessions / cells
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
for iS =length(sess_list):-1:1

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

    data = KA_trialfun(min_fr, plot_dir);

    if isstruct(data) && ~isempty(data.S)
        save([inter_dir filesep sess_list{iS} '_maze_data.mat'], 'data', '-v7.3')
        success(iS) = 1;
    else
        success(iS) = -10;
    end



end

fprintf('<strong>%0.0f total sessions, %0.2f had good cells, %0.0f omitted, %0.0f no spike data, %0.0f too short</strong>\n', length(success), sum(success==1), sum(success==99), sum(success==404), sum(success==-10))
%% Process each cell within a session

rng(1111, 'twister')
cd(inter_dir)
sess_list = dir([inter_dir filesep 'C*.mat']);

load([inter_dir filesep 'ephys_dwell_times_all_sessions.mat'])

phase =[]; sub = []; cell_id = [];
app_out= [];
rew_out = [];
cell_id = {};
spd_mod = [];
spd_p = [];
spd_corr = [];
stats = [];
spd_data= [];
all_dwell= []; 
z_peta = []; 
pre_mean = []; post_mean = []; 
k = 0;
for iS = 1:length(sess_list)

    load([inter_dir filesep sess_list(iS).name])

    for iC = 1:length(data.S.t)

        if ismember([sess_list(iS).name(1:end-4) '_' data.S.label{iC}], omit_cells)
            continue
        end
        k = k+1;

        cell_id{k} = [sess_list(iS).name(1:end-4) '_' data.S.label{iC}];

        % isolate the cell of interest in the session (if there are
        this_S = KA_isolate_S(data.S, data.S.label{iC});

        % get some basic cell stats
        stats{k} = KA_Cell_stats(this_S, data.pos, 0);
        % pause
        % speed modulation ( add spd and acc MI later)
        data.spd_mod = KA_spd_mod([], this_S, data.velo_smooth);
        spd_mod(k) = data.spd_mod.spd_mod;
        spd_p(k) = data.spd_mod.p_val;
        spd_corr(k) = data.spd_mod.spd_corr;
        spd_z(k) = data.spd_mod.z_mod;

        spd_data{k}.FR = data.spd_mod.FR_velo_int;
        spd_data{k}.tvec =  data.velo_smooth.tvec;
        spd_data{k}.spd = data.velo_smooth.data;

        % get the dwell times and correct for NLX offset. 
        dwell_iv = KA_get_dwell(dwellExport.masterTable, cell_id{k}, 2.5); 
        dwell_iv.tstart = dwell_iv.tstart - spd_data{k}.tvec(1); 
        dwell_iv.tend= dwell_iv.tend - spd_data{k}.tvec(1); 


        % summary for plotting
        cfg_peth = [];
        cfg_peth.window = [-3 3];
        cfg_peth. plot_type = 'raw';
        cfg_peth.dt = 0.05;
        cfg_peth.gauss_sd = .1;
        cfg_peth.plot = 'on'; 
         % get the rate
        S_vec = MS_spike2rate(this_S, spd_data{k}.tvec, cfg_peth.dt, cfg_peth.gauss_sd); 
        app_t{iS}= data.app;



        % PETA method
        S_vec.tvec = S_vec.tvec - S_vec.tvec(1); 
        % zscore the data relative to the pseudobaseline. Add 2.5 to keep
        % the reward consumption separate.
        % dwell_S_vec = restrict(S_vec, dwell_iv.tstart+2.5, dwell_iv.tend);
        % all_dwell.mean(k) =  mean(dwell_S_vec.data, 'omitmissing');
        % all_dwell.std(k) =  std(dwell_S_vec.data, 'omitmissing');

        % S_vec.data = (S_vec.data - mean(dwell_S_vec.data, 'omitmissing'))./std(dwell_S_vec.data, 'omitmissing'); 
        for jj = 1:4
            this_idx = data.rew.in == jj;
            if sum(data.rew.t(this_idx)) <1
                this_peta = nan(size(all_peta,1));
            else
                this_peta = KA_PETA(S_vec, data.rew.t(this_idx) - spd_data{k}.tvec(1), cfg_peth.window);
                for iShuff = 500:-1:1
                    this_peta_s(:,iShuff) = mean(KA_PETA(S_vec, MS_randn_range(1,100,cfg_peth.window(1), S_vec.tvec(end)-cfg_peth.window(end)), cfg_peth.window),1, 'omitmissing');
                end
            end
            all_peta(:,jj, k) = mean(this_peta,1, 'omitmissing');
            z_peta(:,jj,k) = mean((this_peta' - mean(this_peta_s, 2, 'omitmissing'))./std(this_peta_s, [], 2, 'omitmissing'),2); 
        end
        [this_peta, tvec_peta] = KA_PETA(S_vec, data.rew.t - spd_data{k}.tvec(1), cfg_peth.window);
        all_peta(:,5, k) = mean(this_peta,1, "omitmissing");
        z_peta(:,5,k) = mean((this_peta' - mean(this_peta_s, 2, 'omitmissing'))./std(this_peta_s, [], 2, 'omitmissing'),2); 

        % test for increased/decreased activity realtive to shuffle
        r_idx = nearest_idx3([0 2], tvec_peta); % get the mean response z score in the reward window
        rew_out.z_mean_fr(k,:) = mean(z_peta(r_idx(1):r_idx(2),:,k)); 
        rew_mean_fr = mean(all_peta(r_idx(1):r_idx(2),:, k),"omitmissing");
        rew_s_mu = mean(this_peta_s(r_idx(1):r_idx(2),:), 1,"omitmissing");
        rew_out.z_mean_fr_z(k,:) = mean((rew_mean_fr' - mean(rew_s_mu, 2, 'omitmissing'))...
            ./std(rew_s_mu, [], 2, 'omitmissing'),2, 'omitmissing');
        rew_out.z_mean_h(k,:) =  abs(rew_out.z_mean_fr(k,:)) > 1.96; 
        rew_out.z_h(k,:) =  sum(abs(z_peta(r_idx(1):r_idx(2),:,k)) > 1.96) > 0; 


        % same but for the approach
        for jj = 1:4
            this_idx = data.app.in == jj;
            if sum(data.app.t(this_idx)) <1
                this_peta = nan(size(all_peta,1));
            else
                this_peta = KA_PETA(S_vec, data.app.t(this_idx) - spd_data{k}.tvec(1), cfg_peth.window);
                for iShuff = 500:-1:1
                    this_peta_s(:,iShuff) = mean(KA_PETA(S_vec, MS_randn_range(1,100,cfg_peth.window(1), S_vec.tvec(end)-cfg_peth.window(end)), cfg_peth.window),1, 'omitmissing');
                end
            end
            all_peta_app(:,jj, k) = mean(this_peta,1, 'omitmissing');
            z_peta_app(:,jj,k) = mean((this_peta' - mean(this_peta_s, 2, 'omitmissing'))./std(this_peta_s, [], 2, 'omitmissing'),2);
        end
        [this_peta, tvec_peta] = KA_PETA(S_vec, data.app.t - spd_data{k}.tvec(1), cfg_peth.window);
        all_peta_app(:,5, k) = mean(this_peta,1, "omitmissing");
        z_peta_app(:,5,k) = mean((this_peta' - mean(this_peta_s, 2, 'omitmissing'))./std(this_peta_s, [], 2, 'omitmissing'),2);

       % test for increased/decreased activity realtive to shuffle
        r_idx = nearest_idx3([-1 1], tvec_peta); % get the mean response z score in the reward window
        app_out.z_mean_fr(k,:) = mean(z_peta_app(r_idx(1):r_idx(2),:,k)); 
        app_mean_fr = mean(all_peta_app(r_idx(1):r_idx(2),:, k),"omitmissing"); 
        app_s_mu = mean(this_peta_s(r_idx(1):r_idx(2),:), 1,"omitmissing"); 
        app_out.z_mean_fr_z(k,:) = mean((app_mean_fr' - mean(app_s_mu, 2, 'omitmissing'))...
            ./std(app_s_mu, [], 2, 'omitmissing'),2, 'omitmissing');
        app_out.z_mean_h(k,:) =  abs(app_out.z_mean_fr(k,:)) > 1.96; 
        app_out.z_h(k,:) =  sum(abs(z_peta_app(r_idx(1):r_idx(2),:,k)) > 1.96) > 0; 




        phase{k} = sess_list(iS).name(6:7);
        sub{k} = str2double(sess_list(iS).name([2 4]));

    end
end

for ii = size(rew_out.h,2):-1:1
    rew_mod(ii) = sum(sum(rew_out.h(:,ii),2)>0);

    app_mod(ii) = sum(sum(app_out.h(:,ii),2)>0);
end

rew_no_mod = sum(sum(rew_out.h(:,1:5),2) == 0);
app_no_mod = sum(sum(app_out.h(:,1:5),2) == 0);

% save([parent_path filesep 'Spd_data.mat'], 'spd_data')

%% sort cells in the output matrix by phase

[phase_sort, phase_sort_idx] = sort(phase); % lucky that they are already in a sortable format. 

phase_sort_idx = flipud(fliplr(phase_sort_idx));
phase_sort = flipud(fliplr(phase_sort));

c_idx = find(contains(phase_sort, 'C'));
oe_idx = find(contains(phase_sort, {'O1', 'O2', 'O3'}));
ol_idx = find(contains(phase_sort, { 'O4','O5', 'O6', 'O7'}));
r_idx = find(contains(phase_sort, 'R'));

% sort the speed vectors
spd_z = spd_z(phase_sort_idx); 
spd_corr = spd_corr(phase_sort_idx); 
spd_p = spd_p(phase_sort_idx); 
spd_mod = spd_mod(phase_sort_idx); 

spd_data = spd_data(phase_sort_idx);
% sort the stats
stats = stats(phase_sort_idx); 

% sort the cell IDS
cell_id = cell_id(phase_sort_idx); 

% phase colours
c = linspecer(5);

c_ord = c; % one for each session type C, O, E, R.
c_ord(2,:) =  c(4,:);
c_ord(3,:) = c(2,:);
c_ord(4,:) = c(5,:);

blues = parula(16); reds = jet(64); oranges = autumn(16);
sess_cord = [flipud(blues(3:2:7,:));(reds(end-7:end-1,:)); c_ord(3,:); flipud(oranges(end-9:end-6,:))];

c_orange = [255 150 0]/255;
c_l_orange = [255 213 153]/255;
c_purple = [98 66 158]/255;
c_l_purple = [198 183 225]/255;

c_blue = [ 0.3639    0.5755    0.7484]; 
c_red = [0.9153    0.2816    0.2878]; 

%% generate histograms of the maximal and miniaml firing rates that exceed +/- 1.96sd (using PETH and PETA)
kk = 5; 
ft_size = 12; % fontsize
% test it
figure(10); clf

this_all_peta = squeeze(z_peta_app(:,5,phase_sort_idx));
[~, pos_max_idx] = max(this_all_peta,[], 1, 'omitmissing');
[~, pos_min_idx] = min(this_all_peta,[], 1, 'omitmissing');
% sig_p_idx = sum(this_all_peta> 1.96,1) > 0; 
% sig_n_idx = sum(this_all_peta < -1.96,1) > 0; 
sig_p_idx = app_out.z_mean_fr(phase_sort_idx,5) > 1.96;
sig_n_idx = app_out.z_mean_fr(phase_sort_idx,5) < -1.96; 

subplot(6,3,[1 4 7]); cla
hold on

 imagesc(tvec_peta, 1:size(all_peta,3), this_all_peta')
% plot(tvec_peta(pos_max_idx), 1:size(all_peth,3), 'x')
plot(tvec_peta(pos_max_idx(sig_p_idx)), find(sig_p_idx), 'x', 'MarkerSize',8, 'LineWidth',2)
plot(tvec_peta(pos_min_idx(sig_n_idx)), find(sig_n_idx), 'xr', 'MarkerSize',8, 'LineWidth',2)
xlim([cfg_peth.window])
clim([-5 5])
ylim([.5 size(all_peta,3)+.5])
xline(0, '--k', 'LineWidth',2); 
title('Shuffle Normalized PETA: approach')
h = line([cfg_peth.window(end)+.1 cfg_peth.window(end)+.1], [c_idx(1) c_idx(end)], 'Color',c_ord(1,:), 'linewidth', 4); h.Clipping = 'off';
h = line([cfg_peth.window(end)+.1 cfg_peth.window(end)+.1], [oe_idx(1), oe_idx(end)], 'Color',c_ord(2,:), 'linewidth', 4); h.Clipping = 'off';
h = line([cfg_peth.window(end)+.1 cfg_peth.window(end)+.1], [ol_idx(1),ol_idx(end)], 'Color',c_ord(3,:), 'linewidth', 4); h.Clipping = 'off';
h = line([cfg_peth.window(end)+.1 cfg_peth.window(end)+.1], [r_idx(1) r_idx(end)], 'Color',c_ord(4,:), 'linewidth', 4); h.Clipping = 'off';
ylabel('Cells')
set(gca, "XTick", [],'ytick', [1 size(all_peta,3)])


subplot(6,3,10); cla
hold on
s_mat = this_all_peta(:,sig_p_idx | sig_n_idx); 
imagesc(tvec_peta, 1:size(s_mat,2), s_mat')
clim([-5 5])
ylim([.5 size(s_mat,2)+.5])
xlim([cfg_peth.window])
xline(0, '--k', 'LineWidth',2); 
% title('Sig Approach ONLY')
ylabel('Responsive cells')
set(gca, "XTick", [])
set(gca, 'ytick', [1 size(s_mat,2)])

subplot(6,3,13); cla
hold on
ylim([-1.5 1.5])
plot(tvec_peta, mean(this_all_peta(:,sig_p_idx)./max(this_all_peta(:,sig_p_idx)), 2,"omitmissing"), '-', 'color', c_blue, 'LineWidth',2)
plot(tvec_peta, mean(this_all_peta(:,sig_n_idx)./max(this_all_peta(:,sig_n_idx)), 2,"omitmissing"),'-', 'color', c_red, 'LineWidth',2)
xlim([cfg_peth.window])
xline(0, '--k', 'LineWidth',2); 
legend({'Pos mean', 'Neg mean', ''}, 'Box','off')
ylabel('norm FR')
set(gca, 'ytick', [-1.5 0 1.5]);
set(gca, "XTick", [])


subplot(6,3,16); cla; hold on;

histogram(tvec_peta(pos_max_idx((sig_p_idx))), cfg_peth.window(1):.25:cfg_peth.window(2), 'FaceColor',c_blue,'facealpha',.7,'edgecolor','none', 'Normalization', 'percentage');
histogram(tvec_peta(pos_min_idx((sig_n_idx))), cfg_peth.window(1):.25:cfg_peth.window(2), 'FaceColor',c_red,'facealpha',.3,'edgecolor','none', 'Normalization','percentage');
[hy, hx] = hist(tvec_peta(pos_max_idx),  cfg_peth.window(1):.1:cfg_peth.window(2)); 
plot(hx, (hy./sum(hy))*100, 'Color', [c_blue .4], 'LineWidth',2); 
[hy, hx] = hist(tvec_peta(pos_min_idx),  cfg_peth.window(1):.1:cfg_peth.window(2)); 
plot(hx, (hy./sum(hy))*100, 'Color', [c_red, .4], 'LineWidth',2); 
ylim([0 25])
xlim([cfg_peth.window])
xlabel('Time from reward (s)')
xline(0, '--k', 'LineWidth',2); 
legend({ ['Sig Pos ' num2str((sum(sig_p_idx)/length(sig_p_idx))*100,3) '%'], ['Sig Neg ' num2str((sum(sig_n_idx)/length(sig_n_idx))*100,3) '%'],'Pop maxima', 'Pop minimum', ''}, 'Box', 'off')
ylabel('% cells')


% Reward PETA
this_all_peta = squeeze(z_peta(:,5,phase_sort_idx));

[~, pos_max_idx] = max(this_all_peta,[], 1, 'omitmissing');
[~, pos_min_idx] = min(this_all_peta,[], 1, 'omitmissing');
sig_p_idx = rew_out.z_mean_fr(phase_sort_idx,5) > 1.96;
sig_n_idx = rew_out.z_mean_fr(phase_sort_idx,5) < -1.96; 


subplot(6,3,[2 5 8]); cla
hold on
 imagesc(tvec_peta, 1:size(all_peta,3), this_all_peta')
% plot(tvec_peta(pos_max_idx), 1:size(all_peth,3), 'x')
plot(tvec_peta(pos_max_idx(sig_p_idx)), find(sig_p_idx), 'x', 'MarkerSize',8, 'LineWidth',2)
plot(tvec_peta(pos_min_idx(sig_n_idx)), find(sig_n_idx), 'xr', 'MarkerSize',8, 'LineWidth',2)
xlim([cfg_peth.window])
clim([-5 5])
ylim([.5 size(all_peta,3)+.5])
xline(0, '--k', 'LineWidth',2); 
title('Shuffle Normalized PETA: Reward')
h = line([cfg_peth.window(end)+.1 cfg_peth.window(end)+.1], [c_idx(1) c_idx(end)], 'Color',c_ord(1,:), 'linewidth', 4); h.Clipping = 'off';
h = line([cfg_peth.window(end)+.1 cfg_peth.window(end)+.1], [oe_idx(1), oe_idx(end)], 'Color',c_ord(2,:), 'linewidth', 4); h.Clipping = 'off';
h = line([cfg_peth.window(end)+.1 cfg_peth.window(end)+.1], [ol_idx(1),ol_idx(end)], 'Color',c_ord(3,:), 'linewidth', 4); h.Clipping = 'off';
h = line([cfg_peth.window(end)+.1 cfg_peth.window(end)+.1], [r_idx(1) r_idx(end)], 'Color',c_ord(4,:), 'linewidth', 4); h.Clipping = 'off';
set(gca, 'ytick', [1 size(all_peta,3)])
cb = colorbar; 
cb.Position(1) = cb.Position(1) + .05;
% cb.Label.String = 'Mean speed (zscore)';
cb.Label.FontSize = ft_size;
cb.FontSize = ft_size;
cb.Ticks = [-5 -2.5 0 2.5 5]; cb.Label.String = 'mean activity (zscore)';
set(gca, "XTick", [])


subplot(6,3,11); cla
hold on
s_mat = this_all_peta(:,sig_p_idx | sig_n_idx); 
imagesc(tvec_peta, 1:size(s_mat,2), s_mat')
clim([-5 5])
ylim([.5 size(s_mat,2)+.5])
xlim([cfg_peth.window])
xline(0, '--k', 'LineWidth',2); 
set(gca, "XTick", [])
set(gca, 'ytick', [1 size(s_mat,2)])

% title('Sig Reward ONLY')

subplot(6,3,14); cla
hold on
plot(tvec_peta, mean(this_all_peta(:,sig_p_idx)./max(this_all_peta(:,sig_p_idx)), 2,"omitmissing"), '-', 'color', c_blue, 'LineWidth',2)
plot(tvec_peta, mean(this_all_peta(:,sig_n_idx)./max(this_all_peta(:,sig_n_idx)), 2,"omitmissing"),'-', 'color', c_red, 'LineWidth',2)
xlim([cfg_peth.window])
ylim([-1.5 1.5])
xline(0, '--k', 'LineWidth',2); 
legend({'Pos mean', 'Neg mean', ''}, 'Box','off')
set(gca, 'ytick', [-1.5 0 1.5]);
set(gca, "XTick", [])

subplot(6,3,17); cla; hold on;
histogram(tvec_peta(pos_max_idx((sig_p_idx))), cfg_peth.window(1):.25:cfg_peth.window(2), 'FaceColor',c_blue,'facealpha',.7,'edgecolor','none', 'Normalization', 'percentage');
histogram(tvec_peta(pos_min_idx((sig_n_idx))), cfg_peth.window(1):.25:cfg_peth.window(2), 'FaceColor',c_red,'facealpha',.3,'edgecolor','none', 'Normalization','percentage');
[hy, hx] = hist(tvec_peta(pos_max_idx),  cfg_peth.window(1):.1:cfg_peth.window(2)); 
plot(hx, (hy./sum(hy))*100, 'Color', [c_blue .4], 'LineWidth',2); 
[hy, hx] = hist(tvec_peta(pos_min_idx),  cfg_peth.window(1):.1:cfg_peth.window(2)); 
plot(hx, (hy./sum(hy))*100, 'Color', [c_red, .4], 'LineWidth',2); 

ylim([0 25])
xlim([cfg_peth.window])
xlabel('Time from reward (s)')
xline(0, '--k', 'LineWidth',2); 
legend({ ['Sig Pos ' num2str((sum(sig_p_idx)/length(sig_p_idx))*100,3) '%'], ['Sig Neg ' num2str((sum(sig_n_idx)/length(sig_n_idx))*100,3) '%'],'Pop maxima', 'Pop minimum', ''}, 'Box', 'off')

% colormap(plasma(512))
cfg_fig = []; 
cfg_fig.ft_size = ft_size; 
SetFigure(cfg_fig, gcf)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,[parent_path filesep 'PETA.pdf'],'-dpdf','-r300')
%% classify cells based on waveform properties

rng(1111, 'twister') % set the rng seed again in case something else got run. 


fr = []; bur_idx = []; s_w = []; pt_r = []; rfint = []; wave_dur = []; wave_forms = []; isi = []; pt_sym = []; 
for ii = length(stats):-1:1
    fr(ii) = stats{ii}.firing_rate;
    bur_idx(ii) = stats{ii}.burst_idx;
    s_w(ii) = stats{ii}.spike_width*1000;
    s_r(ii) = stats{ii}.slopes_ratio; 
    isi(ii) = stats{ii}.ISI_cv;
    pt_sym(ii) = abs(stats{ii}.pt_sym_r(1,2));
    pt_r(ii) = (stats{ii}.pt_ratio);
    rfint(ii) = stats{ii}.rise_fall_inter;
    wave_dur(ii) = stats{ii}.wave_dur*1000;
    wave_forms(:,ii) = stats{ii}.wave(stats{ii}.best_wave_idx,:);
end
figure(808); clf
subplot(2, 2, 1)

% [g_idx, n_val] = MS_kmean_scatter([(fr)', bur_idx',wave_dur'], 2, [1,2,3], 50);
g_idx = clusterdata([(fr)', bur_idx',rfint'], maxclust=8, Linkage="centroid");
scatter3(fr', bur_idx',rfint', 100, g_idx, "filled")
xlabel('Firing rate (Hz)');
ylabel('burst index')
zlabel('spike width')

% try clustering from a different method
% subplot(2, 2, 4)



colors = linspecer(3);  %or any other way of creating the colormap
% set(gca, 'XScale', 'log')
% fprintf('Clustering returned %0.0f groups based on firing rate, burst index, and spike width <strong>G1: %0.0f%% G2: %0.0f%% G3: %0.0f%%</strong>\n',...
%     length(unique(g_idx)), (n_val(1)/length(g_idx))*100, (n_val(2)/length(g_idx))*100, (n_val(3)/length(g_idx))*100)

FS_idx = g_idx == 2;
PC_idx = ~FS_idx;

f_name  = fieldnames(rew_out);
rew_FS = []; rew_PC = [];
app_FS = []; app_PC = [];

for ii = 1:length(f_name)

    rew_FS.(f_name{ii}) = rew_out.(f_name{ii})(FS_idx);
    rew_PC.(f_name{ii}) = rew_out.(f_name{ii})(PC_idx);

    app_FS.(f_name{ii}) = app_out.(f_name{ii})(FS_idx);
    app_PC.(f_name{ii}) = app_out.(f_name{ii})(PC_idx);

end

% plot the mean waveforms for the types of cells
x_range = reshape(this_S.waves{1}.xrange, 1, 32*4);

for ii = 1:length(unique(g_idx))
    subplot(2,2,ii+1)
    cla
    hold on
    for jj = 1
        this_wave = wave_forms(:, g_idx == ii)./max(wave_forms(:, g_idx == ii)); 
        x_range = (this_S.waves{1}.xrange(:,jj) - this_S.waves{1}.xrange(1,jj))./32;
        % plot(x_range , nanmean(this_wave,2))
        % errorbar(this_S.waves{1}.xrange(:,jj), nanmean(wave_forms(jj+1,:, g_idx == ii),3) +std(wave_forms(jj+1,:, g_idx == ii),[],3))
        % plot(x_range, this_wave)
        s = shadedErrorBar(x_range, mean(this_wave,2),std(this_wave,[],2, 'omitmissing') ); 
        s.mainLine.LineWidth = 3; 
        s.mainLine.Color = colors(ii,:);
        s.patch.FaceColor = colors(ii,:);
        s.patch.FaceAlpha = 0.2;
        s.edge(1).Color = colors(ii,:);
        s.edge(2).Color = colors(ii,:);

    end
    title(['Group ' num2str(ii) ' | n=' num2str(sum(g_idx == ii))  ' | FR:' num2str(round(mean(fr(g_idx == ii)),2)) '\pm' num2str(round(std(fr(g_idx == ii)),2))], 'Interpreter','tex')
        % xlim([x_range(1) x_range(end)])
ylabel('voltage (normalized)'); xlabel('time (ms)')
set(gca, 'xtick', [0 .5 1])
end
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf,[parent_path filesep 'wave_feature.pdf'],'-dpdf','-r300')

%% collect the responses over the phases
p_types = unique(phase);
Phase = [];rew_mat = []; app_mat = [];

this_rew_h = abs(rew_out.z_mean_fr) > 1.96; 
this_app_h = abs(app_out.z_mean_fr) > 1.96; 

rew_z_sort = rew_out.z_mean_fr(phase_sort_idx,:);
app_z_sort = app_out.z_mean_fr(phase_sort_idx,:); 


for iP = length(p_types):-1:1
    this_idx = strcmpi(p_types{iP}, phase);
    for ii = size(rew_out.h,2):-1:1
        Phase.rew_mod{iP}(:,ii) = sum(sum(this_rew_h(this_idx,ii),2)>0)/sum(this_idx);
        Phase.app_mod{iP}(:,ii) = sum(sum(this_app_h(this_idx,ii),2)>0)/sum(this_idx);

    end
    rew_mat(iP,:) = Phase.rew_mod{iP};
    app_mat(iP,:) = Phase.app_mod{iP};
end

% get the 'any modulation'
for iP = length(p_types):-1:1
    this_idx = strcmpi(p_types{iP}, phase);

    % any mod
        Phase.rew_mod{iP}(:,6) = sum(sum(this_rew_h(this_idx,1:4),2)>0)/sum(this_idx);
        Phase.app_mod{iP}(:,6) = sum(sum(this_app_h(this_idx,1:4),2)>0)/sum(this_idx);

    rew_mat(iP,6) = Phase.rew_mod{iP}(:,6);
    app_mat(iP,6) = Phase.app_mod{iP}(:,6);
    
    % Large reward only
    Phase.rew_mod{iP}(:,7) = sum(sum(this_rew_h(this_idx,1:2),2)>0)/sum(this_idx);
    Phase.app_mod{iP}(:,7) = sum(sum(this_app_h(this_idx,1:2),2)>0)/sum(this_idx);

    rew_mat(iP,7) = Phase.rew_mod{iP}(:,7);
    app_mat(iP,7) = Phase.app_mod{iP}(:,7);

        % small reward only
    Phase.rew_mod{iP}(:,8) = sum(sum(this_rew_h(this_idx,3:4),2)>0)/sum(this_idx);
    Phase.app_mod{iP}(:,8) = sum(sum(this_app_h(this_idx,3:4),2)>0)/sum(this_idx);

    rew_mat(iP,8) = Phase.rew_mod{iP}(:,8);
    app_mat(iP,8) = Phase.app_mod{iP}(:,8);
end


figure(1099)
clf

subplot(3,2,[1 3])
imagesc(1:8, 1:length(p_types), rew_mat*100);
set(gca, 'xtick', [1:8], 'XTickLabel', {'North', 'West', 'South', 'East', 'overall', 'any', 'large', 'small'});
set(gca,'YTick', 1:length(p_types),  'YTickLabel', p_types)
caxis([0 100]); c = colorbar('Location', 'eastoutside');
c.Ticks = [0 25 50 75 100]; c.Label.String = '% modulated cells';
title('Reward')

subplot(3,2,[2 4]); cla;
imagesc(1:8, 1:length(p_types), app_mat*100);
set(gca, 'xtick', [1:8], 'XTickLabel', {'North', 'West', 'South', 'East', 'overall', 'any', 'large', 'small'});
set(gca,'YTick', 1:length(p_types),  'YTickLabel', p_types)
caxis([0 100]); c = colorbar('Location', 'eastoutside');
c.Ticks = [0 25 50 75 100]; c.Label.String = '% modulated cells';
title('Approach')


% overalls
subplot(3,2,5); cla; hold on
MS_bar_w_err(rew_mat(1:3,6)'*100,rew_mat(4:7,6)'*100, [c_ord(1,:); c_ord(2,:)], 1,'ttest2', 1:2);
MS_bar_w_err(rew_mat(8:10,6)'*100,rew_mat(11:13,6)'*100, [c_ord(3,:); c_ord(4,:)],1, 'ttest2', 3:4);
set(gca, 'xtick', 1:4, 'XTickLabel', {'C', 'O_e', 'O_l', 'R'});
xlim([0.5 4.5])
ylabel('% modulated cells');
title('Reward any')

subplot(3,2,6); cla; hold on
MS_bar_w_err(app_mat(1:3,6)'*100,app_mat(4:7,6)'*100, [c_ord(1,:); c_ord(2,:)], 1,'ttest2', 1:2);
MS_bar_w_err(app_mat(8:10,6)'*100,app_mat(11:13,6)'*100, [c_ord(3,:); c_ord(4,:)],1, 'ttest2', 3:4);
set(gca, 'xtick', 1:4, 'XTickLabel', {'C', 'O_e', 'O_l', 'R'});
xlim([0.5 4.5])
ylabel('% modulated cells');
title('Approach any')

figure(10909)
clf
subplot(2,2,1); cla; hold on
[~, ~, ~,t_p, t_stats] = MS_bar_w_err(rew_mat(:,7)'*100,rew_mat(:,8)'*100, [c_ord(1,:); c_ord(2,:)], 1,'ttest', 1:2);
set(gca, 'xtick', 1:2, 'XTickLabel', {'large', 'small'});
xlim([0.5 2.5])
ylabel('% modulated cells');
title('Reward any')

subplot(2,2,2); cla; hold on
[~, ~, ~,a_t_p, a_t_stats] = MS_bar_w_err(app_mat(:,7)'*100,app_mat(:,8)'*100, [c_ord(1,:); c_ord(2,:)], 1,'ttest', 1:2);
set(gca, 'xtick', 1:2, 'XTickLabel', {'large', 'small'});
xlim([0.5 2.5])
ylabel('% modulated cells');
title('Approach any')


%%%%%%%% receate original plots but with the zscores
figure(2001)
clf
subplot(2,3,1); cla; hold on
[~, ~, ~,t_p, t_stats] = MS_bar_w_err(rew_z_sort(c_idx,5),rew_z_sort(oe_idx,5), [c_ord(1,:); c_ord(2,:)], 1,'ttest2', 1:2); 
[~, ~, ~,t_p, t_stats] = MS_bar_w_err(rew_z_sort(ol_idx,5),rew_z_sort(r_idx,5), [c_ord(3,:); c_ord(4,:)], 1,'ttest2', 3:4);
[~, ~, ~,t_p, t_stats] = MS_bar_w_err(rew_z_sort(oe_idx,5),rew_z_sort(ol_idx,5), [c_ord(2,:); c_ord(3,:)], 0,'ttest2', 2:3);
[~, ~, ~,t_p, t_stats] = MS_bar_w_err(rew_z_sort(c_idx,5),rew_z_sort(ol_idx,5), [c_ord(1,:); c_ord(3,:)], 0,'ttest2', [1 3]);
[~, ~, ~,t_p, t_stats] = MS_bar_w_err(rew_z_sort(c_idx,5),rew_z_sort(r_idx,5), [c_ord(1,:); c_ord(4,:)], 0,'ttest2', [1 4]);
[~, ~, ~,t_p, t_stats] = MS_bar_w_err(rew_z_sort(oe_idx,5),rew_z_sort(r_idx,5), [c_ord(2,:); c_ord(3,:)], 0,'ttest2', [2 4]);

set(gca, 'xtick', [1:4], 'XTickLabel', {'Critia', 'Early Overtrain', 'Late overtrain', 'Post-devaluation'}, 'XTickLabelRotation', 45);
xlim([0.5 4.5])
ylim([-8 14]);
ylabel('Reward response (zscore)');


subplot(2,3,2); cla; hold on
[~, ~, ~,t_p, t_stats] = MS_bar_w_err([rew_z_sort(c_idx,1); rew_z_sort(c_idx,3)] ,[rew_z_sort(c_idx,2); rew_z_sort(c_idx,4)], [c_purple; c_orange], 1,'ttest2', 1:2); 
[~, ~, ~,t_p, t_stats] = MS_bar_w_err([rew_z_sort(oe_idx,1); rew_z_sort(oe_idx,3)] ,[rew_z_sort(oe_idx,2); rew_z_sort(oe_idx,4)], [c_purple; c_orange], 1,'ttest2', 4:5); 
[~, ~, ~,t_p, t_stats] = MS_bar_w_err([rew_z_sort(ol_idx,1); rew_z_sort(ol_idx,3)] ,[rew_z_sort(ol_idx,2); rew_z_sort(ol_idx,4)], [c_purple; c_orange], 1,'ttest2', 7:8); 
[~, ~, ~,t_p, t_stats] = MS_bar_w_err([rew_z_sort(r_idx,1); rew_z_sort(r_idx,3)] ,[rew_z_sort(r_idx,2); rew_z_sort(r_idx,4)], [c_purple; c_orange], 1,'ttest2', 10:11); 

set(gca, 'xtick', [1.5 4.5 7.5 10.5], 'XTickLabel', {'Critia', 'Early Overtrain', 'Late overtrain', 'Post-devaluation'}, 'XTickLabelRotation', 45);
xlim([0.5 11.5])
ylim([-10 25]);
ylabel('Reward response (zscore)');
title('Reward any')

subplot(2,3,3); cla; hold on
[~, ~, ~,t_p, t_stats] = MS_bar_w_err([rew_z_sort(c_idx,3); rew_z_sort(c_idx,4)] ,[rew_z_sort(c_idx,1); rew_z_sort(c_idx,2)], [.8 .8 .8; .5 .5 .5], 1,'ttest2', 1:2); 
[~, ~, ~,t_p, t_stats] = MS_bar_w_err([rew_z_sort(oe_idx,3); rew_z_sort(oe_idx,4)] ,[rew_z_sort(oe_idx,1); rew_z_sort(oe_idx,2)], [.8 .8 .8; .5 .5 .5], 1,'ttest2', 4:5); 
[~, ~, ~,t_p, t_stats] = MS_bar_w_err([rew_z_sort(ol_idx,3); rew_z_sort(ol_idx,4)] ,[rew_z_sort(ol_idx,1); rew_z_sort(ol_idx,2)], [.8 .8 .8; .5 .5 .5], 1,'ttest2', 7:8); 
[~, ~, ~,t_p, t_stats] = MS_bar_w_err([rew_z_sort(r_idx,3); rew_z_sort(r_idx,4)] ,[rew_z_sort(r_idx,1); rew_z_sort(r_idx,2)], [.8 .8 .8; .5 .5 .5], 1,'ttest2', 10:11); 

set(gca, 'xtick', [1.5 4.5 7.5 10.5], 'XTickLabel', {'Critia', 'Early Overtrain', 'Late overtrain', 'Post-devaluation'}, 'XTickLabelRotation', 45);
xlim([0.5 11.5])
ylim([-10 25]);
ylabel('Reward response (zscore)');
title('Reward any')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,[parent_path filesep 'zscore_rew_resp.pdf'],'-dpdf','-r300')


%%  collect data for stats
cell_type = cell(size(g_idx)); 
cell_type(g_idx == 2) = {'FS'}; 
cell_type(g_idx == 1) = {'RS'}; 

for ii = 1:length(sub)
    sub{ii} = num2str(sub{ii}); 
end
tbl_out = table(sub', phase_sort', cell_type, rew_z_sort(:,1), rew_z_sort(:,2), rew_z_sort(:,3), rew_z_sort(:,4), rew_z_sort(:,5),...
    app_z_sort(:,1), app_z_sort(:,2), app_z_sort(:,3), app_z_sort(:,4), app_z_sort(:,5),...
     'VariableNames',{'Subject', 'Phase', 'Cell_type','G_large_z', 'O_large_z', 'G_small_z', 'O_small_z', 'R_overall_z',...
    'app_G_large_z', 'app_O_large_z', 'app_G_small_z', 'app_O_small_z', 'app_R_overall_z',...
    }); 

writetable(tbl_out, [parent_path filesep 'PETA_tbl.csv'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%  Speed coding   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% try some linear decoding using the Lopes et al. 2025 PlosOne method with a 3rd order polynomial fit.
z_err = NaN(length(spd_data),1);
R2 = z_err;
S_R2 = cell(size(R2));

parfor kk = 1:length(spd_data)
    disp(kk)
    [z_err(kk), R2(kk), S_R2{kk}] = KA_lin_decode([], spd_data{kk}.FR, spd_data{kk}.spd);
end
% example cell
ex_idx = 140;
[~, ~, ~, plt_mat, plt_s_mat] = KA_lin_decode([], spd_data{ex_idx}.FR, spd_data{ex_idx}.spd);

nan_idx = isnan(spd_data{ex_idx}.FR) | isnan(spd_data{ex_idx}.spd); % remove nans to match plt_mat

% get the mean shuffle R2 per cell.
S_R2_means = [];
for ii = length(S_R2):-1:1
    S_R2_means(ii,1) = mean(S_R2{ii});
end

% append the speed decoding to the tbl_out 

tbl_out.spd_z = spd_z'; 
tbl_out.spd_mod = spd_mod'; 
tbl_out.spd_corr = spd_corr'; 
tbl_out.spd_p = spd_p'; 

% speed decoding
tbl_out.spd_lin_z_err = z_err; 
tbl_out.spd_lin_R2= R2; 
tbl_out.spd_lin_S_R2 = S_R2; 

writetable(tbl_out, [parent_path filesep 'PETA_tbl.csv'])

%% reproduce the Lopes et al. Fig 7.
figure(1011)
clf
subplot(2,3,1:2); cla
hold on
plot(spd_data{ex_idx}.tvec(~nan_idx) - spd_data{ex_idx}.tvec(1), spd_data{ex_idx}.spd(~nan_idx), 'k', 'linewidth', 3);
plot(spd_data{ex_idx}.tvec(~nan_idx)- spd_data{ex_idx}.tvec(1), mean(plt_mat, 2, 'omitnan'),'-', 'color', c_ord(2,:), 'linewidth', 3);
plot(spd_data{ex_idx}.tvec(~nan_idx)- spd_data{ex_idx}.tvec(1), mean(plt_s_mat, 2, 'omitnan'), 'color', [.7 .7 .7], 'linewidth', 3);
xlim([260 360])
ylim([0 80])
legend('Actual', 'Decoded', 'Shuffle', 'box',  'off')


subplot(2,3,3)
[~,~,~,p, R2_stats] = MS_bar_w_err(R2(~FS_idx' & logical(spd_mod),1)', R2(FS_idx' & logical(spd_mod),1)', [c_ord(1,:); c_ord(2,:)] , 1, 'ttest2', 1:2);
set(gca, 'xticklabel', {'Pyr Speed Cells' 'FS Speed Cells'}, 'XTickLabelRotation', 45)
ylabel('Decoding accuracy log R^2')
set(gca, 'YScale', 'log')
xlim([-0.2 4.2])
set(gca, 'YTick', [0.0001 0.001 0.01 0.1 1])
ylim([ 0.0001 1])


fprintf('<strong>Decoding accuracy: Regularly spiking Speed cells (%0.3f +/- %0.3f) vs Fastspiking speed cells (%0.3f +/- %0.3f); t(%d) = %0.3f, p = %0.3f</strong>\n', ...
    mean(R2(~FS_idx' & logical(spd_mod),1)'),  MS_SEM(R2(~FS_idx' & logical(spd_mod),1)'), mean( R2(FS_idx' & logical(spd_mod),1)'), MS_SEM( R2(FS_idx' & logical(spd_mod),1)'), ...
    R2_stats.df, R2_stats.tstat, p)
% axis square


subplot(2,3,4)
[~,~,~,p, R2_stats] = MS_bar_w_err(R2(spd_z > 1.98,1)', R2(spd_z<-1.98,1)', [c_ord(4,:)*1.5; c_ord(5,:)*2] , 1, 'ttest2', 1:2);
set(gca, 'xticklabel', {'Pos Speed Cells' 'Neg Speed Cells'}, 'XTickLabelRotation', 45)
ylabel('Decoding accuracy log R^2')
set(gca, 'YScale', 'log')
% axis square
fprintf('<strong>Decoding accuracy: Pos Speed cells (%0.3f +/- %0.3f) vs Neg-Speed cells (%0.3f +/- %0.3f); t(%d) = %0.3f, p = %0.3f</strong>\n', ...
    mean(R2(spd_z > 1.98,1)'),  MS_SEM(R2(spd_z > 1.98,1)'), mean(R2(spd_z<-1.98,1)'), MS_SEM(R2(spd_z<-1.98,1)'), ...
    R2_stats.df, R2_stats.tstat, p)
xlim([-0.2 4.2])
set(gca, 'YTick', [0.0001 0.001 0.01 0.1 1])
ylim([ 0.0001 1])

subplot(2,3,5)
[~,~,~,p, R2_stats] = MS_bar_w_err3(R2(logical(spd_mod),1)', R2(~logical(spd_mod),1)',S_R2_means', [c_ord(3,:); .7 .7 .7; .3 .3 .3], 1, 'anova1', 1:3);
set(gca, 'xticklabel', {'Speed Cells' 'Non-Speed Cells' 'shuffle'}, 'XTickLabelRotation', 45)
set(gca, 'YScale', 'log')
% axis square
fprintf('<strong>Decoding accuracy: Speed Cells (%0.3f +/- %0.3f) vs Non-Speed (%0.3f +/- %0.3f) vs shuffles (%0.3f +/- %0.3f); F(%d,%d) = %0.3f, p = %0.3f</strong>\n', ...
    mean(R2(logical(spd_mod),1)'),MS_SEM(R2(logical(spd_mod),1)'), mean(R2(~logical(spd_mod),1)'), MS_SEM(R2(~logical(spd_mod),1)'), mean(S_R2_means'), MS_SEM(S_R2_means'), ...
    R2_stats.a_tbl{2,3}, R2_stats.a_tbl{end,3}, R2_stats.stats.F, R2_stats.a_tbl{2, end})
xlim([-0.2 4.2])
set(gca, 'YTick', [0.0001 0.001 0.01 0.1 1])
ylim([ 0.0001 1])


% maximize
set(gcf,'Units','Inches');
pos = get(gcf,'Position');

set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
pause(2)

% print(gcf,[parent_path filesep 'lin_dec_R2_log.pdf'],'-dpdf','-r300', '-bestfit')
% saveas(gcf,[parent_path filesep 'R2.fig']);

%% collect the speed mod responses

% export as csv
mat_out = [];
mat_out = array2table([logical(spd_mod); round(spd_z,4); spd_p]', 'VariableNames',{'spd_mod' 'spd_zscore' 'spd_pval' });

c_type = [];
for ii  = length(g_idx):-1:1
    if g_idx(ii) ==1
        c_type{ii} = 'RS';
    elseif g_idx(ii) == 2
        c_type{ii} = 'FS';
    end
end

mat_out.R2 = round(R2, 3);
mat_out.S_R2 = round(S_R2, 3);

mat_out.cell_type = c_type';
mat_out.phase = phase';




% add in the reward modulation per cell
mat_out.rew_N_p = rew_out.p(:,1);
mat_out.rew_E_p = rew_out.p(:,2);
mat_out.rew_S_p = rew_out.p(:,3);
mat_out.rew_W_p = rew_out.p(:,4);
mat_out.rew_A_p = rew_out.p(:,5);

mat_out.app_N_p = app_out.p(:,1);
mat_out.app_E_p = app_out.p(:,2);
mat_out.app_S_p = app_out.p(:,3);
mat_out.app_W_p = app_out.p(:,4);
mat_out.app_A_p = app_out.p(:,5);


writetable(mat_out, [parent_path filesep 'Spd_mod.csv'])
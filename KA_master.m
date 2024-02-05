
%% KA screener master script.

% load data
if ismac
    addpath(genpath('/Users/jericcarmichael/Documents/Github/vandermeerlab/code-matlab/shared'))
    addpath(genpath('/Users/jericcarmichael/Documents/Github/EC_State'));
    addpath(genpath('/Users/jericcarmichael/Documents/Github/KA_analyses'));
    data_dir = '/Users/jericcarmichael/Desktop/KA_Data/for_eric_only'; % where all the NLX data is.
    % inter_dir = '/Users/jericcarmichael/Dropbox/KA_Data/inter_';  % where to save the outputs.
    inter_dir = '/Users/jericcarmichael/Desktop/KA_Data/inter_data';
    plot_dir = '/Users/jericcarmichael/Desktop/KA_Data/Behav_plots';
elseif ispc
    % load data
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'))
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\EC_State'));
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\KA_analyses'));
    data_dir = 'C:\Users\ecarm\Desktop\for_eric_only'; % where all the NLX data is.
    inter_dir = 'J:\KA_Data\inter_reward_23';
    inter_dir_app = 'J:\KA_Data\inter_reward_23_approach';
    
    
else
    addpath(genpath('/home/ecar/Github/vandermeerlab/code-matlab/shared'))
    addpath(genpath('/home/ecar/Github/EC_State'));
    addpath(genpath('/home/ecar/Github/KA_analyses'));
    data_dir = '/lustre06/project/6064766/ecar/for_eric_only'; % where all the NLX data is.
    % inter_dir = '/Users/jericcarmichael/Dropbox/KA_Data/inter_';  % where to save the outputs.
    inter_dir = '/lustre06/project/6064766/ecar/KA_Data/inter_approach_2p5_new_sig';
    
    %     % load data
    %     addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'))
    %     addpath(genpath('/home/williamslab/Documents/Github/EC_State'));
    %     addpath(genpath('/home/williamslab/Documents/Github/KA_analyses'));
    %     data_dir = '/home/williamslab/Desktop/for_eric_only'; % where all the NLX data is.
    %     inter_dir = '/media/williamslab/Fenrir/KA_Data/inter_reward_simple';
    %     inter_dir_app = '/media/williamslab/Fenrir/KA_Data/inter_reward_simple';
    
end

cd(data_dir); % move to the data dir.

% make an intermediate directory if it doesn't exist.
if ~exist(inter_dir,'dir')
    mkdir(inter_dir)
end

% if ~exist(inter_dir_app,'dir')
%     mkdir(inter_dir_app)
% end

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
    
    if isstruct(data)
        save([inter_dir filesep sess_list{iS} '_maze_data.mat'], 'data', '-v7.3')
        success(iS) = 1;
    else
        success(iS) = -10;
    end
    
    
    
end

fprintf('<strong>%0.0f total sessions, %0.2f had good cells, %0.0f omitted, %0.0f no spike data, %0.0f too short</strong>\n', length(success), sum(success==1), sum(success==99), sum(success==404), sum(success==-10))
%% Proess each cell within a session
cd(inter_dir)
sess_list = dir([inter_dir filesep '*.mat']);

phase =[]; sub = []; cell_id = [];
app_out= [];
rew_out = [];
cell_id = {};
spd_mod = [];
spd_p = [];
spd_corr = [];
stats = []; 
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
        
        % summary for plotting
        cfg_peth = [];
        cfg_peth.window = [-5 5]; 
        cfg_peth. plot_type = 'raw';
                      cfg_peth.dt = 0.05; 
            % [~, ~,this_gau] = SpikePETH_Shuff(cfg_peth, this_S, data.rew.t ); 


              % PETH for plotting
              cfg_peth = [];
              cfg_peth.window = [-5 5];
              cfg_peth. plot_type = 'raw';
              cfg_peth.dt = 0.05; 
              % [peth{k}] = KA_spike_hist(cfg_peth, this_S,data.rew.t, data.rew.in); 

        % reward centered.
        cfg_wcx_r = [];
        cfg_wcx_r.win = [-1 1]; % window
        cfg_wcx_r.baseline = [-3 -1];

        % get the response using the Wilcoxon from Frazer 2023
        [rew_out.p(k,:), rew_out.h(k,:), rew_out.base_fr(k,:), rew_out.rew_fr(k,:)] = KA_react_WCX(cfg_wcx_r, this_S, data.rew.t, data.rew.in);
        
        
        cfg_wcx_a = [];
        cfg_wcx_a.win = [-1 1]; % window
        cfg_wcx_a.baseline = [-3 -1];

        % get the response using the Wilcoxon from Frazer 2023
        [app_out.p(k,:), app_out.h(k,:), app_out.base_fr(k,:), app_out.rew_fr(k,:)] = KA_react_WCX(cfg_wcx_a, this_S, data.app.t, data.app.in);
        
        
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

%% classify cells based on waveform properties
fr = []; bur_idx = []; s_w = []; pt_r = []; rfint = []; wave_dur = []; wave_forms = []
for ii = length(stats):-1:1
    fr(ii) = stats{ii}.firing_rate; 
    bur_idx(ii) = stats{ii}.burst_idx; 
    s_w(ii) = stats{ii}.spike_width*1000; 
    pt_r(ii) = abs(stats{ii}.pt_ratio); 
    rfint(ii) = stats{ii}.rise_fall_inter; 
    wave_dur(ii) = stats{ii}.wave_dur*1000; 
    wave_forms(:,:,ii) = stats{ii}.wave; 

end
figure(808)
subplot(2, 2, 1)

[g_idx, n_val] = MS_kmean_scatter([fr', bur_idx', s_w'], 2, [1,2,3], 50); 
xlabel('Firing rate (Hz)');
ylabel('burst index')
zlabel('spike width')

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
    for jj = 1:4
        plot(this_S.waves{1}.xrange(:,jj), nanmean(wave_forms(jj+1,:, g_idx == ii),3))
    end
    title(['Group ' num2str(ii) ' | n=' num2str(sum(g_idx == ii))  ' | FR:' num2str(mean(fr(g_idx == ii)))])

end

%% check for cells with the 25ms spike width
% N_25_idx = find(s_w == .25); 
% x_range = reshape(this_S.waves{1}.xrange, 1, 32*4); 
% 
% for ii =1:length(N_25_idx)
% 
% this_idx = N_25_idx(ii); 
% 
%     figure(ii)
%     clf
%     subplot(2,2,1)
% 
% 
%     cla
%     hold on
%     for jj = 1:4
%         plot(this_S.waves{1}.xrange(:,jj), wave_forms(jj+1,:,this_idx))
%     end
% 
%     subplot(2,2,2)
%     cla
%     histogram(log10(stats{this_idx}.ISI),50,'Normalization', 'probability')
% %     set(gca, 'XScale', 'log')
%     set(gca, 'XTickLabel', (10.^get(gca, 'xtick'))*1000)
%     xlabel('ISI (ms)');
%     xline(prctile(log10(stats{this_idx}.ISI), 10), '--', '10th prctl')
%     xlim(log10([0.001 100]))
% 
%     subplot(2,2,4)
%     bar(stats{this_idx}.auto_corr.xbin, stats{this_idx}.auto_corr.ac)
% 
%     subplot(2,2,3)
%     sess_n = strsplit(cell_id{this_idx}, 'DONE_maze_data_'); 
%     text(0, .5, [strrep(sess_n{1}, '_', '-') ' ' sess_n{2}(1:end-4)])
%     text(0, .3, num2str(fr(this_idx)) )
% 
% 
% end
% 
% % %% simple plot for app vs rew FR
% % 
% % figure(108)
% % clf
% % % subplot(2,2,1)
% % % cla
% % scatter(rew_PC.rew_fr./rew_PC.base_fr, app_PC.rew_fr./app_PC.base_fr);
% % xlim([0 50])
% % ylim([0 50])


%% simple sumamry plots
c_ord = linspecer(5);
c_map = [.05 .05 0.05 ; c_ord]; %[6, 25, 34; 249, 160, 27]/255; % SUNS colour b/c KA.

c_orange = [255 150 0]/255;
c_l_orange = [255 213 153]/255;
c_purple = [98 66 158]/255;
c_l_purple = [198 183 225]/255;

% plot the significant reward response array.
figure(1010)
clf
subplot(3,4,[1 5 9])
cla
rew_sig_mat = double(rew_out.h); 
rew_sig_mat(rew_sig_mat(:,1) > 0, 1)  = 1; 
rew_sig_mat(rew_sig_mat(:,2) > 0, 2)  = 2; 
rew_sig_mat(rew_sig_mat(:,3) > 0, 3)  = 3; 
rew_sig_mat(rew_sig_mat(:,4) > 0, 4)  = 4; 
rew_sig_mat(rew_sig_mat(:,5) > 0, 5)  = 5; 

imagesc( 1:5, 1:k, rew_sig_mat)
set(gca, 'xtick', 1:5, 'XTickLabel', {'N', 'E', 'S', 'W', 'overall'})

colormap(c_map)
title(['Reward modulation (' num2str(cfg_wcx_r.win(1)) ' : ' num2str(cfg_wcx_r.win(2)) ')'])

% feeder type mod
subplot(3,4,2)
b = bar(1:6, [rew_mod, rew_no_mod]./size(rew_out.h,1)*100, 'FaceColor', 'flat');
b.CData(1:5,:) = c_ord(:,:);
b.CData(6,:) = [.6 .6 .6];

set(gca, 'xticklabel',{'North (3G)', 'West (3O)', 'South (1G)', 'East (1O)','overall', 'no mod'});
ylabel('% Reward mod');

% Flavour
subplot(3,4,6)
banana = sum(logical(sum(rew_out.h(:,1),2)) & logical(sum(rew_out.h(:,3),2)));
grain = sum(logical(sum(rew_out.h(:,2),2)) & logical(sum(rew_out.h(:,4),2)));

bf = bar(1:2, [rew_mod(1)-banana rew_mod(3)-banana, banana; rew_mod(2)-grain rew_mod(4)-grain, grain]./size(rew_out.h,1)'*100, 'stacked', 'FaceColor','flat');
bf(1).CData(1,:) = c_purple;
bf(2).CData(1,:) = c_l_purple;
bf(3).CData(1,:) = median([c_purple; c_l_purple],1);

bf(1).CData(2,:) = c_orange;
bf(2).CData(2,:) = c_l_orange;
bf(3).CData(2,:) = median([c_orange; c_l_orange],1);

xtips1 = bf(1).XEndPoints;
ytips1 = bf(1).YEndPoints;
% labels1 = string(bf(1).YData);
text(xtips1,ytips1,'high','HorizontalAlignment','center','VerticalAlignment','cap')

xtips1 = bf(2).XEndPoints;
ytips1 = bf(2).YEndPoints;
text(xtips1,ytips1,'low','HorizontalAlignment','center','VerticalAlignment','cap')

xtips1 = bf(3).XEndPoints;
ytips1 = bf(3).YEndPoints;
text(xtips1,ytips1,'both','HorizontalAlignment','center','VerticalAlignment','cap')

set(gca, 'xticklabel',{'Grape', 'Orange'});
% legend({'high', 'low', 'both'}, 'box', 'off');
ylabel('% Reward mod');

% Magnitude
subplot(3,4,10)
high = sum(logical(sum(rew_out.h(:,1),2)) & logical(sum(rew_out.h(:,2),2)));
low = sum(logical(sum(rew_out.h(:,3),2)) & logical(sum(rew_out.h(:,4),2)));

bf = bar(1:2, [rew_mod(1)-high rew_mod(2)-high, high; rew_mod(3)-low rew_mod(4)-low, low]./size(rew_out.h,1)'*100, 'stacked', 'FaceColor','flat');
set(gca, 'xticklabel',{'High', 'Low'});

ylabel('% Reward mod');

bf(1).CData(1,:) = c_purple;
bf(2).CData(1,:) = c_orange;
bf(3).CData(1,:) = median([c_purple; c_orange],1);

bf(1).CData(2,:) = c_l_purple;
bf(2).CData(2,:) = c_l_orange;
bf(3).CData(2,:) = median([c_l_purple ; c_l_orange],1);

xtips1 = bf(1).XEndPoints;
ytips1 = bf(1).YEndPoints;
% labels1 = string(bf(1).YData);
text(xtips1,ytips1,'grape','HorizontalAlignment','center','VerticalAlignment','cap')

xtips1 = bf(2).XEndPoints;
ytips1 = bf(2).YEndPoints;
text(xtips1,ytips1,'orange','HorizontalAlignment','center','VerticalAlignment','cap')

xtips1 = bf(3).XEndPoints;
ytips1 = bf(3).YEndPoints;
text(xtips1,ytips1,'both','HorizontalAlignment','center','VerticalAlignment','cap')

% legend({'high', 'low', 'both'}, 'box', 'off');
ylabel('% Reward mod');


% plot the significant approach response array.
subplot(3,4,[3 7 11])
app_h_mat = double(app_out.h);
app_h_mat(app_h_mat(:,2) == 1,2) = 2;
app_h_mat(app_h_mat(:,3) == 1,3) = 3;
app_h_mat(app_h_mat(:,4) == 1,4) = 4;
app_h_mat(app_h_mat(:,5) == 1,5) = 5;
imagesc( 1:5, 1:k, app_h_mat)
set(gca, 'xtick', 1:5, 'XTickLabel', {'N', 'E', 'S', 'W', 'overall'})

colormap(c_map)
title(['Approach modulation (' num2str(cfg_wcx_a.win(1)) ' : ' num2str(cfg_wcx_a.win(2)) ')'])


% feeder type mod
subplot(3,4,4)
b = bar(1:6, [app_mod, app_no_mod]./size(app_out.h,1)*100, 'FaceColor', 'flat');
b.CData(1:5,:) = c_ord(:,:);
b.CData(6,:) = [.6 .6 .6];

set(gca, 'xticklabel',{'North (3G)', 'West (3O)', 'South (1G)', 'East (1O)','overall', 'no mod'});
ylabel('% Approach mod');

% Flavour
subplot(3,4,8)
banana = sum(logical(sum(app_out.h(:,1),2)) & logical(sum(app_out.h(:,3),2)));
grain = sum(logical(sum(app_out.h(:,2),2)) & logical(sum(app_out.h(:,4),2)));

bf = bar(1:2, [app_mod(1)-banana app_mod(3)-banana, banana; app_mod(2)-grain app_mod(4)-grain, grain]./size(app_out.h,1)'*100, 'stacked', 'FaceColor','flat');
bf(1).CData(1,:) = c_purple;
bf(2).CData(1,:) = c_l_purple;
bf(3).CData(1,:) = median([c_purple; c_l_purple],1);

bf(1).CData(2,:) = c_orange;
bf(2).CData(2,:) = c_l_orange;
bf(3).CData(2,:) = median([c_orange; c_l_orange],1);


xtips1 = bf(1).XEndPoints;
ytips1 = bf(1).YEndPoints;
% labels1 = string(bf(1).YData);
text(xtips1,ytips1,'high','HorizontalAlignment','center','VerticalAlignment','cap')

xtips1 = bf(2).XEndPoints;
ytips1 = bf(2).YEndPoints;
text(xtips1,ytips1,'low','HorizontalAlignment','center','VerticalAlignment','cap')

xtips1 = bf(3).XEndPoints;
ytips1 = bf(3).YEndPoints;
text(xtips1,ytips1,'both','HorizontalAlignment','center','VerticalAlignment','cap')

set(gca, 'xticklabel',{'Grape', 'Orange'});
% legend({'high', 'low', 'both'}, 'box', 'off');
ylabel('% Approach mod');

% Magnitude
subplot(3,4,12)
high = sum(logical(sum(app_out.h(:,1),2)) & logical(sum(app_out.h(:,2),2)));
low = sum(logical(sum(app_out.h(:,3),2)) & logical(sum(app_out.h(:,4),2)));

bf = bar(1:2, [app_mod(1)-high app_mod(2)-high, high; app_mod(3)-low app_mod(4)-low, low]./size(app_out.h,1)'*100, 'stacked', 'FaceColor','flat');
set(gca, 'xticklabel',{'High', 'Low'});
% legend({'banana', 'grain', 'both'}, 'box', 'off');
ylabel('% Approach mod');

bf(1).CData(1,:) = c_purple;
bf(2).CData(1,:) = c_orange;
bf(3).CData(1,:) = median([c_purple; c_orange],1);

bf(1).CData(2,:) = c_l_purple;
bf(2).CData(2,:) = c_l_orange;
bf(3).CData(2,:) = median([c_l_purple ; c_l_orange],1);

xtips1 = bf(1).XEndPoints;
ytips1 = bf(1).YEndPoints;
% labels1 = string(bf(1).YData);
text(xtips1,ytips1,'grape','HorizontalAlignment','center','VerticalAlignment','cap')

xtips1 = bf(2).XEndPoints;
ytips1 = bf(2).YEndPoints;
text(xtips1,ytips1,'orange','HorizontalAlignment','center','VerticalAlignment','cap')

xtips1 = bf(3).XEndPoints;
ytips1 = bf(3).YEndPoints;
text(xtips1,ytips1,'both','HorizontalAlignment','center','VerticalAlignment','cap')

% legend({'high', 'low', 'both'}, 'box', 'off');
ylabel('% Approach mod');


%% split data into session phases
P = [];
P.C_idx = contains(phase, 'C');
P.OE_idx = contains(phase, 'O1') | contains(phase, 'O2') | contains(phase, 'O3');
P.OL_idx = contains(phase, 'O4') | contains(phase, 'O5') | contains(phase, 'O6') | contains(phase, 'O7');
P.R_idx = contains(phase, 'R');
% P.E_idx = contains(phase, 'E');

p_list = fieldnames(P);

s_idx = reshape(1:length(p_list)*2,2, length(p_list))'; % for subplots

figure(1019)
clf

for iP = 1:length(p_list)
    this_rew_mod = [];
    this_app_mod = [];
    for ii = size(rew_out.h,2):-1:1
        sig_rew = rew_out.h(P.(p_list{iP})); 
        this_rew_mod(ii) = sum(sum(rew_out.h(P.(p_list{iP}),ii),2)>0);

        pos_mod = (rew_out.base_fr(P.(p_list{iP}), ii) < rew_out.rew_fr(P.(p_list{iP}),ii))   & rew_out.h(P.(p_list{iP}),ii); 
        neg_mod = (rew_out.base_fr(P.(p_list{iP}), ii) > rew_out.rew_fr(P.(p_list{iP}),ii))   & rew_out.h(P.(p_list{iP}),ii); 

        this_rew_mod_pos(ii) = sum(sum(rew_out.h(P.(p_list{iP}),ii),2)>0);


        this_app_mod(ii) = sum(sum(app_out.h(P.(p_list{iP}),ii),2)>0);

    end
    
    this_rew_no_mod = sum(sum(rew_out.h(P.(p_list{iP}),1:5),2) == 0);
    this_app_no_mod = sum(sum(app_out.h(P.(p_list{iP}),1:5),2) == 0);
    
    
    subplot(length(p_list),2,s_idx(iP,1))
    
    b = bar(1:6, ([this_rew_mod, this_rew_no_mod]./ size(rew_out.h(P.(p_list{iP})),2))*100, 'FaceColor', 'flat');
    b.CData(1:5,:) = c_ord(:,:);
    b.CData(6,:) = [.6 .6 .6];
    ylim([0 100])
    
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(round(b(1).YData));
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','top')
    
    set(gca, 'xticklabel',{'North (3G)', 'West (3O)', 'South (1G)', 'East (1O)','all rew', 'no mod'});
    ylabel('% Reward mod');
    title([p_list{iP}(1:end-4) ' (n = ' num2str(sum(P.(p_list{iP}))) ')'])
    
    
    subplot(length(p_list),2,s_idx(iP,2))
    
    b = bar(1:6, ([this_app_mod, this_app_no_mod]./ size(app_out.h(P.(p_list{iP})),2))*100, 'FaceColor', 'flat');
    b.CData(1:5,:) = c_ord(:,:);
    b.CData(6,:) = [.6 .6 .6];
    ylim([0 100])
    
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(round(b(1).YData));
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','top')
    
    set(gca, 'xticklabel',{'North (3G)', 'West (3O)', 'South (1G)', 'East (1O)','all app', 'no mod'});
    ylabel('% Approach mod');
    title([p_list{iP}(1:end-4) ' (n = ' num2str(sum(P.(p_list{iP}))) ')'])
end


%%
p_types = unique(phase);
Phase = [];rew_mat = []; app_mat = []; 

for iP = length(p_types):-1:1
    this_idx = strcmpi(p_types{iP}, phase);
    for ii = size(rew_out.h,2):-1:1
        Phase.rew_mod{iP}(:,ii) = sum(sum(rew_out.h(this_idx,ii),2)>0)/sum(this_idx);
        Phase.app_mod{iP}(:,ii) = sum(sum(app_out.h(this_idx,ii),2)>0)/sum(this_idx);
    end
    rew_mat(iP,:) = Phase.rew_mod{iP}; 
    app_mat(iP,:) = Phase.app_mod{iP};
%     labels{iP} = [p_types{iP}(1) str2double(p_types{iP}(end))]); 
end

figure(1099)
clf

subplot(1,2,1)
imagesc(1:5, 1:length(p_types), rew_mat*100); 
set(gca, 'xtick', [1:5], 'XTickLabel', {'North', 'West', 'South', 'East', 'overall'});
set(gca,'YTick', 1:length(p_types),  'YTickLabel', p_types)
caxis([0 75]); c = colorbar('Location', 'eastoutside');
c.Ticks = [0 25 50 75]; c.Label.String = '% modulated cells'; 
title('Reward')

subplot(1,2,2); cla;
imagesc(1:5, 1:length(p_types), app_mat*100); 
set(gca, 'xtick', [1:5], 'XTickLabel', {'North', 'West', 'South', 'East', 'overall'});
set(gca,'YTick', 1:length(p_types),  'YTickLabel', p_types)
caxis([0 75]); c = colorbar('Location', 'eastoutside');
c.Ticks = [0 25 50 75]; c.Label.String = '% modulated cells'; 
title('Approach')


%% collect the percentage response
for ii = size(rew_out.h,1):-1:1

    for jj = size(rew_out.rew_fr, 2):-1:1
        rew_prct(ii,jj) = rew_out.rew_fr(ii,jj)./ rew_out.base_fr(ii,jj);

        app_prct(ii,jj) = app_out.rew_fr(ii,jj) / app_out.base_fr(ii,jj);

    end
end

C = rew_prct(P.C_idx,5);
O_e = rew_prct(P.OE_idx,5);
O_l = rew_prct(P.OL_idx,5);
R = rew_prct(P.R_idx,5);

% export as csv
mat_out = NaN(4, max([length(C),length(O_e),length(O_l),length(R)]));


mat_out(1,1:length(C)) = C';
mat_out(2,1:length(O_e)) = O_e';
mat_out(3,1:length(O_l)) = O_l';
mat_out(4,1:length(R)) = R';

mat_out= mat_out';

csvwrite([parent_path filesep 'all_Rew_perc.csv'], mat_out)


aC = app_prct(P.C_idx,5);
aO_e = app_prct(P.OE_idx,5);
aO_l = app_prct(P.OL_idx,5);
aR = app_prct(P.R_idx,5);


% export as csv
mat_out = NaN(4, max([length(aC),length(aO_e),length(aO_l),length(aR)]));


mat_out(1,1:length(aC)) = aC';
mat_out(2,1:length(aO_e)) = aO_e';
mat_out(3,1:length(aO_l)) = aO_l';
mat_out(4,1:length(aR)) = aR';

mat_out= mat_out';

csvwrite([parent_path filesep 'all_App_perc.csv'], mat_out)



%% collect the percentage response for significantly modulated cells
sig_idx = (rew_out.h(:,5) ==1)';
C_sig  = rew_out.rew_fr(P.C_idx & sig_idx,5)./ rew_out.base_fr(P.C_idx & sig_idx,5);
O_e_sig  = rew_out.rew_fr(P.OE_idx & sig_idx,5) ./ rew_out.base_fr(P.OE_idx & sig_idx,5);
O_l_sig  = rew_out.rew_fr(P.OL_idx & sig_idx,5) ./ rew_out.base_fr(P.OL_idx & sig_idx,5);
R_sig  = rew_out.rew_fr(P.R_idx & sig_idx,5) ./ rew_out.base_fr(P.R_idx & sig_idx,5);

% export as csv
mat_out = NaN(4, max([length(Cs_sig),length(O_e_sig),length(O_l_sig),length(R_sig)]));


mat_out(1,1:length(C_sig)) = C_sig';
mat_out(2,1:length(O_e_sig)) = O_e_sig';
mat_out(3,1:length(O_l_sig)) = O_l_sig';
mat_out(4,1:length(R_sig)) = R_sig';

mat_out= mat_out';

csvwrite([parent_path filesep 'Sig_Rew_perc.csv'], mat_out)


sig_idx = (app_out.h(:,5) ==1)';
aC_sig  = app_out.rew_fr(P.C_idx & sig_idx,5) ./ app_out.base_fr(P.C_idx & sig_idx,5);
aO_e_sig  = app_out.rew_fr(P.OE_idx & sig_idx,5) ./ app_out.base_fr(P.OE_idx & sig_idx,5);
aO_l_sig  = app_out.rew_fr(P.OL_idx & sig_idx,5) ./ app_out.base_fr(P.OL_idx & sig_idx,5);
aR_sig  = app_out.rew_fr(P.R_idx & sig_idx,5) ./ app_out.base_fr(P.R_idx & sig_idx,5);

% export as csv
mat_out = NaN(4, max([length(aCs_sig),length(aO_e_sig),length(aO_l_sig),length(aR_sig)]));


mat_out(1,1:length(aC_sig)) = aC_sig';
mat_out(2,1:length(aO_e_sig)) = aO_e_sig';
mat_out(3,1:length(aO_l_sig)) = aO_l_sig';
mat_out(4,1:length(aR_sig)) = aR_sig';

mat_out= mat_out';

csvwrite([parent_path filesep 'Sig_App_perc.csv'], mat_out)

%% collect all response % 
for ii = size(rew_out.h,1):-1:1

    for jj = size(rew_out.rew_fr, 2):-1:1
        rew_prct(ii,jj) = rew_out.rew_fr(ii,jj)
        
        z_base = (rew_out.base_fr(ii,jj) - stats{ii}.firing_rate) / stats{ii}.std     


    end
end

C = rew_prct(P.C_idx,5);
O_e = rew_prct(P.OE_idx,5);
O_l = rew_prct(P.OL_idx,5);
R = rew_prct(P.R_idx,5);

% export as csv
mat_out = NaN(4, max([length(C),length(O_e),length(O_l),length(R)]));


mat_out(1,1:length(C)) = C';
mat_out(2,1:length(O_e)) = O_e';
mat_out(3,1:length(O_l)) = O_l';
mat_out(4,1:length(R)) = R';

mat_out= mat_out';

csvwrite([parent_path filesep 'all_Rew_perc.csv'], mat_out)



%% collect the percentage response for significantly modulated cells
sig_idx = (rew_out.h(:,5) ==1)';
C_sig  = rew_out.rew_fr(P.C_idx & sig_idx,5)./ rew_out.base_fr(P.C_idx & sig_idx,5);
O_e_sig  = rew_out.rew_fr(P.OE_idx & sig_idx,5) ./ rew_out.base_fr(P.OE_idx & sig_idx,5);
O_l_sig  = rew_out.rew_fr(P.OL_idx & sig_idx,5) ./ rew_out.base_fr(P.OL_idx & sig_idx,5);
R_sig  = rew_out.rew_fr(P.R_idx & sig_idx,5) ./ rew_out.base_fr(P.R_idx & sig_idx,5);

% export as csv
mat_out = NaN(4, max([length(Cs_sig),length(O_e_sig),length(O_l_sig),length(R_sig)]));


mat_out(1,1:length(C_sig)) = C_sig';
mat_out(2,1:length(O_e_sig)) = O_e_sig';
mat_out(3,1:length(O_l_sig)) = O_l_sig';
mat_out(4,1:length(R_sig)) = R_sig';

mat_out= mat_out';

csvwrite([parent_path filesep 'Sig_Rew_perc.csv'], mat_out)












%%

%     for iT = 1:length(cells_to_process)
%
%         parts = strsplit(cells_to_process{iT}, filesep);
%         this_file = parts{end};
% %         This_cell = KA_screener_zscore(this_file);
% %         This_cell = KA_screener_approach(this_file);
% %         This_cell = KA_screener(this_file);
% %         This_cell = KA_screener_feeder(this_file);
% %         This_cell = KA_screener_v2(this_file);
%         This_cell = KA_screener_simple(this_file);
%
%
% %         This_cell = KA_screener_pseudo_baseline(this_file);
%
%         % if there were too few spikes in the .t then skip this file.
%         if ischar(This_cell) || isempty(This_cell)
%             fprintf('<strong>%s</strong>: Minimum requirments not met: %s  -   <strong>%s</strong>\n', mfilename, this_file, This_cell)
%             if strcmpi(This_cell, 'too short')
%                 success(iS) = -1;
%                 continue
%             end
%         elseif isnumeric(This_cell)
%                 success(iS) = -3;
%                 FR(length(success)) = This_cell;
%             continue
%         end
%         success(iS) = 1;
%
%         FR(length(success)) = length(This_cell.S.t{1})/(This_cell.pos.tvec(end) - This_cell.pos.tvec(1));
%
%         parts = strsplit(sess_list{iS}, '_');
%         This_cell.subject = [parts{1} '_' parts{2}]; % get the subject ID
%         This_cell.session = parts{3}; % get the session type and number Acquisition, Criteria, Overtrain, Extinction, Reacquisition. nSessions: A ? C 3, O 7, E 1, R 3-5
%         This_cell.date = parts{end};
%
%         save([inter_dir filesep sess_list{iS} '_' this_file(1:strfind(this_file, '.')-1) '_Feeder.mat'], 'This_cell')
%
%         close all
%     end
%
% end
%
% % summarize the files
% fprintf('<strong>Total Sessions: %2.0f\nnSucess: %2.0f (%2.2f%%)\nToo short: %2.0f (%2.2f%%)\nFR too low: %2.0f (%2.2f%%)\n</strong>',...
%     length(success),sum(success==1), ((sum(success==1))/length(success))*100,sum(success==2), ((sum(success==2))/length(success))*100,...
%     sum(success==3),((sum(success==3))/length(success))*100)
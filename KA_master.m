
%% KA screener master script.

% load data
if ismac

usr_name = char(java.lang.System.getProperty('user.name')); 

    addpath(genpath(['/Users/' usr_name '/Documents/Github/vandermeerlab/code-matlab/shared']))
    addpath(genpath(['/Users/' usr_name '/Documents/Github/EC_State']));
    addpath(genpath(['/Users/' usr_name '/Documents/Github/KA_analyses']));
        % addpath(genpath(['/Users/' usr_name '/Documents/Github/NeuroLearning/M1_compiled_loaders_functions/releaseDec2015/binaries'])); % Apple silicon compiled Nlx Loaders

        data_dir = ['/Users/' usr_name '/Desktop/KA_Data/for_eric_only']; % where all the NLX data is.
    % inter_dir = '/Users/jericcarmichael/Dropbox/KA_Data/inter_';  % where to save the outputs.
    inter_dir = ['/Users/' usr_name '/Desktop/KA_Data/inter_data'];
    plot_dir = ['/Users/' usr_name '/Desktop/KA_Data/Behav_plots'];
elseif ispc
    % load data
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'))
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\EC_State'));
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\KA_analyses'));
    data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\for_eric_only'; % where all the NLX data is.
    inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\KA_Data\inter_reward_23';
    inter_dir_app = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\KA_Data\inter_reward_23_approach';
    plot_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\KA_Data\Behav_plots';


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

f_id = {'North', 'West', 'South',  'East', 'Overall'};
flav_id = {'Grape x3','Orange x3', 'Grape x1', 'Orange x1', '   '};

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
spd_data= []; 
k = 0;
for iS = 1:length(sess_list)

    load([inter_dir filesep sess_list(iS).name])

    % get the mean velocity when the animal is moving.
    % mVelo(iS) = mean(data.velo_smooth.data(data.velo_smooth.data>5));

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


        % summary for plotting
        cfg_peth = [];
        cfg_peth.window = [-6 6];
        cfg_peth. plot_type = 'raw';
        cfg_peth.dt = 0.05;
        cfg_peth.gauss_sd = .2;
        for jj = unique(data.rew.in)
            this_idx = data.rew.in == jj;
            [~,~,this_peth] = SpikePETH_Shuff(cfg_peth, this_S, data.rew.t(this_idx) );
            all_peth(:,jj, k) = nanmean(this_peth,2);
        end

        [~,tvec,this_peth] = SpikePETH_Shuff(cfg_peth, this_S, data.rew.t );
        all_peth(:,5, k) = nanmean(this_peth,2);

        app_t{iS}= data.app; 


        % PETH for plotting
        % cfg_peth = [];
        % cfg_peth.window = [-5 5];
        % cfg_peth. plot_type = 'raw';
        % cfg_peth.dt = 0.05;
        % [peth{k}] = KA_spike_hist(cfg_peth, this_S,data.rew.t, data.rew.in);

        % reward centered.
        cfg_wcx_r = [];
        cfg_wcx_r.win = [0 2]; % window
        cfg_wcx_r.baseline = [-3 -1];

        % get the response using the Wilcoxon from Frazer 2023
        [rew_out.p(k,:), rew_out.h(k,:), rew_out.base_fr(k,:), rew_out.rew_fr(k,:)] = KA_react_WCX(cfg_wcx_r, this_S, data.rew.t, data.rew.in);


        cfg_wcx_a = [];
        cfg_wcx_a.win = cfg_wcx_r.win ; % window
        cfg_wcx_a.baseline = cfg_wcx_r.baseline;

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

save([parent_path filesep 'Spd_data.mat'], 'spd_data')
%% classify cells based on waveform properties
fr = []; bur_idx = []; s_w = []; pt_r = []; rfint = []; wave_dur = []; wave_forms = []; 
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
        shadedErrorBar(this_S.waves{1}.xrange(:,jj), nanmean(wave_forms(jj+1,:, g_idx == ii),3),std(wave_forms(jj+1,:, g_idx == ii),[],3) )
        plot(this_S.waves{1}.xrange(:,jj), nanmean(wave_forms(jj+1,:, g_idx == ii),3))
        % errorbar(this_S.waves{1}.xrange(:,jj), nanmean(wave_forms(jj+1,:, g_idx == ii),3) +std(wave_forms(jj+1,:, g_idx == ii),[],3))
    end
    title(['Group ' num2str(ii) ' | n=' num2str(sum(g_idx == ii))  ' | FR:' num2str(round(mean(fr(g_idx == ii)),2)) '\pm' num2str(round(std(fr(g_idx == ii)),2))])

end

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,[parent_path filesep 'wave_feature.pdf'],'-dpdf','-r300')



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

    fprintf('%s modulation: Rew positive %0.2f%%  Rew negative %0.2f%%  No mod %0.2f%%\n',p_list{iP},...
        ((this_rew_mod_pos(5)-sum(neg_mod))./length(pos_mod))*100,...
        ((this_rew_mod_pos(5)-sum(pos_mod))./length(pos_mod))*100, (this_rew_no_mod/length(pos_mod))*100)

    subplot(length(p_list)+1,2,s_idx(iP,1))

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


    subplot(length(p_list)+1,2,s_idx(iP,2))

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

% overall
this_rew_mod = [];
this_app_mod = [];
this_rew_mod_pos = [];
this_rew_mod_neg = [];
this_rew_no_mod = [];

for ii = size(rew_out.h,2):-1:1
    sig_rew = rew_out.h(:,ii);


    pos_mod = (rew_out.base_fr(:, ii) < rew_out.rew_fr(:,ii))   & sig_rew;
    neg_mod = (rew_out.base_fr(:, ii) > rew_out.rew_fr(:,ii))   & sig_rew;

    this_rew_mod_pos(ii) = sum(pos_mod>0)./length(sig_rew);
    this_rew_mod_neg(ii) = sum(neg_mod>0)./length(sig_rew);
    this_rew_mod(ii) = sum(sig_rew>0)./length(sig_rew);
end

this_rew_no_mod = sum(rew_out.h(:,5)==0)./length(rew_out.h(:,5));


subplot(length(p_list)+1,2,((length(p_list)+1)*2)-1)

b = bar(1:6, [this_rew_mod, this_rew_no_mod]*100, 'FaceColor', 'flat');
b.CData(1:5,:) = c_ord(:,:);
b.CData(6,:) = [.6 .6 .6];
ylim([0 100])

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData));
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','top')

set(gca, 'xticklabel',{'North (3G)', 'West (3O)', 'South (1G)', 'East (1O)','all rew', 'no mod'});
ylabel('% Reward mod');
title(['All cells (n = ' num2str(length(rew_out.h(:,5))) ')'])

fprintf('Overall modulation: Rew positive %0.2f%%  Rew negative %0.2f%%  No mod %0.2f%%\n', this_rew_mod_pos(5)*100, this_rew_mod_neg(5)*100, this_rew_no_mod*100)

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
k = 0; 
rew_prct = NaN(size(rew_out.h,1), size(rew_out.rew_fr, 2)); 
app_prct = rew_prct; 
sub_id = cell(1,size(rew_out.h,1)); 
rew_all = NaN(size(rew_out.h,1)*(size(rew_out.rew_fr, 2)-1),1); 
rew_id = cell(size(rew_out.h,1)*(size(rew_out.rew_fr, 2)-1),1); 
rew_feeder = rew_id; 
rew_phase = rew_id; 
R_idx = rew_all; 
for ii = 1:size(rew_out.h,1)

    for jj = 1:size(rew_out.rew_fr, 2)
        rew_prct(ii,jj) = rew_out.rew_fr(ii,jj)./ rew_out.base_fr(ii,jj);

        app_prct(ii,jj) = app_out.rew_fr(ii,jj) / app_out.base_fr(ii,jj);

        sub_id{ii} = cell_id{ii}(1:4);

        if jj~=5
                    k = k+1;

            rew_all(k) = rew_out.rew_fr(ii,jj)./ rew_out.base_fr(ii,jj);
            rew_id{k} = cell_id{ii}(1:4); 
            rew_feeder{k} = f_id{jj}; 
            rew_phase{k} = cell_id{ii}(6:7); 
        end
        if contains(cell_id{ii}, '_R')
            R_idx(k) = 1; 
        else
            R_idx(k) = 0; 
        end
    end
end

R_idx = logical(R_idx); 

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

R_tab = table(rew_id(R_idx),rew_phase(R_idx), rew_feeder(R_idx), rew_all(R_idx), 'VariableNames', {'ID','Phase' 'Arm', 'Response'}); 

writetable(R_tab,[parent_path filesep 'R_only_table.csv'] )



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
%% collect responses based on feeder locations

for ii = 1:4

C = rew_prct(P.C_idx,ii);
O_e = rew_prct(P.OE_idx,ii);
O_l = rew_prct(P.OL_idx,ii);
R = rew_prct(P.R_idx,ii);

% export as csv
mat_out = NaN(4, max([length(C),length(O_e),length(O_l),length(R)]));


mat_out(1,1:length(C)) = C';
mat_out(2,1:length(O_e)) = O_e';
mat_out(3,1:length(O_l)) = O_l';
mat_out(4,1:length(R)) = R';

mat_out= mat_out';

csvwrite([parent_path filesep 'all_Rew_perc_' f_id{ii} '.csv'], mat_out)

end

%% collect the percentage response for significantly modulated cells
sig_idx = (rew_out.h(:,5) ==1)';
C_sig  = rew_out.rew_fr(P.C_idx & sig_idx,5)./ rew_out.base_fr(P.C_idx & sig_idx,5);
O_e_sig  = rew_out.rew_fr(P.OE_idx & sig_idx,5) ./ rew_out.base_fr(P.OE_idx & sig_idx,5);
O_l_sig  = rew_out.rew_fr(P.OL_idx & sig_idx,5) ./ rew_out.base_fr(P.OL_idx & sig_idx,5);
R_sig  = rew_out.rew_fr(P.R_idx & sig_idx,5) ./ rew_out.base_fr(P.R_idx & sig_idx,5);

% export as csv
mat_out = NaN(4, max([length(C_sig),length(O_e_sig),length(O_l_sig),length(R_sig)]));


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
mat_out = NaN(4, max([length(aC_sig),length(aO_e_sig),length(aO_l_sig),length(aR_sig)]));


mat_out(1,1:length(aC_sig)) = aC_sig';
mat_out(2,1:length(aO_e_sig)) = aO_e_sig';
mat_out(3,1:length(aO_l_sig)) = aO_l_sig';
mat_out(4,1:length(aR_sig)) = aR_sig';

mat_out= mat_out';

csvwrite([parent_path filesep 'Sig_App_perc.csv'], mat_out)

% collect all response %
for ii = size(rew_out.h,1):-1:1

    for jj = size(rew_out.rew_fr, 2):-1:1

        z_rew = (rew_out.rew_fr(ii,jj) - stats{ii}.firing_rate) / stats{ii}.std;

        z_base = (rew_out.base_fr(ii,jj) - stats{ii}.firing_rate) / stats{ii}.std     ;

        rew_idx(ii,jj) = z_rew - z_base;
    end
end

C = rew_idx(P.C_idx,5);
O_e = rew_idx(P.OE_idx,5);
O_l = rew_idx(P.OL_idx,5);
R = rew_idx(P.R_idx,5);

% export as csv
mat_out = NaN(4, max([length(C),length(O_e),length(O_l),length(R)]));


mat_out(1,1:length(C)) = C';
mat_out(2,1:length(O_e)) = O_e';
mat_out(3,1:length(O_l)) = O_l';
mat_out(4,1:length(R)) = R';

mat_out= mat_out';

csvwrite([parent_path filesep 'all_Rew_idx.csv'], mat_out)


% same thing but only the signficant ones.

sig_idx = (rew_out.h(:,5) ==1)';

Cs = rew_idx(P.C_idx & sig_idx,5);
O_es = rew_idx(P.OE_idx & sig_idx,5);
O_ls = rew_idx(P.OL_idx & sig_idx,5);
Rs = rew_idx(P.R_idx & sig_idx,5);

% export as csv
mat_out_s = NaN(4, max([length(Cs),length(O_es),length(O_ls),length(Rs)]));


mat_out_s(1,1:length(Cs)) = Cs';
mat_out_s(2,1:length(O_es)) = O_es';
mat_out_s(3,1:length(O_ls)) = O_ls';
mat_out_s(4,1:length(Rs)) = Rs';

mat_out_s= mat_out_s';

csvwrite([parent_path filesep 'all_Rew_sig_idx.csv'], mat_out_s)

nRew_S_pos = sum(mat_out_s > 0, 'all') ./ length(rew_out.h);
nRew_S_neg = sum(mat_out_s < 0, 'all') ./ length(rew_out.h);

% collect the percentage response for significantly modulated cells
sig_idx = (rew_out.h(:,5) ==1)';
C_sig  = rew_out.rew_fr(P.C_idx & sig_idx,5)./ rew_out.base_fr(P.C_idx & sig_idx,5);
O_e_sig  = rew_out.rew_fr(P.OE_idx & sig_idx,5) ./ rew_out.base_fr(P.OE_idx & sig_idx,5);
O_l_sig  = rew_out.rew_fr(P.OL_idx & sig_idx,5) ./ rew_out.base_fr(P.OL_idx & sig_idx,5);
R_sig  = rew_out.rew_fr(P.R_idx & sig_idx,5) ./ rew_out.base_fr(P.R_idx & sig_idx,5);

% export as csv
mat_out = NaN(4, max([length(C_sig),length(O_e_sig),length(O_l_sig),length(R_sig)]));


mat_out(1,1:length(C_sig)) = C_sig';
mat_out(2,1:length(O_e_sig)) = O_e_sig';
mat_out(3,1:length(O_l_sig)) = O_l_sig';
mat_out(4,1:length(R_sig)) = R_sig';

mat_out= mat_out';

csvwrite([parent_path filesep 'Sig_Rew_perc.csv'], mat_out)







%% %%%%% the 'Good Figure' %%%%%%%%%%%

%prepare the data
all_p_id = NaN(size(cell_id));
% convert the cell_id to a numerical
for ii = length(cell_id):-1:1
    P_id{ii}=  cell_id{ii}(6:7);
end

S.C_idx = contains(P_id, 'C1') | contains(P_id, 'C2') | contains(P_id, 'C3');
S.OE_idx = contains(P_id, 'O1') | contains(P_id, 'O2') | contains(P_id, 'O3');
S.OL_idx = contains(P_id, 'O4') | contains(P_id, 'O5') | contains(P_id, 'O6') | contains(P_id, 'O7');
S.R_idx = contains(P_id, 'R');

lump_p_id = all_p_id  ;
lump_p_id(S.C_idx) = 1;
lump_p_id(S.OE_idx) = 2;
lump_p_id(S.OL_idx) = 3;
lump_p_id(S.R_idx) = 4;

all_p_id(contains(P_id, 'C1')) = 1;
all_p_id(contains(P_id, 'C2')) = 2;
all_p_id(contains(P_id, 'C3')) = 3;
all_p_id(contains(P_id, 'O1')) = 4;
all_p_id(contains(P_id, 'O2')) = 5;
all_p_id(contains(P_id, 'O3')) = 6;
all_p_id(contains(P_id, 'O4')) = 7;
all_p_id(contains(P_id, 'O5')) = 8;
all_p_id(contains(P_id, 'O6')) = 9;
all_p_id(contains(P_id, 'O7')) = 10;
all_p_id(contains(P_id, 'R1')) = 11;
all_p_id(contains(P_id, 'R2')) = 12;
all_p_id(contains(P_id, 'R3')) = 13;

[s_all_p_id, sort_idx] = sort(all_p_id);

s_P_id = P_id(sort_idx);

s_d_idx = logical([1 (diff(s_all_p_id) > 0)]);



s_lump_p_id = lump_p_id(sort_idx);
s_cell_id = cell_id(sort_idx);

s_g_idx = g_idx(sort_idx);

s_H = rew_out.h(sort_idx,:);
s_rFR = rew_out.rew_fr(sort_idx,:);
s_bFR = rew_out.base_fr(sort_idx,:);

t = linspecer(5);

c_ord = linspecer(5); % one for each session type C, O, E, R.
c_ord(2,:) =  c_ord(4,:);
c_ord(3,:) = t(2,:);
c_ord(4,:) = t(5,:);

blues = parula(16); reds = jet(64); oranges = autumn(16);
sess_cord = [flipud(blues(3:2:7,:));(reds(end-7:end-1,:)); c_ord(3,:); flipud(oranges(end-9:end-6,:))];


%% z score the peth gaussian outputs
for k = size(all_peth,3):-1:1
    for ii  = 1:5
        z_peth(:,ii,k) = (all_peth(:,ii, k) - mean(stats{k}.firing_rate))./ stats{k}.std;
    end
end

z_peth_s = z_peth(:,:, sort_idx);

Z_sig_min = -5;
Z_sig_max = 5;

f_id = {'North', 'West', 'South',  'East', 'Overall'};
flav_id = {'Grape x3','Orange x3', 'Grape x1', 'Orange x1', '   '};

%% plot
figure(999)
cla

for ii = 1:5

    subplot(1,5,ii)

    imagesc(tvec, 1:size(all_peth,3), squeeze(z_peth_s(:,ii,:))')

    title({f_id{ii}; flav_id{ii}}, 'fontsize', 20)
    set(gca, 'ytick', [])
    set(gca, 'xtick', -5:2.5:5)
    % set(gca,'xticklabel', [get(gca, 'XTicklabel') ; num2str(tvec(end),0)]   )
    vl = xline(0, '-k');
    vl.LineWidth = 2;

    vl = xline(cfg_wcx_r.win(1), '--r');
    vl.LineWidth = 2;
    vl = xline(cfg_wcx_r.win(2), '--r');
    vl.LineWidth = 2;

    vl = xline(cfg_wcx_r.baseline(1), '--r');
    vl.LineWidth = 2;
    vl = xline(cfg_wcx_r.baseline(2), '--r');
    vl.LineWidth = 2;

    caxis([Z_sig_min Z_sig_max]);


    hl = hline(find(diff(s_lump_p_id))+.5, {'k', 'k', 'k'});
    for kk = 1:length(hl)
        hl(kk).LineWidth = 3;
        hl(kk).Color = c_ord(kk+1,:);
    end

    % add in sig markers
    for jj = 1:length(s_H)
        if s_d_idx(jj) && ii == 1
            text(tvec(1)-mode(diff(tvec))*1, jj, s_P_id{jj}, 'fontweight', 'bold','fontsize', 20, 'color', c_ord(s_lump_p_id(jj),:), 'HorizontalAlignment', 'right', 'VerticalAlignment','cap');

        end

        if s_H(jj, ii) && (s_bFR(jj,ii) > s_rFR(jj,ii))
            text(tvec(end)+mode(diff(tvec))*5, jj, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(s_lump_p_id(jj),:), 'HorizontalAlignment', 'left','Interpreter','tex');
        elseif s_H(jj, ii) && (s_bFR(jj,ii) < s_rFR(jj,ii))
            text(tvec(end)+mode(diff(tvec))*10, jj, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(s_lump_p_id(jj),:), 'HorizontalAlignment', 'left','Interpreter','tex');
        end
        if s_g_idx(jj) == 2
            text(tvec(end)+mode(diff(tvec))*15, jj, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(s_lump_p_id(jj),:), 'HorizontalAlignment', 'left','Interpreter','tex');
        end
    end
    xlim([-5 5])
end
cb=colorbar;
cb.Position(1) = cb.Position(1) + .075;
cb.Label.String = 'Mean activity (zscore)';
cb.Label.FontSize = 12;
cb.FontSize = 12;
maximize
pause(2)
%%
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,[parent_path filesep 'Good_figure.pdf'],'-dpdf','-r300')



%% save the mean velocity per condition
s_phase= []; mVelo = []; all_velo_mean = []; all_velo_max = []; 

thresh = 2.5;

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



sess_list(find(contains(sess_list, 'C6_3_E1_2021-10-13_DONE'))) = [];

for iS = length(sess_list):-1:1
    cd([data_dir filesep sess_list{iS}])


    disp([ num2str(iS) '_' sess_list{iS}])
    s_phase{iS}= sess_list{iS}(6:7);

    if ~isempty(dir('*VT*.zip')) && isempty(dir('*.nvt'))
        unzip('VT1.zip')
    end

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
        mVelo(iS) = NaN;
        continue
    end


    cfg_pos.convFact = [560/142 480/142];
    data.pos = LoadPos(cfg_pos);
    data.pos = restrict(data.pos, evt.t{s_rec_idx}(task_rec_idx_s), evt.t{e_rec_idx}(task_rec_idx_e)); % restrict position to time on track.

    cfg_spd.verbose = 0;
    data.velo = getLinSpd(cfg_spd, data.pos);
    data.velo.data = interp1(data.velo.tvec,data.velo.data(1,:),data.pos.tvec,'linear');

    % smooth speed over 0.5 seconds
    data.velo_smooth = data.velo;
    data.velo_smooth.data = smooth(data.velo.data, round(1/mode(diff(data.pos.tvec)))*.5)'; % smooth with gaussian



    % get the mean velocity when the animal is moving.
    mVelo(iS) = mean(data.velo_smooth.data(data.velo_smooth.data>thresh));

%     % get the mean/max velocity per trial type. 
%     trial_mean_velo = NaN(length(data.app.t), 4); % fill with NaNs. 
%     trial_max_velo = NaN(length(data.app.t), 4); 
% 
%     for iT = length(data.app.t):-1:1
%         if (data.app.t(iT) < data.velo_smooth.tvec(1)) || (data.rew.t(iT) > data.velo_smooth.tvec(end))
%             continue
%         end
%         A_velo = restrict(data.velo_smooth, data.app.t(iT), data.rew.t(iT)); 
%         fprintf('Approach dur: %0.2f  velo: %0.2f \n', A_velo.tvec(end) - A_velo.tvec(1), mean(A_velo.data))
% 
%         if A_velo.tvec(end) - A_velo.tvec(1) <= 10 % skip indirect trials that took more than 10s. 
% 
% 
%             trial_mean_velo(iT,data.app.in(iT)) = mean(A_velo.data); 
%             trial_max_velo(iT,data.app.in(iT)) = max(A_velo.data);
%         end
% 
%     end
% 
% all_velo_mean(iS,:) = nanmean(trial_mean_velo); 
% all_velo_max(iS,:) = nanmean(trial_mean_velo); 


end
%%
S = [];
S.C_idx = contains(s_phase, 'C1') | contains(s_phase, 'C2') | contains(s_phase, 'C3');
S.OE_idx = contains(s_phase, 'O1') | contains(s_phase, 'O2') | contains(s_phase, 'O3');
S.OL_idx = contains(s_phase, 'O4') | contains(s_phase, 'O5') | contains(s_phase, 'O6') | contains(s_phase, 'O7');
S.R_idx = contains(s_phase, 'R');

C1 = mVelo(contains(s_phase, 'C1'));
C2 = mVelo(contains(s_phase, 'C2'));
C3 = mVelo(contains(s_phase, 'C3'));
O1 = mVelo(contains(s_phase, 'O1'));
O2 = mVelo(contains(s_phase, 'O2'));
O3 = mVelo(contains(s_phase, 'O3'));
O4 = mVelo(contains(s_phase, 'O4'));
O5 = mVelo(contains(s_phase, 'O5'));
O6 = mVelo(contains(s_phase, 'O6'));
O7 = mVelo(contains(s_phase, 'O7'));
R1 = mVelo(contains(s_phase, 'R1'));
R2 = mVelo(contains(s_phase, 'R2'));
R3 = mVelo(contains(s_phase, 'R3'));



% export as csv
mat_out = NaN(13, max([length(C1),length(C2),length(C3),...
    length(O1), length(O2),length(O3),length(O4),length(O5),length(O6),length(O7),...
    length(R1),length(R2),length(R3)]));


mat_out(1,1:length(C1)) = C1';
mat_out(2,1:length(C2)) = C2';
mat_out(3,1:length(C3)) = C3';
mat_out(4,1:length(O1)) = O1';
mat_out(5,1:length(O2)) = O2';
mat_out(6,1:length(O3)) = O3';
mat_out(7,1:length(O4)) = O4';
mat_out(8,1:length(O5)) = O5';
mat_out(9,1:length(O6)) = O6';
mat_out(10,1:length(O7)) = O7';
mat_out(11,1:length(R1)) = R1';
mat_out(12,1:length(R2)) = R2';
mat_out(13,1:length(R3)) = R3';


mat_out= mat_out';

csvwrite([parent_path filesep 'Mean_Velo_min' strrep(num2str(thresh), '.', 'p') '.csv'], mat_out)

%% check the mean response for modulated cells

C_idx = find(contains(s_P_id, 'C'));
% C_idx = find(contains(s_P_id, 'O1') | contains(s_P_id, 'O2') | contains(s_P_id, 'O3'));


figure(666)
cla
hold on
for ii = length(C_idx):-1:1
    if s_H(C_idx(ii), 5) && (s_bFR(ii,5) < s_rFR(ii,5))
        plot(tvec, z_peth_s(:,5,C_idx(ii)), 'color', c_ord(1,:))
    end

end


OL6_idx = find(contains(s_P_id, 'O5') | contains(s_P_id, 'O6') | contains(s_P_id, 'O7'));


for ii = length(OL6_idx):-1:1
    if s_H(OL6_idx(ii), 5) && (s_bFR(ii,5) < s_rFR(ii,5))
        plot(tvec, z_peth_s(:,5,OL6_idx(ii)), 'color', c_ord(3,:))
    end

end



%%  example PETHS
% 27 36 111 144 72 44 23inter
% 27 36 111 144 72
ex_cells ={s_cell_id{[23 27 36 44 111 144 72]}};%{'C5_2_C1_2021-04-20_DONE_maze_data_TT2_02.t64'};
c_ord = linspecer(5); 
  
for iC = 1:length(ex_cells)
close all

    c_idx = find(contains(s_cell_id, ex_cells{iC})); 

    parts = strsplit(ex_cells{iC}, '_TT');
    load([parts{1} '.mat'])

    % make a color map
    for ii = length(data.rew.in):-1:1
        Zone_cord(ii,:) = c_ord(data.rew.in(ii),:);
    end


    this_S =  KA_isolate_S(data.S,['TT' parts{2}]);

    cfg_peth = [];
    cfg_peth.window = [-6 6];
    cfg_peth. plot_type = 'raw';
    cfg_peth.dt = 0.05;
    cfg_peth.gauss_sd = .2;
    cfg_peth.waves = this_S.waves{1};
    cfg_peth.evt_color_mat = Zone_cord; 

    example_peth = [];  

    for jj = unique(data.rew.in)
        this_idx = data.rew.in == jj;
        [~,outputIT{jj},this_peth] = SpikePETH_Shuff(cfg_peth, this_S, data.rew.t(this_idx) );
        example_peth(:,jj) = nanmean(this_peth,2);
    end

         cfg_peth.evt_color_mat = Zone_cord;
         cfg_peth.markersize = 15; 
        [~,outputIT{5},this_peth] = SpikePETH_Shuff(cfg_peth, this_S, data.rew.t);
        example_peth(:,5) = nanmean(this_peth,2);

    % collect data from PETHs

    subplot(211)
    title([strrep(parts{1}(1:18),'_', '-') ' TT' strrep(parts{2}(1:end-4), '_', '-')])
        xlim([-5 5])

    subplot(212)
    cla
    hold on
    for ii =1: size(example_peth,2)
        if s_H(c_idx, ii)
            plot(outputIT{5}, example_peth(:,ii),'-',  'color', [c_ord(ii,:), 1], 'linewidth', 3)
        else
            plot(outputIT{5}, example_peth(:,ii),'--',  'color', [c_ord(ii,:), .8], 'linewidth', 2)
        end
    end

    y_lim = [min(example_peth, [], 'all') max(example_peth, [], 'all')];
    ylim(y_lim);
        leg = legend([ fliplr(f_id(1:4)), 'All', ], 'Orientation', 'horizontal', 'Location', 'northwest','FontSize',14);
    set(leg, 'box', 'off')

    % move the pre post means up
    chil = get(gca, 'Children');

    % make the text finder adaptive.
    for ii = length(chil):-1:1
        type{ii} = chil(ii).Type;
    end
    % text_idx = find(contains(type, 'text'));
    % rec_idx = find(contains(type, 'rectangle'));
    % 
    % chil(text_idx(1)).Position = [chil(text_idx(1)).Position(1) max(example_peth, [], 'all')*.3 chil(text_idx(1)).Position(3)];
    % chil(text_idx(2)).Position = [chil(text_idx(2)).Position(1) max(example_peth, [], 'all')*.5 chil(text_idx(2)).Position(3)];
    % % adjust the vertical line at 0 which is a rectanlge.
    % chil(rec_idx).Position = [chil(rec_idx).Position(1) y_lim(1) chil(rec_idx).Position(3) y_lim(2) - y_lim(1)];

    % add velocity
    velo_window = [cfg_peth.window(1)*floor(1/mode(diff(data.velo_smooth.tvec))), cfg_peth.window(2)*floor(1/mode(diff(data.velo_smooth.tvec)))];
    all_velo = NaN(length(data.rew.t), (abs(velo_window(1)) + abs(velo_window(2)) +1));
    for ii = length(data.rew.t):-1:1 
        this_idx = nearest_idx3(data.rew.t(ii), data.velo_smooth.tvec);

        if this_idx < abs(velo_window(1)) || velo_window(2)+this_idx > length(data.velo_smooth.data)
            continue
        end
        all_velo(ii,:) = data.velo_smooth.data((velo_window(1)+this_idx):(velo_window(2)+this_idx));
    end

    velo_mean = nanmedian(all_velo, 1);
    velo_SEM = nanstd(all_velo)./sqrt(length(all_velo));
    velo_tvec = cfg_peth.window(1) : 1/floor(1/mode(diff(data.velo_smooth.tvec))):cfg_peth.window(2);

    yyaxis right
    hv = shadedErrorBar(velo_tvec,velo_mean,nanstd(all_velo) / sqrt(size(all_velo,1)));
    hv.mainLine.Color =  [.5 .5 .5 .3];
    hv.edge(1).Color =  [.5 .5 .5 .3];
    hv.edge(2).Color =  [.5 .5 .5 .3];

    hv.patch.FaceAlpha =  .2;
    ax = gca;
    ax.YAxis(1).Color = [.5 .5 .5];
    ax.YAxis(2).Color = [.5 .5 .5];
    ylabel('speed (cm/s)')


    subplot(212)
    y_lim = ylim;
    % rectangle('position', [0 y_lim(1) 0.001  y_lim(2) - y_lim(1)], 'facecolor', [[4,172,218]./255 0.5], 'edgecolor', [[4,172,218]./255 0.5])

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = [.5 .5 .5];

    leg = legend([f_id(1:4),'All' , 'speed'], 'Orientation', 'horizontal', 'Location', 'northwest','FontSize',14);
    set(leg, 'box', 'off')

    SetFigure([], gcf);
    xlim([-5 5])
% maximize
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
pause(2)

print(gcf,[parent_path filesep strrep(parts{1}(1:18), '-', '_') ' TT' strrep(parts{2}(1:end-4), '-', '_') '_PETH.pdf'],'-dpdf','-r300')


end
%% collect the position data for E and R sessions. 
E.map = []; R1.map = []; R2.map = []; R3.map = []; 
E.dur = []; R1.dur = []; R2.dur = []; R3.dur = []; 
EG.map = []; R1G.map = []; R2G.map = []; R3G.map = []; 

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

keep_idx = (contains(sess_list, 'E1') | contains(sess_list, 'R')); 

sess_list(~keep_idx) = [];

for iS = length(sess_list):-1:1
    cd([data_dir filesep sess_list{iS}])

    disp([ num2str(iS) '_' sess_list{iS}])
    s_phase{iS}= sess_list{iS}(6:7);

    if ~isempty(dir('*VT*.zip')) && isempty(dir('*.nvt'))
        unzip('VT1.zip')
    end

    evt = LoadEvents([]);
    % add in check for multiple recording periods.  Some seem to have a pre and post recoding.
    s_rec_idx = find(contains(evt.label, 'Starting Record'));
    e_rec_idx = find(contains(evt.label, 'Stopping Record'));
    if length(evt.t{s_rec_idx}) > length(evt.t{e_rec_idx})
        evt.t{s_rec_idx}  = evt.t{s_rec_idx}(1:length(evt.t{e_rec_idx})); 

    end
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
    % if max(rec_dur)/60 < 10
    %     fprintf('<strong>Minumum recording duration (10mins) not met: %2.1fmins</strong>\n', max(rec_dur)/60);
    %     continue
    % end


    cfg_pos.convFact = [560/142 480/142];
    data.pos = LoadPos(cfg_pos);
    data.pos = restrict(data.pos, evt.t{s_rec_idx}(task_rec_idx_s), evt.t{e_rec_idx}(task_rec_idx_e)); % restrict position to time on track.

% occupancy histogram. 
    SET_xmin = 0; SET_ymin = 0; % set up bins
    SET_xmax = 180; SET_ymax = 180;
    SET_xBinSz = (SET_xmax - SET_xmin)/30; SET_yBinSz = (SET_xmax - SET_xmin)/30;


    x_edges = SET_xmin:SET_xBinSz:SET_xmax;
    y_edges = SET_ymin:SET_yBinSz:SET_ymax;

    % set up gaussian
    kernel = gausskernel([8 8],4); % 2d gaussian in bins

    % compute occupancy
    occ_hist = hist3(data.pos.data(1:2,:)', 'edges', {y_edges x_edges});
    %     occ_hist = histcn(out.pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
    low_occ_idx = find(occ_hist <7.5 );
    occ_hist(low_occ_idx) = 0;

    %     occ_hist = conv2(occ_hist,kernel,'same');

    no_occ_idx = find(occ_hist == 0 ); % NaN out bins never visited
    occ_hist(no_occ_idx) = NaN;
    %
    occ_hist = occ_hist .* mode(diff(data.pos.tvec)); % convert samples to seconds using video frame rate (30 Hz)


    if contains(sess_list{iS}, 'C1_1') || contains(sess_list{iS}, 'C4_3') || contains(sess_list{iS}, 'C5_3') || contains(sess_list{iS}, 'C6_4')
        if contains(sess_list{iS}, 'E1')
            E.map = cat(3,E.map, occ_hist);
            E.dur = [E.dur rec_dur];
        elseif contains(sess_list{iS}, 'R1')
            R1.map = cat(3,R1.map, occ_hist);
        elseif contains(sess_list{iS}, 'R2')
            R2.map = cat(3,R2.map, occ_hist);
        elseif contains(sess_list{iS}, 'R3')
            R3.map = cat(3,R3.map, occ_hist);
        end

    elseif contains(sess_list{iS}, 'C3_2') || contains(sess_list{iS}, 'C3_3') || contains(sess_list{iS}, 'C5_2') || contains(sess_list{iS}, 'C6_2')

        if contains(sess_list{iS}, 'E1')
            EG.map = cat(3,EG.map, occ_hist);
            E.dur = [E.dur rec_dur];
        elseif contains(sess_list{iS}, 'R1')
            R1G.map = cat(3,R1G.map, occ_hist);
        elseif contains(sess_list{iS}, 'R2')
            R2G.map = cat(3,R2G.map, occ_hist);
        elseif contains(sess_list{iS}, 'R3')
            R3G.map = cat(3,R3G.map, occ_hist);
        end
    else

        disp(sess_list{iS})
    end




end

figure(880)
cla
subplot(2,4,1)
imagesc(nanmean(E.map, 3))
title('E Orange (E/W)')

subplot(2,4,2)
imagesc(nanmean(R1.map, 3))
title('R1')

subplot(2,4,3)
imagesc(nanmean(R2.map, 3))
title('R2')

subplot(2,4,4)
imagesc(nanmean(R3.map, 3))
title('R3')


subplot(2,4,5)
imagesc(nanmean(EG.map, 3))
title('E grape (N/S)')

subplot(2,4,6)
imagesc(nanmean(R1G.map, 3))
title('R1')

subplot(2,4,7)
imagesc(nanmean(R2G.map, 3))
title('R2')

subplot(2,4,8)
imagesc(nanmean(R3G.map, 3))
title('R3')

cmap = [1 1 1; parula(256)]; 

colormap(cmap)


%% count the feeder fires 

for iS = length(sess_list):-1:1
    cd([data_dir filesep sess_list{iS}])




        this_data = KA_trialfun_noS(plot_dir);

        for ii = unique(this_data.rew.in)
            this_idx = this_data.rew.in == ii; 
            this_ratio(ii) = sum(this_idx) / sum(this_data.PM.FeedersFired == ii); 
        end



end


%% speed plots

figure(868)
clf
subplot(6,3,[1 4 7])
bins = -.6:0.01:.6;
hold on
histogram(spd_corr(~logical(spd_mod)), bins, 'FaceColor',[.7 .7 .7])
histogram(spd_corr(logical(spd_mod)), bins, 'FaceColor',c_ord(2,:))
ylabel('count')
xlabel('speed corr')

subplot(6,3,[10 13 16])

hold on
scatter(spd_corr(logical(spd_mod)), fr(logical(spd_mod)),  35, 'markerfacecolor', c_ord(2,:), 'markeredgecolor', c_ord(2,:))
scatter(spd_corr(~logical(spd_mod)), fr(~logical(spd_mod)),  35, 'markerfacecolor', [.7 .7 .7], 'markeredgecolor', [.7 .7 .7])
xlabel('speed corr')

% bins = -10:0.25:10; 
% hold on; 
% histogram(spd_z(abs(spd_z) <2.58) , bins, 'FaceColor',[.7 .7 .7])
% histogram(spd_z(abs(spd_z) >=2.58) , bins, 'FaceColor',c_ord(2,:))
% ylabel('count')
% xlabel('z score corr')


% rank the spd corr; 
[~, spd_mod_max_idx] = sort(spd_corr, 'descend'); 
[~, spd_mod_min_idx] = sort(abs(spd_corr), 'descend'); 

ex_idx = [spd_mod_max_idx(1) spd_mod_max_idx(3) spd_mod_max_idx(end-3) spd_mod_max_idx(end) spd_mod_min_idx(end-3) spd_mod_min_idx(end) ]; 
% spd_c_ord  = [0.2238    0.4408    0.7226;0.2069    0.5754    0.7336; ...
%     0.7598    0.1518    0.3010; 0.9745    0.5140    0.2853;...
%      .6769    0.4447    0.7114 ; .6769    0.4447    0.7114];

spd_c_ord = [67, 127, 151; 132 147 35; 255 179 13; 253 22 26; 80 80 80; 80 80 80]/255;
splt_idx = 2:3:17; 

for ii =1:length(ex_idx)

subplot(6,3,[1 4 7])
scatter(spd_corr(ex_idx(ii)), 5, 35, 'markerfacecolor', spd_c_ord(ii,:), 'markeredgecolor', spd_c_ord(ii,:))

subplot(6,3,[10 13 16])
scatter(spd_corr(ex_idx(ii)), fr(ex_idx(ii)), 55, 'd',  'markerfacecolor', spd_c_ord(ii,:), 'markeredgecolor', spd_c_ord(ii,:))



subplot(6, 3, [splt_idx(ii) splt_idx(ii)+1])
yyaxis right
plot(spd_data{ex_idx(ii)}.tvec - spd_data{ex_idx(ii)}.tvec(1), spd_data{ex_idx(ii)}.spd, 'LineWidth',1,'color',  [.5 .5 .5])
ylabel('speed (cm/s)')

yyaxis left

plot(spd_data{ex_idx(ii)}.tvec - spd_data{ex_idx(ii)}.tvec(1), spd_data{ex_idx(ii)}.FR, 'LineWidth',2, 'color', spd_c_ord(ii,:))
ylabel('firing rate (Hz)')

xlim([60 180])

ax = gca; 
ax.YAxis(1).Color = spd_c_ord(ii,:); 
ax.YAxis(2).Color = [.5 .5 .5];

if ii ~= length(ex_idx)
    set(gca, 'xtick', [])
else
    set(gca, 'xtick', [0:60:180])
    xlabel('times (sec)')
end

title([ 'Speed Corr: ' num2str(spd_corr(ex_idx(ii)),'%4.2f') ' | ' strrep([cell_id{ex_idx(ii)}(1:7) ' | ' cell_id{ex_idx(ii)}(end-9:end-4)], '_', ' ') ])
end


subplot(6,3,[10 13 16])
set(gca, 'YScale', 'log')
ylabel('log firing rate (Hz)')

%% try some linear decoding using the Lopes et al. 2025 PlosOne method with a 3rd order polynomial fit. 
%load([parent_dir filesep 'Spd_data.mat']); 

z_err = NaN(length(spd_data),1);
R2 = z_err; 
S_R2 = cell(size(R2)); 

parfor kk = 1:length(spd_data)
   disp(kk)

    [z_err(kk), R2(kk), S_R2{kk}] = KA_lin_decode([], spd_data{kk}.FR, spd_data{kk}.spd); 

end

    [~, ~, ~, plt_mat, plt_s_mat] = KA_lin_decode([], spd_data{27}.FR, spd_data{27}.spd); 

    nan_idx = isnan(spd_data{27}.FR) | isnan(spd_data{27}.spd); % remove nans to match plt_mat

%% reproduce the Lopes et al. Fig 7. 
figure(1011)
clf
subplot(2,3,1:2); cla
hold on
plot(spd_data{27}.tvec(~nan_idx) - spd_data{27}.tvec(1), spd_data{27}.spd(~nan_idx), 'k', 'linewidth', 3); 
plot(spd_data{27}.tvec(~nan_idx)- spd_data{27}.tvec(1), mean(plt_mat, 2, 'omitnan'),'-', 'color', c_ord(2,:), 'linewidth', 3); 
plot(spd_data{27}.tvec(~nan_idx)- spd_data{27}.tvec(1), mean(plt_s_mat, 2, 'omitnan'), 'color', [.7 .7 .7], 'linewidth', 3); 
xlim([260 360])
ylim([0 inf])
legend('Actual', 'Decoded', 'Shuffle', 'box',  'off')


subplot(2,3,3)
MS_bar_w_err(R2(~FS_idx' & logical(spd_mod),1)', R2(FS_idx' & logical(spd_mod),1)', [c_ord(1,:); c_ord(2,:)] , 1, 'ttest2', 1:2);
set(gca, 'xticklabel', {'Pyr Speed Cells' 'FS Speed Cells'}, 'XTickLabelRotation', 45)
ylabel('Decoding accuracy R^2')
set(gca, 'YScale', 'log')

% axis square


subplot(2,3,4)
MS_bar_w_err(R2(spd_z > 1.98,1)', R2(spd_z<-1.98,1)', [c_ord(4,:)*1.5; c_ord(5,:)*2] , 1, 'ttest2', 1:2);
set(gca, 'xticklabel', {'Pos Speed Cells' 'Neg Speed Cells'}, 'XTickLabelRotation', 45)
ylabel('Decoding accuracy R^2')

set(gca, 'YScale', 'log')
% axis square


subplot(2,3,5:6)
MS_bar_w_err3(R2(logical(spd_mod),1)', R2(~logical(spd_mod),1)',S_R2(:,1)', [c_ord(3,:); .7 .7 .7; .3 .3 .3], 1, 'anova1', 1:3);
set(gca, 'xticklabel', {'Speed Cells' 'Non-Speed Cells' 'shuffle'}, 'XTickLabelRotation', 45)
set(gca, 'YScale', 'log')
% axis square

% MS_bar_w_err(R2(:,1)', S_R2(:,1)', [.25 .25 .25 ; .7 .7 .7] , 1, 'ttest', 1:2);
% set(gca, 'xticklabel', {'All Cells' 'Shuffle'})
% axis square
% ylabel('Decoding accuracy R^2')

% subplot(2,3,4)
% MS_bar_w_err(R2(:,1)', S_R2(:,1)', [.25 .25 .25 ; .7 .7 .7] , 1, 'ttest', 1:2);
% set(gca, 'xticklabel', {'All Cells' 'Shuffle'})
% set(gca, 'YScale', 'log')
% axis square
% ylabel('Decoding accuracy log R^2')



% subplot(2,3,3)
% MS_bar_w_err(R2(logical(spd_mod),1)', R2(~logical(spd_mod),1)', [c_ord(1,:); c_ord(1,:)*2] , 1, 'ttest2', 1:2);
% set(gca, 'xticklabel', {'Speed Cells' 'Other'})
% set(gca, 'YScale', 'log')
% axis square




% subplot(2,3,6)
% MS_bar_w_err(R2(~FS_idx' & logical(spd_mod),1)', R2(FS_idx' & logical(spd_mod),1)', [c_ord(1,:); c_ord(2,:)] , 1, 'ttest2', 1:2);
% set(gca, 'xticklabel', {'Pyr Speed Cells' 'FS Speed Cells'})
% set(gca, 'YScale', 'log')
% axis square

% maximize
set(gcf,'Units','Inches');
pos = get(gcf,'Position');

set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
pause(2)

print(gcf,[parent_path filesep 'lin_dec_R2_log.pdf'],'-dpdf','-r300', '-bestfit')

saveas(gcf,[parent_path filesep 'R2.fig']);
% print(gcf,[parent_path filesep 'R2.pdf'],'-dpdf','-r300')

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

%% get the trial velocity per condition

cd(inter_dir)
sess_list = dir([inter_dir filesep '*.mat']);

s_phase= []; mVelo = []; all_velo_mean = NaN(length(sess_list), 5); all_velo_max = NaN(length(sess_list), 5); 


for iS = length(sess_list):-1:1


    load([inter_dir filesep sess_list(iS).name])


    % get the mean/max velocity per trial type.
    trial_mean_velo = NaN(length(data.app.t), 4); % fill with NaNs.
    trial_max_velo = NaN(length(data.app.t), 4);

    for iT = length(data.app.t):-1:1
        if (data.app.t(iT) < data.velo_smooth.tvec(1)) || (data.rew.t(iT) > data.velo_smooth.tvec(end)) || ((data.rew.t(iT) - data.app.t(iT))<0)
            continue
        end
        A_velo = restrict(data.velo_smooth, data.app.t(iT)-2.5, data.rew.t(iT));
        % fprintf('Approach dur: %0.2f  velo: %0.2f \n', A_velo.tvec(end) - A_velo.tvec(1), mean(A_velo.data))

        if A_velo.tvec(end) - A_velo.tvec(1) <= 10 % skip indirect trials that took more than 10s.

            trial_mean_velo(iT,data.app.in(iT)) = mean(A_velo.data);
            trial_max_velo(iT,data.app.in(iT)) = max(A_velo.data);
        end

    end

    all_velo_mean(iS,1:4) = nanmean(trial_mean_velo);
    all_velo_mean(iS,5)  = nanmean(trial_mean_velo, 'all'); 

    all_velo_max(iS,1:4) = nanmean(trial_max_velo); 
    all_velo_max(iS,5)  = nanmean(trial_max_velo, 'all');

    all_velo(iS) = mean(data.velo_smooth.data); 
    all_mvelo(iS) = mean(data.velo_smooth.data(data.velo_smooth.data>2.5)); 

    s_ID{iS} = sess_list(iS).name(6:7); 
    c_ID{iS} = sess_list(iS).name(1:4); 
    F_id{iS} = data.PM.Feeder_type; 
    F_val(iS,1:4) = data.PM.Feeder_mag; 

end

% export as csv
mat_out = [];
mat_out = array2table([all_velo_mean';  all_velo_max'; all_velo; all_mvelo]', 'VariableNames',{'Mean_B3', 'Mean_G3', 'Mean_B1', 'Mean_G1', 'Mean_all',...
    'Max_B3', 'Max_G3', 'Max_B1', 'Max_G1', 'Max_all', 'overall_mean', 'overall_mean_mov_only'}); 

c_type = []; 

mat_out.sess_id = s_ID'; 
mat_out.sub_id = c_ID'; 

% 
writetable(mat_out, [parent_path filesep 'Sess_spd_2p5.csv'])


%% try some information theory modeling
% data
% x = spd_data{27}.FR; 
% y = spd_data{27}.spd; 
% 
% k = glmfit(x, y, 'Poisson')
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
%
% %% old code
% set(gca, 'ytick', [])
% title({'East'; 'grain x 1'})
% % set(gca, 'ytick', 1:size(mean_E_mat,1), 'yTickLabel', sess_label)
% set(gca, 'ytick', [])
% set(gca, 'xtick', -5:2.5:5)
% vl = vline(0, '--k');
% vl.LineWidth = 2;
% caxis([Z_sig_min Z_sig_max]);
% hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
% for ii = 1:length(hl)
%     hl(ii).LineWidth = 3;
%     hl(ii).Color = c_ord(ii+1,:);
% end
% for ii = 1:length(all_mat_order)
%     if east_mat_sig_pre(ii)
%         text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
%     elseif east_mat_sig_post(ii)
%         text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
%     end
%     if E_deval_mat(ii)
%         text(max(tvec.C1)+mode(diff(tvec.C1))*6, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
%     end
% end

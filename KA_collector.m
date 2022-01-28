function [all_out] = KA_collector(inter_dir, smooth_val)

%% KA_collector
%
%   Collects the output of KA_master + KA_screener files stored in inter_dir.  Compiles the cells and
%   splits them into each session type.
%
%  Session type code:
%       Acquisition, Criteria, Overtrain, Extinction, Reacquisition.
%
%  nSessions: A ? C 3, O 7, E 1, R 3-5
%

%% initialize
if nargin == 0
    inter_dir = cd;
end


% custom colormap
purple = [18 13 44]./255;
yellow = [246 207 70]./255;
cmapsize = 256;   %or as passed in
NumPoints = floor(cmapsize * 15/100) - floor(cmapsize * 0/100);
SUNS_map = [linspace(purple(1) , yellow(1), NumPoints).', ...
    linspace(purple(2) , yellow(2), NumPoints).', ...
    linspace(purple(3) , yellow(3), NumPoints).'];
% sunsmap
sunsmap = [0.2422    0.1504    0.6603
    0.2504    0.1650    0.7076
    0.2578    0.1818    0.7511
    0.2647    0.1978    0.7952
    0.2706    0.2147    0.8364
    0.2751    0.2342    0.8710
    0.2783    0.2559    0.8991
    0.2803    0.2782    0.9221
    0.2813    0.3006    0.9414
    0.2810    0.3228    0.9579
    0.2795    0.3447    0.9717
    0.2760    0.3667    0.9829
    0.2699    0.3892    0.9906
    0.2602    0.4123    0.9952
    0.2440    0.4358    0.9988
    0.2206    0.4603    0.9973
    0.1963    0.4847    0.9892
    0.1834    0.5074    0.9798
    0.1786    0.5289    0.9682
    0.1764    0.5499    0.9520
    0.1687    0.5703    0.9359
    0.1540    0.5902    0.9218
    0.1460    0.6091    0.9079
    0.1380    0.6276    0.8973
    0.1248    0.6459    0.8883
    0.1113    0.6635    0.8763
    0.0952    0.6798    0.8598
    0.0689    0.6948    0.8394
    0.0297    0.7082    0.8163
    0.0036    0.7203    0.7917
    0.0067    0.7312    0.7660
    0.0433    0.7411    0.7394
    0.0964    0.7500    0.7120
    0.1408    0.7584    0.6842
    0.1717    0.7670    0.6554
    0.1938    0.7758    0.6251
    0.2161    0.7843    0.5923
    0.2470    0.7918    0.5567
    0.2906    0.7973    0.5188
    0.3406    0.8008    0.4789
    0.3909    0.8029    0.4354
    0.4456    0.8024    0.3909
    0.5044    0.7993    0.3480
    0.5616    0.7942    0.3045
    0.6174    0.7876    0.2612
    0.6720    0.7793    0.2227
    0.7242    0.7698    0.1910
    0.7738    0.7598    0.1646
    0.8203    0.7498    0.1535
    0.8634    0.7406    0.1596
    0.9035    0.7330    0.1774
    0.9393    0.7288    0.2100
    0.9728    0.7298    0.2394
    0.9956    0.7434    0.2371
    0.9970    0.7659    0.2199
    0.9952    0.7893    0.2028
    0.9892    0.8136    0.1885
    0.9786    0.8386    0.1766
    0.9676    0.8639    0.1643
    0.9610    0.8890    0.1537
    0.9597    0.9135    0.1423
    0.9628    0.9373    0.1265
    0.9691    0.9606    0.1064
    0.9769    0.9839    0.0805];

%% allocate some empty arrays to fill in.
% for ii = 1:15 % Acqusition stage could talke many days. Will trim empty ones later.
%     all_out.(['A' num2str(ii)]) = [];
%     all_sig.(['A' num2str(ii)]) = [];
%
% end
subj_list = {'C1_1', 'C3_2','C3_3', 'C3_4', 'C4_3', 'C5_2', 'C5_3','C6_3', 'C6_4'};
deval_type = {'WE', 'NS', 'NS', ' ' , 'WE', 'WE', 'NS', 'WE', 'NS'};
for ii = 1:3  % Criteria should be 3
    all_out.(['C' num2str(ii)]) = [];
    all_sig.(['C' num2str(ii)]) = [];
end

for ii = 1:7 % Overtaining should b 7
    all_out.(['O' num2str(ii)]) = [];
    all_sig.(['O' num2str(ii)]) = [];
end

for ii = 1:1 % Extinction is only 1
    all_out.(['E' num2str(ii)]) = [];
    all_sig.(['E' num2str(ii)]) = [];
end

for ii = 1:5 % Reacqusition is up to 5
    all_out.(['R' num2str(ii)]) = [];
    all_sig.(['R' num2str(ii)]) = [];
end

% clone for each feeder type
north_out = all_out;
east_out = all_out;
south_out = all_out;
west_out = all_out;

north_sig = all_sig;
east_sig = all_sig;
south_sig = all_sig;
west_sig = all_sig;

all_sig_pre = all_sig;
north_sig_pre = all_sig;
east_sig_pre = all_sig;
south_sig_pre = all_sig;
west_sig_pre = all_sig;

all_sig_post = all_sig;
north_sig_post = all_sig;
east_sig_post = all_sig;
south_sig_post = all_sig;
west_sig_post = all_sig;

all_sig_act = all_sig;
north_sig_act = all_sig;
east_sig_act = all_sig;
south_sig_act = all_sig;
west_sig_act = all_sig;

all_sig_act_pre = all_sig;
north_sig_act_pre = all_sig;
east_sig_act_pre = all_sig;
south_sig_act_pre = all_sig;
west_sig_act_pre = all_sig;

all_sig_act_post = all_sig;
north_sig_act_post = all_sig;
east_sig_act_post = all_sig;
south_sig_act_post = all_sig;
west_sig_act_post = all_sig;

all_sig_rew = all_sig;
north_sig_rew = all_sig;
east_sig_rew = all_sig;
south_sig_rew = all_sig;
west_sig_rew = all_sig;

all_sig_rew_pre = all_sig;
north_sig_rew_pre = all_sig;
east_sig_rew_pre = all_sig;
south_sig_rew_pre = all_sig;
west_sig_rew_pre = all_sig;

all_sig_rew_post = all_sig;
north_sig_rew_post = all_sig;
east_sig_rew_post = all_sig;
south_sig_rew_post = all_sig;
west_sig_rew_post = all_sig;


all_out_perc  = all_out;

all_subj = all_out;

all_FR = all_out;

% cycle and collect

cd(inter_dir);
file_names = FindFiles('*Feeder*');  % get all the files containing 'TT' which should be all the outputs from KA_master.


for iF = 1:length(file_names)
    load(file_names{iF}); % load the data as 'This_cell'
    
    % put additional criteria for inclusion here.
    
    %     if str2double(file_names{iF}(strfind(file_names{iF}, '.mat')-1)) > 3
    %         fprintf('Quality is %d, skipping %s\n', str2double(file_names{iF}(strfind(file_names{iF}, '.mat')-1)) , file_names{iF})
    %         continue
    %     end
    %
    
    %     velo_out
    
    this_sess = This_cell.session; % get the session type.
    all_subj.(this_sess) =[all_subj.(this_sess), find(ismember(subj_list, This_cell.subject))];
    
    all_subj_list(iF) = find(ismember(subj_list, This_cell.subject)); % index value for the subject.
    all_sess_names{iF} = [This_cell.subject '_' This_cell.session '_' This_cell.date];
    
    
    all_FR.(this_sess) = [all_FR.(this_sess) length(This_cell.S.t{1})/(This_cell.pos.tvec(end) - This_cell.pos.tvec(1))];
    
    if exist('smooth_val', 'var')
        all_out.(this_sess) = [all_out.(this_sess) filter(gausswin(smooth_val), 1,This_cell.mean_S_gau{5})];
        north_out.(this_sess) = [north_out.(this_sess) filter(gausswin(smooth_val), 1,This_cell.mean_S_gau{1})];
        west_out.(this_sess) = [west_out.(this_sess) filter(gausswin(smooth_val), 1,This_cell.mean_S_gau{2})];
        south_out.(this_sess) = [south_out.(this_sess) filter(gausswin(smooth_val), 1,This_cell.mean_S_gau{3})];
        east_out.(this_sess) = [east_out.(this_sess) filter(gausswin(smooth_val), 1,This_cell.mean_S_gau{4})];
        
    else
        all_out.(this_sess) = [all_out.(this_sess) This_cell.mean_S_gau{5}];
        north_out.(this_sess) = [north_out.(this_sess) This_cell.mean_S_gau{1}];
        west_out.(this_sess) = [west_out.(this_sess) This_cell.mean_S_gau{2}];
        south_out.(this_sess) = [south_out.(this_sess) This_cell.mean_S_gau{3}];
        east_out.(this_sess) = [east_out.(this_sess) This_cell.mean_S_gau{4}];
    end
    
    if isfield(This_cell, 'H')
        all_sig.(this_sess) =  [all_sig.(this_sess) This_cell.H{5}];
        
        % get sig for two sided ttest (any modulation)
        north_sig.(this_sess) =  [north_sig.(this_sess) This_cell.H{1}];
        north_sig.(this_sess) = north_sig.(this_sess) ==1;
        west_sig.(this_sess) =  [west_sig.(this_sess) This_cell.H{2}];
        west_sig.(this_sess) = west_sig.(this_sess) ==1;
        south_sig.(this_sess) =  [south_sig.(this_sess) This_cell.H{3}];
        south_sig.(this_sess) = south_sig.(this_sess) == 1;
        east_sig.(this_sess) =  [east_sig.(this_sess) This_cell.H{4}];
        east_sig.(this_sess) = east_sig.(this_sess) == 1;
        
        % same but only post reward > pre
        all_sig_post.(this_sess) =  [all_sig_post.(this_sess) This_cell.H_post{5}];
        all_sig_post.(this_sess) = all_sig_post.(this_sess) ==1;
        north_sig_post.(this_sess) =  [north_sig_post.(this_sess) This_cell.H_post{1}];
        north_sig_post.(this_sess) = north_sig_post.(this_sess) ==1;
        west_sig_post.(this_sess) =  [west_sig_post.(this_sess) This_cell.H_post{2}];
        west_sig_post.(this_sess) = west_sig_post.(this_sess) ==1;
        south_sig_post.(this_sess) =  [south_sig_post.(this_sess) This_cell.H_post{3}];
        south_sig_post.(this_sess) = south_sig_post.(this_sess) == 1;
        east_sig_post.(this_sess) =  [east_sig_post.(this_sess) This_cell.H_post{4}];
        east_sig_post.(this_sess) = east_sig_post.(this_sess) == 1;
        
        % same but only pre reward > post
        all_sig_pre.(this_sess) =  [all_sig_pre.(this_sess) This_cell.H_pre{5}];
        all_sig_pre.(this_sess) = all_sig_pre.(this_sess) ==1;
        north_sig_pre.(this_sess) =  [north_sig_pre.(this_sess) This_cell.H_pre{1}];
        north_sig_pre.(this_sess) = north_sig_pre.(this_sess) ==1;
        west_sig_pre.(this_sess) =  [west_sig_pre.(this_sess) This_cell.H_pre{2}];
        west_sig_pre.(this_sess) = west_sig_pre.(this_sess) ==1;
        south_sig_pre.(this_sess) =  [south_sig_pre.(this_sess) This_cell.H_pre{3}];
        south_sig_pre.(this_sess) = south_sig_pre.(this_sess) == 1;
        east_sig_pre.(this_sess) =  [east_sig_pre.(this_sess) This_cell.H_pre{4}];
        east_sig_pre.(this_sess) = east_sig_pre.(this_sess) == 1;
    elseif isfield(This_cell, 'H_act_mod')
        
        all_sig_act.(this_sess) =  [all_sig_act.(this_sess) This_cell.H_act_mod{5}];
        all_sig_rew.(this_sess) =  [all_sig_rew.(this_sess) This_cell.H_rew_mod{5}];
        all_sig.(this_sess) = [all_sig.(this_sess) (This_cell.H_rew_mod{5} || This_cell.H_act_mod{5})];
        
        north_sig_act.(this_sess) = [north_sig_act.(this_sess) This_cell.H_act_mod{1}];
        north_sig_act.(this_sess) = north_sig_act.(this_sess) ==1;
        west_sig_act.(this_sess) =  [west_sig_act.(this_sess) This_cell.H_act_mod{2}];
        west_sig_act.(this_sess) =  west_sig_act.(this_sess) ==1;
        south_sig_act.(this_sess) = [south_sig_act.(this_sess) This_cell.H_act_mod{3}];
        south_sig_act.(this_sess) = south_sig_act.(this_sess) == 1;
        east_sig_act.(this_sess) =  [east_sig_act.(this_sess) This_cell.H_act_mod{4}];
        east_sig_act.(this_sess) =  east_sig_act.(this_sess) == 1;
        
        % same but only post act > pre
        all_sig_act_post.(this_sess) =  [all_sig_act_post.(this_sess) This_cell.H_act{5}];
        all_sig_act_post.(this_sess) = all_sig_act_post.(this_sess) ==1;
        north_sig_act_post.(this_sess) =  [north_sig_act_post.(this_sess) This_cell.H_act{1}];
        north_sig_act_post.(this_sess) = north_sig_act_post.(this_sess) ==1;
        west_sig_act_post.(this_sess) =  [west_sig_act_post.(this_sess) This_cell.H_act{2}];
        west_sig_act_post.(this_sess) = west_sig_act_post.(this_sess) ==1;
        south_sig_act_post.(this_sess) =  [south_sig_act_post.(this_sess) This_cell.H_act{3}];
        south_sig_act_post.(this_sess) = south_sig_act_post.(this_sess) == 1;
        east_sig_act_post.(this_sess) =  [east_sig_act_post.(this_sess) This_cell.H_act{4}];
        east_sig_act_post.(this_sess) = east_sig_act_post.(this_sess) == 1;
        
        % same but only pre act < post
        all_sig_act_pre.(this_sess) =  [all_sig_act_pre.(this_sess) This_cell.H_pre_a{5}];
        all_sig_act_pre.(this_sess) = all_sig_act_pre.(this_sess) ==1;
        north_sig_act_pre.(this_sess) =  [north_sig_act_pre.(this_sess) This_cell.H_pre_a{1}];
        north_sig_act_pre.(this_sess) = north_sig_act_pre.(this_sess) ==1;
        west_sig_act_pre.(this_sess) =  [west_sig_act_pre.(this_sess) This_cell.H_pre_a{2}];
        west_sig_act_pre.(this_sess) = west_sig_act_pre.(this_sess) ==1;
        south_sig_act_pre.(this_sess) =  [south_sig_act_pre.(this_sess) This_cell.H_pre_a{3}];
        south_sig_act_pre.(this_sess) = south_sig_act_pre.(this_sess) == 1;
        east_sig_act_pre.(this_sess) =  [east_sig_act_pre.(this_sess) This_cell.H_pre_a{4}];
        east_sig_act_pre.(this_sess) = east_sig_act_pre.(this_sess) == 1;
        
        
        north_sig_rew.(this_sess) = [north_sig_rew.(this_sess) This_cell.H_rew_mod{1}];
        north_sig_rew.(this_sess) = north_sig_rew.(this_sess) ==1;
        west_sig_rew.(this_sess) =  [west_sig_rew.(this_sess) This_cell.H_rew_mod{2}];
        west_sig_rew.(this_sess) =  west_sig_rew.(this_sess) ==1;
        south_sig_rew.(this_sess) = [south_sig_rew.(this_sess) This_cell.H_rew_mod{3}];
        south_sig_rew.(this_sess) = south_sig_rew.(this_sess) == 1;
        east_sig_rew.(this_sess) =  [east_sig_rew.(this_sess) This_cell.H_rew_mod{4}];
        east_sig_rew.(this_sess) =  east_sig_rew.(this_sess) == 1;
        
        % same but only post reward > pre
        all_sig_rew_post.(this_sess) =  [all_sig_rew_post.(this_sess) This_cell.H_rew{5}];
        all_sig_rew_post.(this_sess) = all_sig_rew_post.(this_sess) ==1;
        north_sig_rew_post.(this_sess) =  [north_sig_rew_post.(this_sess) This_cell.H_rew{1}];
        north_sig_rew_post.(this_sess) = north_sig_rew_post.(this_sess) ==1;
        west_sig_rew_post.(this_sess) =  [west_sig_rew_post.(this_sess) This_cell.H_rew{2}];
        west_sig_rew_post.(this_sess) = west_sig_rew_post.(this_sess) ==1;
        south_sig_rew_post.(this_sess) =  [south_sig_rew_post.(this_sess) This_cell.H_rew{3}];
        south_sig_rew_post.(this_sess) = south_sig_rew_post.(this_sess) == 1;
        east_sig_rew_post.(this_sess) =  [east_sig_rew_post.(this_sess) This_cell.H_rew{4}];
        east_sig_rew_post.(this_sess) = east_sig_rew_post.(this_sess) == 1;
        
        % same but only pre reward < post
        all_sig_rew_pre.(this_sess) =  [all_sig_rew_pre.(this_sess) This_cell.H_pre_r{5}];
        all_sig_rew_pre.(this_sess) = all_sig_rew_pre.(this_sess) ==1;
        north_sig_rew_pre.(this_sess) =  [north_sig_rew_pre.(this_sess) This_cell.H_pre_r{1}];
        north_sig_rew_pre.(this_sess) = north_sig_rew_pre.(this_sess) ==1;
        west_sig_rew_pre.(this_sess) =  [west_sig_rew_pre.(this_sess) This_cell.H_pre_r{2}];
        west_sig_rew_pre.(this_sess) = west_sig_rew_pre.(this_sess) ==1;
        south_sig_rew_pre.(this_sess) =  [south_sig_rew_pre.(this_sess) This_cell.H_pre_r{3}];
        south_sig_rew_pre.(this_sess) = south_sig_rew_pre.(this_sess) == 1;
        east_sig_rew_pre.(this_sess) =  [east_sig_rew_pre.(this_sess) This_cell.H_pre_r{4}];
        east_sig_rew_pre.(this_sess) = east_sig_rew_pre.(this_sess) == 1;
    end
    
    
    
    % % change from pre to post reward
    all_out_perc.(this_sess) = [all_out_perc.(this_sess), This_cell.post_stim_means{5} / This_cell.pre_stim_means{5}];
    
    
    tvec.(this_sess) = This_cell.outputIT{1};
    clear This_cell
end % iF files

%% compile into sessions and get the mean Z score.
sessions = fieldnames(all_out);
% make some empty matricies to fill in with each cell.
all_mat = [];
north_mat = [];
west_mat = [];
south_mat = [];
east_mat = [];

% collect the H for each cell (for all_mat)
all_mat_sig = [];  all_mat_sig_pre = [];  all_mat_sig_post = [];
north_mat_sig = [];  north_mat_sig_pre = [];  north_mat_sig_post = [];
west_mat_sig = [];  west_mat_sig_pre = [];  west_mat_sig_post = [];
south_mat_sig = [];  south_mat_sig_pre = [];  south_mat_sig_post = [];
east_mat_sig = [];  east_mat_sig_pre = [];  east_mat_sig_post = [];

% same for action
all_mat_sig_act = [];  all_mat_sig_act_pre = [];  all_mat_sig_act_post = [];
north_mat_sig_act = [];  north_mat_sig_act_pre = [];  north_mat_sig_act_post = [];
west_mat_sig_act = [];  west_mat_sig_act_pre = [];  west_mat_sig_act_post = [];
south_mat_sig_act = [];  south_mat_sig_act_pre = [];  south_mat_sig_act_post = [];
east_mat_sig_act = [];  east_mat_sig_act_pre = [];  east_mat_sig_act_post = [];

% same for reward
all_mat_sig_rew = [];  all_mat_sig_rew_pre = [];  all_mat_sig_rew_post = [];
north_mat_sig_rew = [];  north_mat_sig_rew_pre = [];  north_mat_sig_rew_post = [];
west_mat_sig_rew = [];  west_mat_sig_rew_pre = [];  west_mat_sig_rew_post = [];
south_mat_sig_rew = [];  south_mat_sig_rew_pre = [];  south_mat_sig_rew_post = [];
east_mat_sig_rew = [];  east_mat_sig_rew_pre = [];  east_mat_sig_rew_post = [];


% collect only significant cells
all_sig_mat = [];
north_sig_mat = [];
west_sig_mat = [];
south_sig_mat = [];
east_sig_mat = [];

% action
all_sig_act_mat = [];
north_sig_act_mat = [];
west_sig_act_mat = [];
south_sig_act_mat = [];
east_sig_act_mat = [];
% reward
all_sig_rew_mat = [];
north_sig_rew_mat = [];
west_sig_rew_mat = [];
south_sig_rew_mat = [];
east_sig_rew_mat = [];


all_subj_mat = [];
all_nCells_label = [];
all_sig_nCells_label = [];
all_sig_act_nCells_label = [];
all_sig_rew_nCells_label = [];

all_nCells_c_ord = []; % for session type colors.

all_FR_mat = [];
% % remove R1 due to oddness.
% R1_idx = find(contains(sessions,'R1'));
% sessions(R1_idx) = [];

for iS = 1:length(sessions)
    if isempty(all_out.(sessions{iS}))
        continue
    end
    
    
    
    if ~exist('all_sig_act', 'var')
        % keep only cells with any sig modulation
        sig_mean_out.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig.(sessions{iS}))),2);
        sig_mean_N_out.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig.(sessions{iS}))),2);
        sig_mean_W_out.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig.(sessions{iS}))),2);
        sig_mean_S_out.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig.(sessions{iS}))),2);
        sig_mean_E_out.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig.(sessions{iS}))),2);
        
        sig_mean_all_mat(iS,:) = sig_mean_out.(sessions{iS}); %#ok<*AGROW>
        sig_mean_N_mat(iS,:) = sig_mean_N_out.(sessions{iS});
        sig_mean_W_mat(iS,:) = sig_mean_W_out.(sessions{iS});
        sig_mean_S_mat(iS,:) = sig_mean_S_out.(sessions{iS});
        sig_mean_E_mat(iS,:) = sig_mean_E_out.(sessions{iS});
        
        % keep only cells with sig modulation after 'post' action.
        sig_mean_out_pre.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig_pre.(sessions{iS}))),2);
        sig_mean_N_out_pre.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig_pre.(sessions{iS}))),2);
        sig_mean_W_out_pre.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig_pre.(sessions{iS}))),2);
        sig_mean_S_out_pre.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig_pre.(sessions{iS}))),2);
        sig_mean_E_out_pre.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig_pre.(sessions{iS}))),2);
        
        sig_mean_all_mat_pre(iS,:) = sig_mean_out_pre.(sessions{iS}); %#ok<*AGROW>
        sig_mean_N_mat_pre(iS,:) = sig_mean_N_out_pre.(sessions{iS});
        sig_mean_W_mat_pre(iS,:) = sig_mean_W_out_pre.(sessions{iS});
        sig_mean_S_mat_pre(iS,:) = sig_mean_S_out_pre.(sessions{iS});
        sig_mean_E_mat_pre(iS,:) = sig_mean_E_out_pre.(sessions{iS});
        
        % keep only cells with sig modulation after 'post' action.
        sig_mean_out_post.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig_post.(sessions{iS}))),2);
        sig_mean_N_out_post.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig_post.(sessions{iS}))),2);
        sig_mean_W_out_post.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig_post.(sessions{iS}))),2);
        sig_mean_S_out_post.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig_post.(sessions{iS}))),2);
        sig_mean_E_out_post.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig_post.(sessions{iS}))),2);
        
        sig_mean_all_mat_post(iS,:) = sig_mean_out_post.(sessions{iS}); %#ok<*AGROW>
        sig_mean_N_mat_post(iS,:) = sig_mean_N_out_post.(sessions{iS});
        sig_mean_W_mat_post(iS,:) = sig_mean_W_out_post.(sessions{iS});
        sig_mean_S_mat_post(iS,:) = sig_mean_S_out_post.(sessions{iS});
        sig_mean_E_mat_post(iS,:) = sig_mean_E_out_post.(sessions{iS});
        
    elseif exist('all_mat_sig_act', 'var')
        sig_mean_out_act.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig_act.(sessions{iS}))),2);
        sig_mean_N_out_act.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig_act.(sessions{iS}))),2);
        sig_mean_W_out_act.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig_act.(sessions{iS}))),2);
        sig_mean_S_out_act.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig_act.(sessions{iS}))),2);
        sig_mean_E_out_act.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig_act.(sessions{iS}))),2);
        
        sig_mean_all_act_mat(iS,:) = sig_mean_out_act.(sessions{iS}); %#ok<*AGROW>
        sig_mean_N_act_mat(iS,:) = sig_mean_N_out_act.(sessions{iS});
        sig_mean_W_act_mat(iS,:) = sig_mean_W_out_act.(sessions{iS});
        sig_mean_S_act_mat(iS,:) = sig_mean_S_out_act.(sessions{iS});
        sig_mean_E_act_mat(iS,:) = sig_mean_E_out_act.(sessions{iS});
        
        % keep only cells with sig modulation after 'post' action.
        sig_mean_out_act_pre.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig_act_pre.(sessions{iS}))),2);
        sig_mean_N_out_act_pre.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig_act_pre.(sessions{iS}))),2);
        sig_mean_W_out_act_pre.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig_act_pre.(sessions{iS}))),2);
        sig_mean_S_out_act_pre.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig_act_pre.(sessions{iS}))),2);
        sig_mean_E_out_act_pre.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig_act_pre.(sessions{iS}))),2);
        
        sig_mean_all_mat_act_pre(iS,:) = sig_mean_out_act_pre.(sessions{iS}); %#ok<*AGROW>
        sig_mean_N_mat_act_pre(iS,:) = sig_mean_N_out_act_pre.(sessions{iS});
        sig_mean_W_mat_act_pre(iS,:) = sig_mean_W_out_act_pre.(sessions{iS});
        sig_mean_S_mat_act_pre(iS,:) = sig_mean_S_out_act_pre.(sessions{iS});
        sig_mean_E_mat_act_pre(iS,:) = sig_mean_E_out_act_pre.(sessions{iS});
        
        % keep only cells with sig modulation after 'post' action.
        sig_mean_out_act_post.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig_act_post.(sessions{iS}))),2);
        sig_mean_N_out_act_post.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig_act_post.(sessions{iS}))),2);
        sig_mean_W_out_act_post.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig_act_post.(sessions{iS}))),2);
        sig_mean_S_out_act_post.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig_act_post.(sessions{iS}))),2);
        sig_mean_E_out_act_post.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig_act_post.(sessions{iS}))),2);
        
        sig_mean_all_mat_act_post(iS,:) = sig_mean_out_act_post.(sessions{iS}); %#ok<*AGROW>
        sig_mean_N_mat_act_post(iS,:) = sig_mean_N_out_act_post.(sessions{iS});
        sig_mean_W_mat_act_post(iS,:) = sig_mean_W_out_act_post.(sessions{iS});
        sig_mean_S_mat_act_post(iS,:) = sig_mean_S_out_act_post.(sessions{iS});
        sig_mean_E_mat_act_post(iS,:) = sig_mean_E_out_act_post.(sessions{iS});
        
        %%%% same for reward
        sig_mean_out_rew.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig_rew.(sessions{iS}))),2);
        sig_mean_N_out_rew.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig_rew.(sessions{iS}))),2);
        sig_mean_W_out_rew.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig_rew.(sessions{iS}))),2);
        sig_mean_S_out_rew.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig_rew.(sessions{iS}))),2);
        sig_mean_E_out_rew.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig_rew.(sessions{iS}))),2);
        
        sig_mean_all_rew_mat(iS,:) = sig_mean_out_rew.(sessions{iS}); %#ok<*AGROW>
        sig_mean_N_rew_mat(iS,:) = sig_mean_N_out_rew.(sessions{iS});
        sig_mean_W_rew_mat(iS,:) = sig_mean_W_out_rew.(sessions{iS});
        sig_mean_S_rew_mat(iS,:) = sig_mean_S_out_rew.(sessions{iS});
        sig_mean_E_rew_mat(iS,:) = sig_mean_E_out_rew.(sessions{iS});
        
        % keep only cells with sig modulation after 'post' rewion.
        sig_mean_out_rew_pre.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig_rew_pre.(sessions{iS}))),2);
        sig_mean_N_out_rew_pre.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig_rew_pre.(sessions{iS}))),2);
        sig_mean_W_out_rew_pre.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig_rew_pre.(sessions{iS}))),2);
        sig_mean_S_out_rew_pre.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig_rew_pre.(sessions{iS}))),2);
        sig_mean_E_out_rew_pre.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig_rew_pre.(sessions{iS}))),2);
        
        sig_mean_all_mat_rew_pre(iS,:) = sig_mean_out_rew_pre.(sessions{iS}); %#ok<*AGROW>
        sig_mean_N_mat_rew_pre(iS,:) = sig_mean_N_out_rew_pre.(sessions{iS});
        sig_mean_W_mat_rew_pre(iS,:) = sig_mean_W_out_rew_pre.(sessions{iS});
        sig_mean_S_mat_rew_pre(iS,:) = sig_mean_S_out_rew_pre.(sessions{iS});
        sig_mean_E_mat_rew_pre(iS,:) = sig_mean_E_out_rew_pre.(sessions{iS});
        
        % keep only cells with sig modulation after 'post' rewion.
        sig_mean_out_rew_post.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig_rew_post.(sessions{iS}))),2);
        sig_mean_N_out_rew_post.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig_rew_post.(sessions{iS}))),2);
        sig_mean_W_out_rew_post.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig_rew_post.(sessions{iS}))),2);
        sig_mean_S_out_rew_post.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig_rew_post.(sessions{iS}))),2);
        sig_mean_E_out_rew_post.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig_rew_post.(sessions{iS}))),2);
        
        sig_mean_all_mat_rew_post(iS,:) = sig_mean_out_rew_post.(sessions{iS}); %#ok<*AGROW>
        sig_mean_N_mat_rew_post(iS,:) = sig_mean_N_out_rew_post.(sessions{iS});
        sig_mean_W_mat_rew_post(iS,:) = sig_mean_W_out_rew_post.(sessions{iS});
        sig_mean_S_mat_rew_post(iS,:) = sig_mean_S_out_rew_post.(sessions{iS});
        sig_mean_E_mat_rew_post(iS,:) = sig_mean_E_out_rew_post.(sessions{iS});
        
        
    end
    
    
    mean_out.(sessions{iS}) = nanmean(all_out.(sessions{iS}),2);
    mean_N_out.(sessions{iS}) = nanmean(north_out.(sessions{iS}),2);
    mean_W_out.(sessions{iS}) = nanmean(west_out.(sessions{iS}),2);
    mean_S_out.(sessions{iS}) = nanmean(south_out.(sessions{iS}),2);
    mean_E_out.(sessions{iS}) = nanmean(east_out.(sessions{iS}),2);
    
    
    mean_all_mat(iS,:) = mean_out.(sessions{iS});
    mean_N_mat(iS,:) = mean_N_out.(sessions{iS});
    mean_W_mat(iS,:) = mean_W_out.(sessions{iS});
    mean_S_mat(iS,:) = mean_S_out.(sessions{iS});
    mean_E_mat(iS,:) = mean_E_out.(sessions{iS});
    
    
    % get all the firing rates
    all_FR_mat = [all_FR_mat,  all_FR.(sessions{iS})];
    
    
    %collect all subject ids
    all_subj_mat = [all_subj_mat, all_subj.(sessions{iS})];
    
    % collect every cell.
    
    all_mat = [all_mat, all_out.(sessions{iS})];
    north_mat = [north_mat, north_out.(sessions{iS})];
    west_mat = [west_mat, west_out.(sessions{iS})];
    south_mat = [south_mat, south_out.(sessions{iS})];
    east_mat = [east_mat, east_out.(sessions{iS})];
    
            sess_label{iS} = sessions{iS};
        nCells(iS) = size(all_out.(sessions{iS}),2);
        
    if ~exist('all_mat_sig_act', 'var')
        % get the corresponding significane value (two sided)
        all_mat_sig = [all_mat_sig, all_sig.(sessions{iS})];
        north_mat_sig = [north_mat_sig, north_sig.(sessions{iS})];
        west_mat_sig = [west_mat_sig, west_sig.(sessions{iS})];
        south_mat_sig = [south_mat_sig, south_sig.(sessions{iS})];
        east_mat_sig = [east_mat_sig, east_sig.(sessions{iS})];
        
        % get the corresponding significane value (left side -> greater
        % before reward)
        all_mat_sig_pre = [all_mat_sig_pre, all_sig_pre.(sessions{iS})];
        north_mat_sig_pre = [north_mat_sig_pre, north_sig_pre.(sessions{iS})];
        west_mat_sig_pre = [west_mat_sig_pre, west_sig_pre.(sessions{iS})];
        south_mat_sig_pre = [south_mat_sig_pre, south_sig_pre.(sessions{iS})];
        east_mat_sig_pre = [east_mat_sig_pre, east_sig_pre.(sessions{iS})];
        
        
        % get the corresponding significane value (right side, greater post
        % reward)
        all_mat_sig_post = [all_mat_sig_post, all_sig_post.(sessions{iS})];
        north_mat_sig_post = [north_mat_sig_post, north_sig_post.(sessions{iS})];
        west_mat_sig_post = [west_mat_sig_post, west_sig_post.(sessions{iS})];
        south_mat_sig_post = [south_mat_sig_post, south_sig_post.(sessions{iS})];
        east_mat_sig_post = [east_mat_sig_post, east_sig_post.(sessions{iS})];
        
        
        % collect every significant cell.
        all_sig_mat = [all_sig_mat, all_out.(sessions{iS})(:,logical(all_sig.(sessions{iS})))];
        north_sig_mat = [north_sig_mat, north_out.(sessions{iS})(:,logical(north_sig.(sessions{iS})))];
        west_sig_mat = [west_sig_mat, west_out.(sessions{iS})(:,logical(west_sig.(sessions{iS})))];
        south_sig_mat = [south_sig_mat, south_out.(sessions{iS})(:,logical(south_sig.(sessions{iS})))];
        east_sig_mat = [east_sig_mat, east_out.(sessions{iS})(:,logical(east_sig.(sessions{iS})))];
        
        %     all_sig_nCells_label = [all_sig_nCells_label, repmat(iS,1, size(all_out.(sessions{iS})(:,logical(all_sig.(sessions{iS}))),2))];

        sig_all_nCells(iS) = size(all_out.(sessions{iS})(:,logical(all_sig.(sessions{iS}))),2);
        sig_N_nCells(iS) = size(north_out.(sessions{iS})(:,logical(north_sig.(sessions{iS}))),2);
        sig_W_nCells(iS) = size(west_out.(sessions{iS})(:,logical(west_sig.(sessions{iS}))),2);
        sig_S_nCells(iS) = size(south_out.(sessions{iS})(:,logical(south_sig.(sessions{iS}))),2);
        sig_E_nCells(iS) = size(east_out.(sessions{iS})(:,logical(east_sig.(sessions{iS}))),2);
        % pre only
        sig_all_nCells_pre(iS) = size(all_out.(sessions{iS})(:,logical(all_sig_pre.(sessions{iS}))),2);
        sig_N_nCells_pre(iS) = size(north_out.(sessions{iS})(:,logical(north_sig_pre.(sessions{iS}))),2);
        sig_W_nCells_pre(iS) = size(west_out.(sessions{iS})(:,logical(west_sig_pre.(sessions{iS}))),2);
        sig_S_nCells_pre(iS) = size(south_out.(sessions{iS})(:,logical(south_sig_pre.(sessions{iS}))),2);
        sig_E_nCells_pre(iS) = size(east_out.(sessions{iS})(:,logical(east_sig_pre.(sessions{iS}))),2);
        % post only
        sig_all_nCells_post(iS) = size(all_out.(sessions{iS})(:,logical(all_sig_post.(sessions{iS}))),2);
        sig_N_nCells_post(iS) = size(north_out.(sessions{iS})(:,logical(north_sig_post.(sessions{iS}))),2);
        sig_W_nCells_post(iS) = size(west_out.(sessions{iS})(:,logical(west_sig_post.(sessions{iS}))),2);
        sig_S_nCells_post(iS) = size(south_out.(sessions{iS})(:,logical(south_sig_post.(sessions{iS}))),2);
        sig_E_nCells_post(iS) = size(east_out.(sessions{iS})(:,logical(east_sig_post.(sessions{iS}))),2);
        
    elseif exist('all_mat_sig_act', 'var')
        % get the corresponding significane value (two sided)
        all_mat_sig_act = [all_mat_sig_act, all_sig_act.(sessions{iS})];
        north_mat_sig_act = [north_mat_sig_act, north_sig.(sessions{iS})];
        west_mat_sig_act = [west_mat_sig_act, west_sig_act.(sessions{iS})];
        south_mat_sig_act = [south_mat_sig_act, south_sig_act.(sessions{iS})];
        east_mat_sig_act = [east_mat_sig_act, east_sig_act.(sessions{iS})];
        
        % get the corresponding significane value (left side -> greater
        % before reward)
        all_mat_sig_act_pre = [all_mat_sig_act_pre, all_sig_act_pre.(sessions{iS})];
        north_mat_sig_act_pre = [north_mat_sig_act_pre, north_sig_act_pre.(sessions{iS})];
        west_mat_sig_act_pre = [west_mat_sig_act_pre, west_sig_act_pre.(sessions{iS})];
        south_mat_sig_act_pre = [south_mat_sig_act_pre, south_sig_act_pre.(sessions{iS})];
        east_mat_sig_act_pre = [east_mat_sig_act_pre, east_sig_act_pre.(sessions{iS})];
        
        
        % get the corresponding significane value (right side, greater post
        % reward)
        all_mat_sig_act_post = [all_mat_sig_act_post, all_sig_act_post.(sessions{iS})];
        north_mat_sig_act_post = [north_mat_sig_act_post, north_sig_act_post.(sessions{iS})];
        west_mat_sig_act_post = [west_mat_sig_act_post, west_sig_act_post.(sessions{iS})];
        south_mat_sig_act_post = [south_mat_sig_act_post, south_sig_act_post.(sessions{iS})];
        east_mat_sig_act_post = [east_mat_sig_act_post, east_sig_act_post.(sessions{iS})];
        
        % collect every significant cell.
        all_sig_act_mat = [all_sig_act_mat, all_out.(sessions{iS})(:,logical(all_sig_act.(sessions{iS})))];
        north_sig_act_mat = [north_sig_act_mat, north_out.(sessions{iS})(:,logical(north_sig_act.(sessions{iS})))];
        west_sig_act_mat = [west_sig_act_mat, west_out.(sessions{iS})(:,logical(west_sig_act.(sessions{iS})))];
        south_sig_act_mat = [south_sig_act_mat, south_out.(sessions{iS})(:,logical(south_sig_act.(sessions{iS})))];
        east_sig_act_mat = [east_sig_act_mat, east_out.(sessions{iS})(:,logical(east_sig_act.(sessions{iS})))];
        
        sig_all_nCells_act(iS) = size(all_out.(sessions{iS})(:,logical(all_sig_act.(sessions{iS}))),2);
        sig_N_nCells_act(iS) = size(north_out.(sessions{iS})(:,logical(north_sig_act.(sessions{iS}))),2);
        sig_W_nCells_act(iS) = size(west_out.(sessions{iS})(:,logical(west_sig_act.(sessions{iS}))),2);
        sig_S_nCells_act(iS) = size(south_out.(sessions{iS})(:,logical(south_sig_act.(sessions{iS}))),2);
        sig_E_nCells_act(iS) = size(east_out.(sessions{iS})(:,logical(east_sig_act.(sessions{iS}))),2);
        % pre only
        sig_all_nCells_act_pre(iS) = size(all_out.(sessions{iS})(:,logical(all_sig_act_pre.(sessions{iS}))),2);
        sig_N_nCells_act_pre(iS) = size(north_out.(sessions{iS})(:,logical(north_sig_act_pre.(sessions{iS}))),2);
        sig_W_nCells_act_pre(iS) = size(west_out.(sessions{iS})(:,logical(west_sig_act_pre.(sessions{iS}))),2);
        sig_S_nCells_act_pre(iS) = size(south_out.(sessions{iS})(:,logical(south_sig_act_pre.(sessions{iS}))),2);
        sig_E_nCells_act_pre(iS) = size(east_out.(sessions{iS})(:,logical(east_sig_act_pre.(sessions{iS}))),2);
        % post only
        sig_all_nCells_act_post(iS) = size(all_out.(sessions{iS})(:,logical(all_sig_act_post.(sessions{iS}))),2);
        sig_N_nCells_act_post(iS) = size(north_out.(sessions{iS})(:,logical(north_sig_act_post.(sessions{iS}))),2);
        sig_W_nCells_act_post(iS) = size(west_out.(sessions{iS})(:,logical(west_sig_act_post.(sessions{iS}))),2);
        sig_S_nCells_act_post(iS) = size(south_out.(sessions{iS})(:,logical(south_sig_act_post.(sessions{iS}))),2);
        sig_E_nCells_act_post(iS) = size(east_out.(sessions{iS})(:,logical(east_sig_act_post.(sessions{iS}))),2);
        
        
        %%% same for reward
         % get the corresponding significane value (two sided)
        all_mat_sig_rew = [all_mat_sig_rew, all_sig_rew.(sessions{iS})];
        north_mat_sig_rew = [north_mat_sig_rew, north_sig.(sessions{iS})];
        west_mat_sig_rew = [west_mat_sig_rew, west_sig_rew.(sessions{iS})];
        south_mat_sig_rew = [south_mat_sig_rew, south_sig_rew.(sessions{iS})];
        east_mat_sig_rew = [east_mat_sig_rew, east_sig_rew.(sessions{iS})];
        
        % get the corresponding significane value (left side -> greater
        % before reward)
        all_mat_sig_rew_pre = [all_mat_sig_rew_pre, all_sig_rew_pre.(sessions{iS})];
        north_mat_sig_rew_pre = [north_mat_sig_rew_pre, north_sig_rew_pre.(sessions{iS})];
        west_mat_sig_rew_pre = [west_mat_sig_rew_pre, west_sig_rew_pre.(sessions{iS})];
        south_mat_sig_rew_pre = [south_mat_sig_rew_pre, south_sig_rew_pre.(sessions{iS})];
        east_mat_sig_rew_pre = [east_mat_sig_rew_pre, east_sig_rew_pre.(sessions{iS})];
        
        
        % get the corresponding significane value (right side, greater post
        % reward)
        all_mat_sig_rew_post = [all_mat_sig_rew_post, all_sig_rew_post.(sessions{iS})];
        north_mat_sig_rew_post = [north_mat_sig_rew_post, north_sig_rew_post.(sessions{iS})];
        west_mat_sig_rew_post = [west_mat_sig_rew_post, west_sig_rew_post.(sessions{iS})];
        south_mat_sig_rew_post = [south_mat_sig_rew_post, south_sig_rew_post.(sessions{iS})];
        east_mat_sig_rew_post = [east_mat_sig_rew_post, east_sig_rew_post.(sessions{iS})];
        
        % collect every significant cell.
        all_sig_rew_mat = [all_sig_rew_mat, all_out.(sessions{iS})(:,logical(all_sig_rew.(sessions{iS})))];
        north_sig_rew_mat = [north_sig_rew_mat, north_out.(sessions{iS})(:,logical(north_sig_rew.(sessions{iS})))];
        west_sig_rew_mat = [west_sig_rew_mat, west_out.(sessions{iS})(:,logical(west_sig_rew.(sessions{iS})))];
        south_sig_rew_mat = [south_sig_rew_mat, south_out.(sessions{iS})(:,logical(south_sig_rew.(sessions{iS})))];
        east_sig_rew_mat = [east_sig_rew_mat, east_out.(sessions{iS})(:,logical(east_sig_rew.(sessions{iS})))];
        
        sig_all_nCells_rew(iS) = size(all_out.(sessions{iS})(:,logical(all_sig_rew.(sessions{iS}))),2);
        sig_N_nCells_rew(iS) = size(north_out.(sessions{iS})(:,logical(north_sig_rew.(sessions{iS}))),2);
        sig_W_nCells_rew(iS) = size(west_out.(sessions{iS})(:,logical(west_sig_rew.(sessions{iS}))),2);
        sig_S_nCells_rew(iS) = size(south_out.(sessions{iS})(:,logical(south_sig_rew.(sessions{iS}))),2);
        sig_E_nCells_rew(iS) = size(east_out.(sessions{iS})(:,logical(east_sig_rew.(sessions{iS}))),2);
        % pre only
        sig_all_nCells_rew_pre(iS) = size(all_out.(sessions{iS})(:,logical(all_sig_rew_pre.(sessions{iS}))),2);
        sig_N_nCells_rew_pre(iS) = size(north_out.(sessions{iS})(:,logical(north_sig_rew_pre.(sessions{iS}))),2);
        sig_W_nCells_rew_pre(iS) = size(west_out.(sessions{iS})(:,logical(west_sig_rew_pre.(sessions{iS}))),2);
        sig_S_nCells_rew_pre(iS) = size(south_out.(sessions{iS})(:,logical(south_sig_rew_pre.(sessions{iS}))),2);
        sig_E_nCells_rew_pre(iS) = size(east_out.(sessions{iS})(:,logical(east_sig_rew_pre.(sessions{iS}))),2);
        % post only
        sig_all_nCells_rew_post(iS) = size(all_out.(sessions{iS})(:,logical(all_sig_rew_post.(sessions{iS}))),2);
        sig_N_nCells_rew_post(iS) = size(north_out.(sessions{iS})(:,logical(north_sig_rew_post.(sessions{iS}))),2);
        sig_W_nCells_rew_post(iS) = size(west_out.(sessions{iS})(:,logical(west_sig_rew_post.(sessions{iS}))),2);
        sig_S_nCells_rew_post(iS) = size(south_out.(sessions{iS})(:,logical(south_sig_rew_post.(sessions{iS}))),2);
        sig_E_nCells_rew_post(iS) = size(east_out.(sessions{iS})(:,logical(east_sig_rew_post.(sessions{iS}))),2);
        
    end
    
    
    
    all_nCells_label = [all_nCells_label, repmat(iS,1, size(all_out.(sessions{iS}),2))];
    
    if strcmp(sessions{iS}(1), 'C')
        this_type = 1;
    elseif strcmp(sessions{iS}(1), 'O')
        this_type = 2;
    elseif strcmp(sessions{iS}(1), 'E')
        this_type = 3;
    elseif strcmp(sessions{iS}(1), 'R')
        this_type = 4;
    end
    
    all_nCells_c_ord = [all_nCells_c_ord, repmat(this_type,1, size(all_out.(sessions{iS}),2))];
    
    
end


%     %gaussian 1D smoothing
%     for jj = size(mean_all_mat, 1):-1:1
%         mean_all_mat(jj,:) = filter(gausswin(10),1,mean_all_mat(jj,:)); % conv(mean_all_mat(:,jj), ones(15,1)/15, 'same');
%     end


Z_min = min(min(mean_all_mat));% mean_N_mat; mean_W_mat; mean_S_mat; mean_W_mat]));
Z_max = max(max(mean_all_mat));% mean_N_mat; mean_W_mat; mean_S_mat; mean_W_mat]));

% if ~exist('all_mat_sig_act', 'var')
%     Z_sig_min = min(min([mean_all_mat; sig_mean_N_mat; sig_mean_W_mat; sig_mean_S_mat; sig_mean_W_mat]));
%     Z_sig_max = max(max([mean_all_mat; sig_mean_N_mat; sig_mean_W_mat; sig_mean_S_mat; sig_mean_W_mat]));
% elseif exist('all_mat_sig_act', 'var')
%     Z_sig_min_act = min(min([sig_mean_all_act_mat; sig_mean_N_mat; sig_mean_W_mat; sig_mean_S_mat; sig_mean_W_mat]));
%     Z_sig_max = max(max([sig_mean_all_act_mat; sig_mean_N_mat; sig_mean_W_mat; sig_mean_S_mat; sig_mean_W_mat]));
%     
%     Z_sig_min = min(min([mean_all_mat; sig_mean_N_mat; sig_mean_W_mat; sig_mean_S_mat; sig_mean_W_mat]));
%     Z_sig_max = max(max([mean_all_mat; sig_mean_N_mat; sig_mean_W_mat; sig_mean_S_mat; sig_mean_W_mat]));
% end


Z_cell_min = min(min(all_mat));
Z_cell_max = max(max(all_mat));

sess_idx = 1:length(sessions);
[all_mat_order, sort_idx] = sort(all_nCells_label, 'ascend');

c_ord = linspecer(4); % one for each session type C, O, E, R.
blues = parula(16); reds = jet(64); oranges = autumn(16);
sess_cord = [flipud(blues(3:2:7,:));(reds(end-7:end-1,:)); c_ord(3,:); flipud(oranges(end-9:end-6,:))];
% figure(1010)
% hold on
% for ii = 1:length(sess_cord)
%     plot(1:10, ones(1,10)*ii, 'color', sess_cord(ii,:), 'linewidth', 4)
% end
%% get an index for the deval type

for ii = length(all_subj_mat):-1:1
    if ismember(all_subj_mat(ii), [1 5 6 8]) && (all_nCells_label(ii) >= 12)
        N_deval_mat(ii) = 0;
        S_deval_mat(ii) = 0;
        E_deval_mat(ii) = 1;
        W_deval_mat(ii) = 1;
    elseif ismember(all_subj_mat(ii), [2, 3, 7, 9]) && (all_nCells_label(ii) >=12)
        N_deval_mat(ii) = 1;
        S_deval_mat(ii) = 1;
        E_deval_mat(ii) = 0;
        W_deval_mat(ii) = 0;
    else
        N_deval_mat(ii) = 0;
        S_deval_mat(ii) = 0;
        E_deval_mat(ii) = 0;
        W_deval_mat(ii) = 0;
    end
end

%% make some simple image imagesc using all cells

% first plot the mean Z score for each session (row)
figure(101)

subplot(1,5,1)
imagesc(tvec.C1, 1:size(mean_N_mat,1), mean_N_mat)
title({'North' ;'banana x 3'})
set(gca, 'ytick', 1:size(mean_N_mat,1), 'yTickLabel', sess_label)
set(gca, 'xtick', -5:2.5:5)
xlabel('feeder time (s)')
ylabel('session')
vline(0, 'k')
for ii = 1:size(mean_N_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_N_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
% colorbar
caxis([Z_min Z_max]);


subplot(1,5,2)
imagesc(tvec.C1, 1:size(mean_W_mat,1), mean_W_mat)
title({'West'; 'grain x 3'})
% set(gca, 'ytick', 1:size(mean_W_mat,1), 'yTickLabel', sess_label);
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
for ii = 1:size(mean_N_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_W_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
caxis([Z_min Z_max]);


subplot(1,5,3)
imagesc(tvec.C1, 1:size(mean_S_mat,1), mean_S_mat)
title({'South'; 'banana x 1'})
% set(gca, 'ytick', 1:size(mean_S_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
for ii = 1:size(mean_N_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_N_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
caxis([Z_min Z_max]);



subplot(1,5,4)
imagesc(tvec.C1, 1:size(mean_E_mat,1), mean_E_mat)
title({'East'; 'grain x 1'})
% set(gca, 'ytick', 1:size(mean_E_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
for ii = 1:size(mean_N_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_E_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
caxis([Z_min Z_max]);


subplot(1,5,5)
imagesc(tvec.C1, 1:size(mean_all_mat,1), mean_all_mat)
title('All feeders')
% set(gca, 'ytick', 1:size(mean_all_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
for ii = 1:size(mean_N_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_all_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
caxis([Z_min Z_max]);
cb=colorbar;
cb.Position(1) = cb.Position(1) + .075;
ylabel(cb, 'mean zscore','Rotation',90)


% add a main title
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 0.98,'Mean Zscore per session')
% title(currentFigure.Children(end), 'Mean Zscore per session');

SetFigure([], gcf)


%% same plot using only cells with sig modulation.
figure(102)

subplot(1,5,1)
imagesc(tvec.C1, 1:size(sig_mean_N_mat,1), sig_mean_N_mat)
title({'North' ;'banana x 3'})
set(gca, 'ytick', 1:size(sig_mean_N_mat,1), 'yTickLabel', sess_label)
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
xlabel('feeder time (s)')
ylabel('session')
vline(0, 'k')
for ii = 1:size(sig_mean_N_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_N_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
% colorbar
caxis([Z_sig_min Z_sig_max]);


subplot(1,5,2)
imagesc(tvec.C1, 1:size(sig_mean_W_mat,1), sig_mean_W_mat)
title({'West'; 'grain x 3'})
% set(gca, 'ytick', 1:size(mean_W_mat,1), 'yTickLabel', sess_label);
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_W_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_W_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);


subplot(1,5,3)
imagesc(tvec.C1, 1:size(sig_mean_S_mat,1), sig_mean_S_mat)
title({'South'; 'banana x 1'})
% set(gca, 'ytick', 1:size(mean_S_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
for ii = 1:size(sig_mean_S_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_S_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);



subplot(1,5,4)
imagesc(tvec.C1, 1:size(sig_mean_E_mat,1), sig_mean_E_mat)
title({'East'; 'grain x 1'})
% set(gca, 'ytick', 1:size(mean_E_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_E_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_E_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);


subplot(1,5,5)
imagesc(tvec.C1, 1:size(sig_mean_all_mat,1), sig_mean_all_mat)
title('All feeders')
% set(gca, 'ytick', 1:size(mean_all_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_all_mat,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_all_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);
cb=colorbar;
cb.Position(1) = cb.Position(1) + .075;
ylabel(cb, 'mean zscore','Rotation',90)


% add a main title
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 0.98,'Mean Zscore per session')
% title(currentFigure.Children(end), 'Mean Zscore per session');

SetFigure([], gcf)


%% %% same plot using only cells with sig modulation. split by pos or negative.
figure(1025)
UP  = char(11016);
DOWN = char(11018);

sub_map = reshape(1:80, 10, 8)';
% north post
subplot(8,10,sub_map(:,1))
imagesc(tvec.C1, 1:size(sig_mean_N_mat_post,1), sig_mean_N_mat_post)
title({'North' ;'banana x 3'; UP})
set(gca, 'ytick', 1:size(sig_mean_N_mat_post,1), 'yTickLabel', sess_label)
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
xlabel('feeder time (s)')
ylabel('session')
vline(0, 'k')
for ii = 1:size(sig_mean_N_mat_post,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_N_nCells_post(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);

% north Pre
subplot(8,10,sub_map(:,2))
imagesc(tvec.C1, 1:size(sig_mean_N_mat_pre,1), sig_mean_N_mat_pre)
title({'North' ;'banana x 3'; DOWN})
set(gca, 'ytick', [], 'yTickLabel', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_N_mat_pre,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_N_nCells_pre(ii)) ], 'fontsize', 12);
end
% colorbar
caxis([Z_sig_min Z_sig_max]);



%%%% WEST post
subplot(8,10,sub_map(:,3))
imagesc(tvec.C1, 1:size(sig_mean_W_mat_post,1), sig_mean_W_mat_post)
title({'West'; 'grain x 3'; UP})
% set(gca, 'ytick', 1:size(mean_W_mat,1), 'yTickLabel', sess_label);
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_W_mat_post,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_W_nCells_post(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);


subplot(8,10,sub_map(:,4))
imagesc(tvec.C1, 1:size(sig_mean_W_mat_pre,1), sig_mean_W_mat_pre)
title({'West'; 'grain x 3'; DOWN})
% set(gca, 'ytick', 1:size(mean_W_mat,1), 'yTickLabel', sess_label);
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_W_mat_pre,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_W_nCells_pre(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);



% SOUTH post
subplot(8,10,sub_map(:,5))
imagesc(tvec.C1, 1:size(sig_mean_S_mat_post,1), sig_mean_S_mat_post)
title({'South'; 'banana x 1'; UP})
% set(gca, 'ytick', 1:size(mean_W_mat,1), 'yTickLabel', sess_label);
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_S_mat_post,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_S_nCells_post(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);

subplot(8,10,sub_map(:,6))
imagesc(tvec.C1, 1:size(sig_mean_S_mat_pre,1), sig_mean_S_mat_pre)
title({'South'; 'banana x 1'; DOWN})
% set(gca, 'ytick', 1:size(mean_S_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_S_mat_pre,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_S_nCells_pre(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);



% EAST post
subplot(8,10,sub_map(:,7))
imagesc(tvec.C1, 1:size(sig_mean_E_mat_post,1), sig_mean_E_mat_post)
title({'East'; 'grain x 1'; UP})
% set(gca, 'ytick', 1:size(mean_W_mat,1), 'yTickLabel', sess_label);
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_E_mat_post,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_E_nCells_post(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);

subplot(8,10,sub_map(:,8))
imagesc(tvec.C1, 1:size(sig_mean_E_mat_pre,1), sig_mean_E_mat_pre)
title({'East'; 'grain x 1'; DOWN})
% set(gca, 'ytick', 1:size(mean_S_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_E_mat_pre,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_E_nCells_pre(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);


% ALL post
subplot(8,10,sub_map(:,9))
imagesc(tvec.C1, 1:size(sig_mean_all_mat_post,1), sig_mean_all_mat_post)
title({'All feeders'; UP})
% set(gca, 'ytick', 1:size(mean_all_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_all_mat_post,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_all_nCells_post(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);


subplot(8,10,sub_map(:,10))
imagesc(tvec.C1, 1:size(sig_mean_all_mat_pre,1), sig_mean_all_mat_pre)
title({'All feeders'; DOWN})
% set(gca, 'ytick', 1:size(mean_all_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:max(tvec.C1)+mode(diff(tvec.C1)))
vline(0, 'k')
for ii = 1:size(sig_mean_all_mat_pre,1)
    text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, [num2str(sig_all_nCells_pre(ii))], 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);

cb=colorbar;
cb.Position(1) = cb.Position(1) + .075;
ylabel(cb, 'mean zscore','Rotation',90)


% add a main title
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 0.98,'Mean Zscore per session')
% title(currentFigure.Children(end), 'Mean Zscore per session');
cfg_fig.ft_size = 14;
SetFigure(cfg_fig, gcf)
%% Again with all of the cells stacked.

figure(103)

subplot(1,5,1)
imagesc(tvec.C1, 1:length(all_mat_order), north_mat(:,sort_idx)')
title({'North' ;'banana x 3'})
% set(gca, 'ytick', 1:length(all_mat_order), 'yTickLabel', sessions(all_mat_order))
set(gca, 'ytick', [])

set(gca, 'xtick', -5:2.5:5)
xlabel('feeder time (s)')
ylabel('session')
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_cell_min Z_cell_max]);
for ii = 1:length(all_mat_order)
    text(min(tvec.C1)-mode(diff(tvec.C1))*5, ii, sessions(all_mat_order(ii)), 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left');
    if north_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif north_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
    if N_deval_mat(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*6, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end


hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end



%
subplot(1,5,2)
imagesc(tvec.C1, 1:length(all_mat_order), west_mat(:,sort_idx)')
% set(gca, 'ytick', 1:length(all_mat_order), 'yTickLabel', sessions(all_mat_order))
set(gca, 'ytick', [])
title({'West'; 'grain x 3'})
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:2.5)
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_cell_min Z_cell_max]);
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
for ii = 1:length(all_mat_order)
    if west_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif west_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
    if W_deval_mat(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*6, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end


subplot(1,5,3)
imagesc(tvec.C1, 1:length(all_mat_order), south_mat(:,sort_idx)')
% set(gca, 'ytick', 1:length(all_mat_order), 'yTickLabel', sessions(all_mat_order))
set(gca, 'ytick', [])
title({'South'; 'banana x 1'})
% set(gca, 'ytick', 1:size(mean_S_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_cell_min Z_cell_max]);
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end

for ii = 1:length(all_mat_order)
    if south_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif south_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
    if S_deval_mat(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*6, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end


subplot(1,5,4)
imagesc(tvec.C1, 1:length(all_mat_order), east_mat(:,sort_idx)')
% set(gca, 'ytick', 1:length(all_mat_order), 'yTickLabel', sessions(all_mat_order))
set(gca, 'ytick', [])
title({'East'; 'grain x 1'})
% set(gca, 'ytick', 1:size(mean_E_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_cell_min Z_cell_max]);
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
for ii = 1:length(all_mat_order)
    if east_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif east_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
    if E_deval_mat(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*6, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color',sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end

subplot(1,5,5)
imagesc(tvec.C1, 1:length(all_mat_order), all_mat(:,sort_idx)')
% set(gca, 'ytick', 1:length(all_mat_order), 'yTickLabel', sessions(all_mat_order))
set(gca, 'ytick', [])
title('All feeders')
% set(gca, 'ytick', 1:size(mean_all_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_cell_min Z_cell_max]);
cb=colorbar;
cb.Position(1) = cb.Position(1) + .075;
ylabel(cb, 'mean zscore','Rotation',90)
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
for ii = 1:length(all_mat_order)
    if all_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif all_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', sess_cord(all_mat_order(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end

% add a main title
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 0.98,'Mean Zscore per session')
% title(currentFigure.Children(end), 'Mean Zscore per session');

cfg_fig.resize = 2;
SetFigure(cfg_fig, gcf)



%% get the feeder activity vs end of matrix
% turn_idx(1) = nearest_idx(-1, tvec.C1);
% turn_idx(2) = nearest_idx(1, tvec.C1);
%
% feeder_idx(1) = nearest_idx(1, tvec.C1);
% feeder_idx(2) = length(tvec.C1);
base_idx(1) = 1;
base_idx(2) = nearest_idx(-1, tvec.C1) -1;

turn_idx(1) = nearest_idx(-1, tvec.C1);
turn_idx(2) = nearest_idx(0.8, tvec.C1);

feeder_idx(1) = nearest_idx(0.8, tvec.C1);
feeder_idx(2) = nearest_idx(1.8, tvec.C1);


% get the mean zscore at the feeder point (3 : 5 s) and subtract the mean
% zscore at the decision point (-1 : 1s)
all_feed_turn_val = (mean(all_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(all_mat(turn_idx(1):turn_idx(2),:),1))/2;
N_feed_turn_val = (mean(north_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(north_mat(turn_idx(1):turn_idx(2),:),1))/2;
E_feed_turn_val = (mean(east_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(east_mat(turn_idx(1):turn_idx(2),:),1))/2;
S_feed_turn_val = (mean(south_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(south_mat(turn_idx(1):turn_idx(2),:),1))/2;
W_feed_turn_val = (mean(west_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(west_mat(turn_idx(1):turn_idx(2),:),1))/2;

all_feed_turn_val  =  - abs(mean(all_mat(turn_idx(1):turn_idx(2),:),1)) + abs(mean(all_mat(feeder_idx(1):feeder_idx(2),:),1)); %abs(mean(all_mat(base_idx(1):base_idx(2),:),1))


figure(1010)
% NORTH
subplot(1,30,1:4)
imagesc(tvec.C1, 1:length(all_mat_order), north_mat(:,sort_idx)')
title({'North' ;'banana x 3'})
% set(gca, 'ytick', 1:length(all_mat_order), 'yTickLabel', sessions(all_mat_order))
set(gca, 'ytick', [])

set(gca, 'xtick', min(tvec.C1):2.5:2.5)
xlabel('feeder time (s)')
ylabel('session')
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_sig_min Z_sig_max]);
for ii = 1:length(all_mat_order)
    text(min(tvec.C1)-(mode(diff(tvec.C1))*5), ii, sessions(all_mat_order(ii)), 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left');
    if north_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif north_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
    if N_deval_mat(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*6, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end

subplot(1,30, 6)
% imagesc(1, 1:length(all_mat_order), N_feed_turn_val(:,sort_idx)')
scatter(N_feed_turn_val(:,sort_idx), (1:length(all_mat_order))-.5,125, N_feed_turn_val(:,sort_idx), 'filled');
% colormap(gca,sunsmap)
set(gca, 'YDir', 'reverse')
title({'Rew' ;  'idx'})
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []);
ylim([.5 length(all_mat_order)-.5]);


% WEST
subplot(1,30, 7:10)
imagesc(tvec.C1, 1:length(all_mat_order), west_mat(:,sort_idx)')
set(gca, 'ytick', [])
title({'West'; 'grain x 3'})
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:2.5)
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_sig_min Z_sig_max]);
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
for ii = 1:length(all_mat_order)
    if west_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif west_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
    
    if W_deval_mat(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*6, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end

subplot(1,30, 12)
% imagesc(1, 1:length(all_mat_order), W_feed_turn_val(:,sort_idx)');
scatter(W_feed_turn_val(:,sort_idx), (1:length(all_mat_order))-.5,125, W_feed_turn_val(:,sort_idx), 'filled')
set(gca, 'YDir', 'reverse')
title({'Rew' ;  'idx'})
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []);
ylim([.5 length(all_mat_order)-.5]);


% SOUTH
subplot(1,30, 13:16)
imagesc(tvec.C1, 1:length(all_mat_order), south_mat(:,sort_idx)')
set(gca, 'ytick', [])
title({'South'; 'banana x 1'})
% set(gca, 'ytick', 1:size(mean_S_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:2.5)
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_sig_min Z_sig_max]);
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end

for ii = 1:length(all_mat_order)
    if south_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif south_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
    if S_deval_mat(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*6, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end

subplot(1,30, 18)
% imagesc(1, 1:length(all_mat_order), S_feed_turn_val(:,sort_idx)')
scatter(S_feed_turn_val(:,sort_idx), (1:length(all_mat_order))-.5,125, S_feed_turn_val(:,sort_idx), 'filled')
set(gca, 'YDir', 'reverse')
title({'Rew' ;  'idx'})
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []);
ylim([.5 length(all_mat_order)-.5]);




% EAST
subplot(1,30, 19:22)
imagesc(tvec.C1, 1:length(all_mat_order), east_mat(:,sort_idx)')
set(gca, 'ytick', [])
title({'East'; 'grain x 1'})
% set(gca, 'ytick', 1:size(mean_E_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_sig_min Z_sig_max]);
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
for ii = 1:length(all_mat_order)
    if east_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif east_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
    if E_deval_mat(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*6, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end

subplot(1,30,24)
% imagesc(tvec.C1, 1:length(all_mat_order), E_feed_turn_val(:,sort_idx)')
scatter(E_feed_turn_val(:,sort_idx),(1:length(all_mat_order))-.5,125, E_feed_turn_val(:,sort_idx), 'filled')
set(gca, 'YDir', 'reverse')
title({'Rew' ;  'idx'})
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []);
ylim([.5 length(all_mat_order)-.5]);


% ALL
subplot(1,30,25:28)
imagesc(tvec.C1, 1:length(all_mat_order), all_mat(:,sort_idx)')
title('All')
set(gca, 'ytick', [])
title('All feeders')
% set(gca, 'ytick', 1:size(mean_all_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', min(tvec.C1):2.5:2.5)
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_sig_min Z_sig_max]);
cb=colorbar;
cb.Position(1) = cb.Position(1) + .105;
ylabel(cb, 'mean zscore','Rotation',90)
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
for ii = 1:length(all_mat_order)
    if all_mat_sig_pre(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*2, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    elseif all_mat_sig_post(ii)
        text(max(tvec.C1)+mode(diff(tvec.C1))*4, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
    end
end

subplot(1,30,30)
% imagesc(1, 1:length(all_mat_order), all_feed_turn_val(:,sort_idx)')
scatter(all_feed_turn_val(:,sort_idx), (1:length(all_mat_order))-.5,125, all_feed_turn_val(:,sort_idx), 'filled')
set(gca, 'YDir', 'reverse')
title({'Rew' ;  'idx'})
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []);
ylim([.5 length(all_mat_order)-.5]);


SetFigure([], gcf)


% save the output data as a table/csv
for ii = 1:length(all_mat_order)
    sess_id(ii) = sessions(all_mat_order(ii));
end

index_out = array2table([N_feed_turn_val(:,sort_idx)', W_feed_turn_val(:,sort_idx)', S_feed_turn_val(:,sort_idx)', E_feed_turn_val(:,sort_idx)', all_feed_turn_val(:,sort_idx)'], 'VariableNames',{'North_index',  'West_index', 'South_index', 'East_index', 'All_index'}) ;

index_out = [cell2table(sess_id'), index_out];

writetable( index_out,'Rew_turn_idx.csv')

%% save all the figures
if exist('smooth_val', 'var')
    figure(101)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_mean_smooth.png'])
    saveas(gcf, [inter_dir filesep 'All_mean_smooth.fig'])
    print(gcf,[inter_dir filesep  'All_mean_smooth'],'-depsc')
    
    figure(102)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_mean_sig_only_smooth.png'])
    saveas(gcf, [inter_dir filesep 'All_mean_sig_only_smooth.fig'])
    print(gcf,[inter_dir filesep  'All_mean_sig_only_smooth'],'-depsc')
    
    figure(1025)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_mean_sig_only_split_smooth.png'])
    saveas(gcf, [inter_dir filesep 'All_mean_sig_only_split_smooth.fig'])
    print(gcf,[inter_dir filesep  'All_mean_sig_only_split_smooth'],'-depsc')
    
    
    figure(103)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_cells_smooth.png'])
    saveas(gcf, [inter_dir filesep 'All_cells_smooth.fig'])
    print(gcf,[inter_dir filesep  'All_cells_smooth'],'-depsc')
    
    figure(1010)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_cells_mod_idx_smooth.png'])
    saveas(gcf, [inter_dir filesep 'All_cells_mod_idx_smooth.fig'])
    print(gcf,[inter_dir filesep  'All_cells_mod_idx_smooth'],'-depsc')
    
else
    figure(101)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_mean.png'])
    saveas(gcf, [inter_dir filesep 'All_mean.fig'])
    print(gcf,[inter_dir filesep  'All_mean'],'-depsc')
    
    figure(102)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_mean_sig_only.png'])
    saveas(gcf, [inter_dir filesep 'All_mean_sig_only.fig'])
    print(gcf,[inter_dir filesep  'All_mean_sig_only'],'-depsc')
    
    figure(1025)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_mean_sig_only_split.png'])
    saveas(gcf, [inter_dir filesep 'All_mean_sig_only_split.fig'])
    print(gcf,[inter_dir filesep  'All_mean_sig_only_split'],'-depsc')
    
    
    figure(103)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_cells.png'])
    saveas(gcf, [inter_dir filesep 'All_cells.fig'])
    print(gcf,[inter_dir filesep  'All_cells'],'-depsc')
    
    figure(1010)
    maximize
    
    saveas(gcf, [inter_dir filesep 'All_cells_mod_idx.png'])
    saveas(gcf, [inter_dir filesep 'All_cells_mod_idx.fig'])
    print(gcf,[inter_dir filesep  'All_cells_mod_idx'],'-depsc')
    
end



%% bar plots for % change

Cs = [all_out_perc.C1(all_sig_pre.C1 | all_sig_post.C1),...
    all_out_perc.C2(all_sig_pre.C2 | all_sig_post.C2),...
    all_out_perc.C3(all_sig_pre.C3 | all_sig_post.C3)];
O_e = [all_out_perc.O1(all_sig_pre.O1 | all_sig_post.O1),...
    all_out_perc.O2(all_sig_pre.O2 | all_sig_post.O2),...
    all_out_perc.O3(all_sig_pre.O3 | all_sig_post.O3),...
    all_out_perc.O4(all_sig_pre.O4 | all_sig_post.O4)];
O_l = [all_out_perc.O5(all_sig_pre.O5 | all_sig_post.O5),...
    all_out_perc.O6(all_sig_pre.O6 | all_sig_post.O6),...
    all_out_perc.O7(all_sig_pre.O7 | all_sig_post.O7)];
Rs = [all_out_perc.R1(all_sig_pre.R1 | all_sig_post.R1),...
    all_out_perc.R2(all_sig_pre.R2 | all_sig_post.R2),...
    all_out_perc.R3(all_sig_pre.R3 | all_sig_post.R3)];


% figure(302)
% % bar([Cs; O_e; O_l; Rs]);
%
% set(gca, 'xticklabel', {'Crit', 'Over Early', 'Over Late', 'Rein'});
% hline(1)

% export as csv
mat_out = NaN(4, max([length(Cs),length(O_e),length(O_l),length(Rs)]));


mat_out(1,1:length(Cs)) = Cs;
mat_out(2,1:length(O_e)) = O_e;
mat_out(3,1:length(O_l)) = O_l;
mat_out(4,1:length(Rs)) = Rs;

mat_out= mat_out';

csvwrite('block_perc.csv', mat_out)


%% same thing for rew-turn index
% fnames = fieldnames(all_out);
% Cs_idx = fnames

Cs = [N_feed_turn_val, W_feed_turn_val, S_feed_turn_val, E_feed_turn_val all_feed_turn_val] ;
O_e = [all_out_perc.O1(all_sig_pre.O1 | all_sig_post.O1),...
    all_out_perc.O2(all_sig_pre.O2 | all_sig_post.O2),...
    all_out_perc.O3(all_sig_pre.O3 | all_sig_post.O3),...
    all_out_perc.O4(all_sig_pre.O4 | all_sig_post.O4)];
O_l = [all_out_perc.O5(all_sig_pre.O5 | all_sig_post.O5),...
    all_out_perc.O6(all_sig_pre.O6 | all_sig_post.O6),...
    all_out_perc.O7(all_sig_pre.O7 | all_sig_post.O7)];
Rs = [all_out_perc.R1(all_sig_pre.R1 | all_sig_post.R1),...
    all_out_perc.R2(all_sig_pre.R2 | all_sig_post.R2),...
    all_out_perc.R3(all_sig_pre.R3 | all_sig_post.R3)];


% figure(302)
% % bar([Cs; O_e; O_l; Rs]);
%
% set(gca, 'xticklabel', {'Crit', 'Over Early', 'Over Late', 'Rein'});
% hline(1)

% export as csv
mat_out = NaN(4, max([length(Cs),length(O_e),length(O_l),length(Rs)]));


mat_out(1,1:length(Cs)) = Cs;
mat_out(2,1:length(O_e)) = O_e;
mat_out(3,1:length(O_l)) = O_l;
mat_out(4,1:length(Rs)) = Rs;

mat_out= mat_out';

csvwrite('block_perc.csv', mat_out)

%%  unrelated Kenny shadded error bar thing

% Initialize variables.
filename = '/Users/jericcarmichael/Dropbox (Personal)/KA_Data/CeA_record_train.csv';
delimiter = ',';
startRow = 2;

% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r','n','UTF-8');
% Skip the BOM (Byte Order Mark).
fseek(fileID, 3, 'bof');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

% Close the text file.
fclose(fileID);

% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,4,5,6,7,8,9,10,11]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


% Split data into numeric and string columns.
rawNumericColumns = raw(:, [2,4,5,6,7,8,9,10,11]);
rawStringColumns = string(raw(:, [1,3]));


% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
for catIdx = [1,2]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end

% Create output variable
CeArecordtrain1 = table;
CeArecordtrain1.rat = categorical(rawStringColumns(:, 1));
CeArecordtrain1.day = cell2mat(rawNumericColumns(:, 1));
CeArecordtrain1.sex = categorical(rawStringColumns(:, 2));
CeArecordtrain1.time = cell2mat(rawNumericColumns(:, 2));
CeArecordtrain1.total = cell2mat(rawNumericColumns(:, 3));
CeArecordtrain1.correct = cell2mat(rawNumericColumns(:, 4));
CeArecordtrain1.incorrect = cell2mat(rawNumericColumns(:, 5));
CeArecordtrain1.outcome = cell2mat(rawNumericColumns(:, 6));
CeArecordtrain1.perc = cell2mat(rawNumericColumns(:, 7));
CeArecordtrain1.min = cell2mat(rawNumericColumns(:, 8));
CeArecordtrain1.speed = cell2mat(rawNumericColumns(:, 9));

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns catIdx idx;


close all

%% plot stuff avec error bars
days = unique(CeArecordtrain1.day);
for ii = 1:length(days)
    
    perc_mean(ii) = mean(CeArecordtrain1.perc(CeArecordtrain1.day == days(ii)));
    perc_SEM(ii) = std(CeArecordtrain1.perc(CeArecordtrain1.day == days(ii))) / sqrt(length(CeArecordtrain1.perc(CeArecordtrain1.day == days(ii))));
    
    
end
figure;
h = shadedErrorBar(days, perc_mean*100, perc_SEM*100);
h.mainLine.Color = c_ord(1,:);
h.mainLine.LineWidth = 3;
h.edge(1).Color = h.patch.FaceColor;
h.edge(2).Color = h.patch.FaceColor;

hline(60)

% bar(CeArecordtrain1.day, CeArecordtrain1.perc')%, {@mean,@std})


%% Clustering.
sc_cord =[c_ord(2,:); [139, 114, 142]/255; c_ord(1,:); c_ord(1,:)];
% generate a mean FR value for each bin.
bins = -2.5:(5/3):2.5;
t_idx = nearest_idx3(bins, tvec.C1);
tvec.C1(t_idx);
t_idx(end) = length(tvec.C1);

binned_vec = [];
for ii = size(all_mat, 2):-1:1
    
    for jj = length(t_idx):-1:2
        binned_vec(ii,jj-1) = mean(all_mat(t_idx(jj-1):t_idx(jj),ii));
    end
    
    scat_cord(ii,:)  = sess_cord(all_mat_order(ii),:);
    %     if (all_mat_order(ii) > 4) && (all_mat_order(ii) < 7)
    %                 scat_g_cord(ii,:) = sc_cord(3, :);
    %     else
    scat_g_cord(ii,:) = sc_cord(all_nCells_c_ord(ii), :);
    %     end
    ratio(ii) = binned_vec(ii, 3) - binned_vec(ii,2);
    act_ratio(ii) = binned_vec(ii,2) - binned_vec(ii,1);
end

% k-means

% k-means
% data_in = [firing_rate', burst_idx', rise_fall_inter'];

% [g_idx, n_idx, frq_idx]= MS_kmean_scatter([binned_vec(:, 1),binned_vec(:, 2),binned_vec(:, 3)], 3,[3,2,1], 50);
%
%
% xref = linspace(floor(min(all_FR_mat)),ceil(max(all_FR_mat)),length(all_FR_mat));
% indexed_x = interp1(all_FR_mat,0:128,'nearest');
% cmap = parula(length(0:128));
% xcolors = cmap(indexed_x,:);
figure(909)
hold on
c_group = [scat_cord(1,:); scat_cord(80,:); scat_cord(109,:); scat_cord(136,:)];
for kk = 1:4
    hL = plot(min(binned_vec(:, 1)):max(binned_vec(:, 1)), nan(2, 1),'o', 'color', sc_cord(kk,:),'MarkerFaceColor', sc_cord(kk,:));
end
% hs = scatter3(ratio(all_nCells_c_ord ~= 4),act_ratio(all_nCells_c_ord ~= 4),all_FR_mat(all_nCells_c_ord ~= 4), 120,scat_g_cord(all_nCells_c_ord ~= 4,:), 'filled');
% hs = scatter(act_ratio(all_nCells_c_ord ~= 4), ratio(all_nCells_c_ord ~= 4),120,scat_g_cord(all_nCells_c_ord ~= 4,:), 'filled');

hs = scatter(all_feed_turn_val(:,sort_idx), all_FR_mat, 120, scat_g_cord, 'filled')
legend({'C1-1'; 'O1-3'; 'O4-7'})
% view([-145, 25])
grid on
% [g_idx, n_idx, frq_idx]= MS_kmean_scatter([all_FR_mat',ratio',binned_vec(:, 3)], 3,[3,2,1], 50);

% zlabel('overall FR (Hz)');
ylabel('rew - act ');
xlabel('act - baseline');

% xlabel('previous (-2.5 : -0.8s) Zscore');
% ylabel('action (-0.8 : 0.8s) Zscore');
% zlabel('reward (0.8 : 2.5s) Zscore');

SetFigure([], gcf)

%
% figure(606)
% subplot(2,2,1)
% scatter3(ratio(all_,act_ratio,all_FR_mat, 120,scat_g_cord, 'filled')

%% percentage of cells with reward mod.
% overall
fprintf('<strong>Overall Reward mod: %2.2f%%</strong>\n', (sum(all_mat_sig_post)/length(all_mat_sig_post))*100);
fprintf('<strong>C1-3 Reward mod: %2.2f%%</strong>\n', (sum(all_mat_sig_post(all_nCells_c_ord ==1))/length(all_mat_sig_post((all_nCells_c_ord ==1))))*100);
O13 = ismember(all_mat_order, [4:6]);
fprintf('<strong>O1-3 Reward mod: %2.2f%%</strong>\n', (sum(all_mat_sig_post(O13))/length(all_mat_sig_post(O13)))*100);
O47 = ismember(all_mat_order, [7 8 9 10]);
fprintf('<strong>O4-7 Reward mod: %2.2f%%</strong>\n',  (sum(all_mat_sig_post(O47))/length(all_mat_sig_post(O47)))*100);
fprintf('<strong>R1-3 Reward mod: %2.2f%%</strong>\n', (sum(all_mat_sig_post(all_nCells_c_ord ==4))/length(all_mat_sig_post((all_nCells_c_ord ==4))))*100);



%%
fprintf('\n')
% overall
any_sig_mod = (all_mat_sig_act | all_mat_sig_rew); 
fprintf('<strong>Overall Reward or Action mod: %2.2f%%</strong>\n', (sum(any_sig_mod)/length(any_sig_mod))*100);
fprintf('<strong>C1-3 Reward or Action mod: %2.2f%%</strong>\n', (sum(any_sig_mod(all_nCells_c_ord ==1))/length(any_sig_mod((all_nCells_c_ord ==1))))*100);
O13 = ismember(all_mat_order, [4:6]);
fprintf('<strong>O1-3 Reward or Action mod: %2.2f%%</strong>\n', (sum(any_sig_mod(O13))/length(any_sig_mod(O13)))*100);
O47 = ismember(all_mat_order, [7 8 9 10]);
fprintf('<strong>O4-7 Reward or Action mod: %2.2f%%</strong>\n',  (sum(any_sig_mod(O47))/length(any_sig_mod(O47)))*100);
fprintf('<strong>R1-3 Reward or Action mod: %2.2f%%</strong>\n', (sum(any_sig_mod(all_nCells_c_ord ==4))/length(any_sig_mod((all_nCells_c_ord ==4))))*100);
fprintf('_____________________\n')


any_sig_mod = (all_mat_sig_act & all_mat_sig_rew); 
fprintf('<strong>Overall Reward & Action mod: %2.2f%%</strong>\n', (sum(any_sig_mod)/length(any_sig_mod))*100);
fprintf('<strong>C1-3 Reward & Action mod: %2.2f%%</strong>\n', (sum(any_sig_mod(all_nCells_c_ord ==1))/length(any_sig_mod((all_nCells_c_ord ==1))))*100);
O13 = ismember(all_mat_order, [4:6]);
fprintf('<strong>O1-3 Reward & Action mod: %2.2f%%</strong>\n', (sum(any_sig_mod(O13))/length(any_sig_mod(O13)))*100);
O47 = ismember(all_mat_order, [7 8 9 10]);
fprintf('<strong>O4-7 Reward & Action mod: %2.2f%%</strong>\n',  (sum(any_sig_mod(O47))/length(any_sig_mod(O47)))*100);
fprintf('<strong>R1-3 Reward & Action mod: %2.2f%%</strong>\n', (sum(any_sig_mod(all_nCells_c_ord ==4))/length(any_sig_mod((all_nCells_c_ord ==4))))*100);

fprintf('_____________________\n')

% overall
fprintf('<strong>Overall Action mod: %2.2f%%</strong>\n', (sum(all_mat_sig_rew)/length(all_mat_sig_rew))*100);
fprintf('<strong>C1-3 Action mod: %2.2f%%</strong>\n', (sum(all_mat_sig_rew(all_nCells_c_ord ==1))/length(all_mat_sig_rew((all_nCells_c_ord ==1))))*100);
O13 = ismember(all_mat_order, [4:6]);
fprintf('<strong>O1-3 Action mod: %2.2f%%</strong>\n', (sum(all_mat_sig_rew(O13))/length(all_mat_sig_rew(O13)))*100);
O47 = ismember(all_mat_order, [7 8 9 10]);
fprintf('<strong>O4-7 Action mod: %2.2f%%</strong>\n',  (sum(all_mat_sig_rew(O47))/length(all_mat_sig_rew(O47)))*100);
fprintf('<strong>R1-3 Action mod: %2.2f%%</strong>\n', (sum(all_mat_sig_rew(all_nCells_c_ord ==4))/length(all_mat_sig_rew((all_nCells_c_ord ==4))))*100);

fprintf('_____________________\n')
% overall
fprintf('<strong>Overall Reward mod: %2.2f%%</strong>\n', (sum(all_mat_sig_act)/length(all_mat_sig_act))*100);
fprintf('<strong>C1-3 Reward mod: %2.2f%%</strong>\n', (sum(all_mat_sig_act(all_nCells_c_ord ==1))/length(all_mat_sig_act((all_nCells_c_ord ==1))))*100);
O13 = ismember(all_mat_order, [4:6]);
fprintf('<strong>O1-3 Reward mod: %2.2f%%</strong>\n', (sum(all_mat_sig_act(O13))/length(all_mat_sig_act(O13)))*100);
O47 = ismember(all_mat_order, [7 8 9 10]);
fprintf('<strong>O4-7 Reward mod: %2.2f%%</strong>\n',  (sum(all_mat_sig_act(O47))/length(all_mat_sig_act(O47)))*100);
fprintf('<strong>R1-3 Reward mod: %2.2f%%</strong>\n', (sum(all_mat_sig_act(all_nCells_c_ord ==4))/length(all_mat_sig_act((all_nCells_c_ord ==4))))*100);


%% counts

fprintf('\n')
% overall
any_sig_mod = (all_mat_sig_act | all_mat_sig_rew); 
fprintf('<strong>Overall Reward or Action mod: %2.0f</strong>\n', (sum(any_sig_mod)));
fprintf('<strong>C1-3 Reward or Action mod: %2.0f</strong>\n', (sum(any_sig_mod(all_nCells_c_ord ==1))));
O13 = ismember(all_mat_order, [4:6]);
fprintf('<strong>O1-3 Reward or Action mod: %2.0f</strong>\n', (sum(any_sig_mod(O13))));
O47 = ismember(all_mat_order, [7 8 9 10]);
fprintf('<strong>O4-7 Reward or Action mod: %2.0f</strong>\n',  (sum(any_sig_mod(O47))));
fprintf('<strong>R1-3 Reward or Action mod: %2.0f</strong>\n', (sum(any_sig_mod(all_nCells_c_ord ==4))));
fprintf('_____________________\n')


any_sig_mod = (all_mat_sig_act & all_mat_sig_rew); 
fprintf('<strong>Overall Reward & Action mod: %2.0f</strong>\n', (sum(any_sig_mod)));
fprintf('<strong>C1-3 Reward & Action mod: %2.0f</strong>\n', (sum(any_sig_mod(all_nCells_c_ord ==1))));
O13 = ismember(all_mat_order, [4:6]);
fprintf('<strong>O1-3 Reward & Action mod: %2.0f</strong>\n', (sum(any_sig_mod(O13))));
O47 = ismember(all_mat_order, [7 8 9 10]);
fprintf('<strong>O4-7 Reward & Action mod: %2.0f</strong>\n',  (sum(any_sig_mod(O47))));
fprintf('<strong>R1-3 Reward & Action mod: %2.0f</strong>\n', (sum(any_sig_mod(all_nCells_c_ord ==4))));

fprintf('_____________________\n')

% overall
fprintf('<strong>Overall Action mod: %2.0f</strong>\n', (sum(all_mat_sig_rew)));
fprintf('<strong>C1-3 Action mod: %2.0f</strong>\n', (sum(all_mat_sig_rew(all_nCells_c_ord ==1))));
O13 = ismember(all_mat_order, [4:6]);
fprintf('<strong>O1-3 Action mod: %2.0f</strong>\n', (sum(all_mat_sig_rew(O13))));
O47 = ismember(all_mat_order, [7 8 9 10]);
fprintf('<strong>O4-7 Action mod: %2.0f</strong>\n',  (sum(all_mat_sig_rew(O47))));
fprintf('<strong>R1-3 Action mod: %2.0f</strong>\n', (sum(all_mat_sig_rew(all_nCells_c_ord ==4))));

fprintf('_____________________\n')
% overall
fprintf('<strong>Overall Reward mod: %2.0f</strong>\n', sum(all_mat_sig_act));
fprintf('<strong>C1-3 Reward mod: %2.0f</strong>\n', (sum(all_mat_sig_act(all_nCells_c_ord ==1))));
O13 = ismember(all_mat_order, [4:6]);
fprintf('<strong>O1-3 Reward mod: %2.0f</strong>\n', (sum(all_mat_sig_act(O13))));
O47 = ismember(all_mat_order, [7 8 9 10]);
fprintf('<strong>O4-7 Reward mod: %2.0f</strong>\n',  (sum(all_mat_sig_act(O47))));
fprintf('<strong>R1-3 Reward mod: %2.0f</strong>\n', (sum(all_mat_sig_act(all_nCells_c_ord ==4))));



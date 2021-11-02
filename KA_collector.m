function [all_out] = KA_collector(inter_dir)

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

%% allocate some empty arrays to fill in. 
% for ii = 1:15 % Acqusition stage could talke many days. Will trim empty ones later. 
%     all_out.(['A' num2str(ii)]) = [];
%     all_sig.(['A' num2str(ii)]) = [];
% 
% end
subj_list = {'C1_1', 'C3_2','C3_3', 'C3_4', 'C4_3', 'C5_2', 'C5_3'}; 
deval_type = {'W', 'N', 'N', ' ' , 'W', 'W', 'N'}; 
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

all_out_perc  = all_out; 

all_subj = all_out; 

% cycle and collect

cd(inter_dir); 
file_names = FindFiles('*TT*');  % get all the files containing 'TT' which should be all the outputs from KA_master. 


for iF = 1:length(file_names)
    load(file_names{iF}); % load the data as 'This_cell'
    
    % put additional criteria for inclusion here. 
    
%     if str2double(file_names{iF}(strfind(file_names{iF}, '.mat')-1)) > 3
%         fprintf('Quality is %d, skipping %s\n', str2double(file_names{iF}(strfind(file_names{iF}, '.mat')-1)) , file_names{iF})
%         continue
%     end
%     
    
    this_sess = This_cell.session; % get the session type. 
    all_subj.(this_sess) =[all_subj.(this_sess), find(ismember(subj_list, This_cell.subject))]; 

    all_out.(this_sess) = [all_out.(this_sess) This_cell.mean_S_gau{5}]; 
    north_out.(this_sess) = [north_out.(this_sess) This_cell.mean_S_gau{1}];
    west_out.(this_sess) = [west_out.(this_sess) This_cell.mean_S_gau{2}];
    south_out.(this_sess) = [south_out.(this_sess) This_cell.mean_S_gau{3}]; 
    east_out.(this_sess) = [east_out.(this_sess) This_cell.mean_S_gau{4}]; 

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

% collect only significant cells
all_sig_mat = [];
north_sig_mat = [];
west_sig_mat = [];
south_sig_mat = [];
east_sig_mat = []; 

all_subj_mat = [];
all_nCells_label = []; 
all_sig_nCells_label = []; 
all_nCells_c_ord = []; % for session type colors. 


% % remove R1 due to oddness. 
% R1_idx = find(contains(sessions,'R1'));
% sessions(R1_idx) = []; 

for iS = 1:length(sessions)
    if isempty(all_out.(sessions{iS}))
        continue
    end
    % keep only cells with sig modulation
    sig_mean_out.(sessions{iS}) = nanmean(all_out.(sessions{iS})(:,logical(all_sig.(sessions{iS}))),2);
    sig_mean_N_out.(sessions{iS}) = nanmean(north_out.(sessions{iS})(:,logical(north_sig.(sessions{iS}))),2);
    sig_mean_W_out.(sessions{iS}) = nanmean(west_out.(sessions{iS})(:,logical(west_sig.(sessions{iS}))),2);
    sig_mean_S_out.(sessions{iS}) = nanmean(south_out.(sessions{iS})(:,logical(south_sig.(sessions{iS}))),2);
    sig_mean_E_out.(sessions{iS}) = nanmean(east_out.(sessions{iS})(:,logical(east_sig.(sessions{iS}))),2);
    
    sig_mean_all_mat(iS,:) = sig_mean_out.(sessions{iS});
    sig_mean_N_mat(iS,:) = sig_mean_N_out.(sessions{iS});
    sig_mean_W_mat(iS,:) = sig_mean_W_out.(sessions{iS});
    sig_mean_S_mat(iS,:) = sig_mean_S_out.(sessions{iS});
    sig_mean_E_mat(iS,:) = sig_mean_E_out.(sessions{iS});


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
    
    %collect all subject ids
    all_subj_mat = [all_subj_mat, all_subj.(sessions{iS})]; 
    
    % collect every cell. 
    all_mat = [all_mat, all_out.(sessions{iS})];
    north_mat = [north_mat, north_out.(sessions{iS})];
    west_mat = [west_mat, west_out.(sessions{iS})];
    south_mat = [south_mat, south_out.(sessions{iS})];
    east_mat = [east_mat, east_out.(sessions{iS})];
    
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

    % collect every significant cell.
    all_sig_mat = [all_sig_mat, all_out.(sessions{iS})(:,logical(all_sig.(sessions{iS})))];
    north_sig_mat = [north_sig_mat, north_out.(sessions{iS})(:,logical(north_sig.(sessions{iS})))];
    west_sig_mat = [west_sig_mat, west_out.(sessions{iS})(:,logical(west_sig.(sessions{iS})))];
    south_sig_mat = [south_sig_mat, south_out.(sessions{iS})(:,logical(south_sig.(sessions{iS})))];
    east_sig_mat = [east_sig_mat, east_out.(sessions{iS})(:,logical(east_sig.(sessions{iS})))];

%     all_sig_nCells_label = [all_sig_nCells_label, repmat(iS,1, size(all_out.(sessions{iS})(:,logical(all_sig.(sessions{iS}))),2))]; 

    
sess_label{iS} = sessions{iS}; 
nCells(iS) = size(all_out.(sessions{iS}),2);
sig_all_nCells(iS) = size(all_out.(sessions{iS})(:,logical(all_sig.(sessions{iS}))),2);
sig_N_nCells(iS) = size(north_out.(sessions{iS})(:,logical(north_sig.(sessions{iS}))),2);
sig_W_nCells(iS) = size(west_out.(sessions{iS})(:,logical(west_sig.(sessions{iS}))),2);
sig_S_nCells(iS) = size(south_out.(sessions{iS})(:,logical(south_sig.(sessions{iS}))),2);
sig_E_nCells(iS) = size(east_out.(sessions{iS})(:,logical(east_sig.(sessions{iS}))),2);

end

Z_min = min(min([mean_all_mat]))%; mean_N_mat; mean_W_mat; mean_S_mat; mean_W_mat])); 
Z_max = max(max([mean_all_mat]))%; mean_N_mat; mean_W_mat; mean_S_mat; mean_W_mat])); 

Z_sig_min = min(min([mean_all_mat; sig_mean_N_mat; sig_mean_W_mat; sig_mean_S_mat; sig_mean_W_mat])); 
Z_sig_max = max(max([mean_all_mat; sig_mean_N_mat; sig_mean_W_mat; sig_mean_S_mat; sig_mean_W_mat])); 

[all_mat_order, sort_idx] = sort(all_nCells_label, 'ascend');

c_ord = linspecer(length(unique(all_nCells_c_ord))); 

%% get an index for the deval type

for ii = 1:length(all_subj_mat)
    if ismember(all_subj_mat(ii), [1 5 6]) && (all_nCells_label(ii) >= 12)
        N_deval_mat(ii) = 0;
        S_deval_mat(ii) = 0;
        E_deval_mat(ii) = 1;
        W_deval_mat(ii) = 1;
    elseif ismember(all_subj_mat(ii), [2, 3, 6]) && (all_nCells_label(ii) >=12)
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
   text(5.5, ii, [num2str(sig_N_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
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
   text(5.5, ii, [num2str(sig_W_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
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
   text(5.5, ii, [num2str(sig_N_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
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
   text(5.5, ii, [num2str(sig_E_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
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
   text(5.5, ii, [num2str(sig_all_nCells(ii)) '/' num2str(nCells(ii))], 'fontsize', 12);
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
set(gca, 'xtick', -5:2.5:5)
xlabel('feeder time (s)')
ylabel('session')
vline(0, 'k')
for ii = 1:size(sig_mean_N_mat,1)
   text(5.5, ii, num2str(sig_N_nCells(ii)), 'fontsize', 12);
end
% colorbar
caxis([Z_min Z_max]);


subplot(1,5,2)
imagesc(tvec.C1, 1:size(sig_mean_W_mat,1), sig_mean_W_mat)
title({'West'; 'grain x 3'})
% set(gca, 'ytick', 1:size(mean_W_mat,1), 'yTickLabel', sess_label);
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
for ii = 1:size(sig_mean_W_mat,1)
   text(5.5, ii, num2str(sig_W_nCells(ii)), 'fontsize', 12);
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
   text(5.5, ii, num2str(sig_S_nCells(ii)), 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);



subplot(1,5,4)
imagesc(tvec.C1, 1:size(sig_mean_E_mat,1), sig_mean_E_mat)
title({'East'; 'grain x 1'})
% set(gca, 'ytick', 1:size(mean_E_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
for ii = 1:size(sig_mean_E_mat,1)
   text(5.5, ii, num2str(sig_E_nCells(ii)), 'fontsize', 12);
end
caxis([Z_sig_min Z_sig_max]);


subplot(1,5,5)
imagesc(tvec.C1, 1:size(sig_mean_all_mat,1), sig_mean_all_mat)
title('All feeders')
% set(gca, 'ytick', 1:size(mean_all_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
for ii = 1:size(sig_mean_all_mat,1)
   text(5.5, ii, num2str(sig_all_nCells(ii)), 'fontsize', 12);
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
caxis([Z_min Z_max]);
for ii = 1:length(all_mat_order)
   text(-5.8, ii, sessions(all_mat_order(ii)), 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left');
   if north_mat_sig_pre(ii)
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif north_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
   if N_deval_mat(ii)
       text(6.25, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
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
   if west_mat_sig_pre(ii)
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif west_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
   
   if W_deval_mat(ii)
       text(6.25, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
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
caxis([Z_sig_min Z_sig_max]);
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end

for ii = 1:length(all_mat_order)
   if south_mat_sig_pre(ii)
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif south_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
   if S_deval_mat(ii)
       text(6.25, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
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
caxis([Z_sig_min Z_sig_max]);
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
for ii = 1:length(all_mat_order)
   if east_mat_sig_pre(ii)
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif east_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
   if E_deval_mat(ii)
       text(6.25, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
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
caxis([Z_sig_min Z_sig_max]);
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
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif all_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
end

% add a main title
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 0.98,'Mean Zscore per session')
% title(currentFigure.Children(end), 'Mean Zscore per session');

cfg_fig.resize = 2;
SetFigure(cfg_fig, gcf)



%% get the feeder activity vs end of matrix 
turn_idx(1) = nearest_idx(-1, tvec.C1);
turn_idx(2) = nearest_idx(1, tvec.C1);

feeder_idx(1) = nearest_idx(3, tvec.C1);
feeder_idx(2) = length(tvec.C1);

% get the mean zscore at the feeder point (3 : 5 s) and subtract the mean
% zscore at the decision point (-1 : 1s)
all_feed_turn_val = (mean(all_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(all_mat(turn_idx(1):turn_idx(2),:),1))/2;
N_feed_turn_val = (mean(north_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(north_mat(turn_idx(1):turn_idx(2),:),1))/2;
E_feed_turn_val = (mean(east_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(east_mat(turn_idx(1):turn_idx(2),:),1))/2;
S_feed_turn_val = (mean(south_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(south_mat(turn_idx(1):turn_idx(2),:),1))/2;
W_feed_turn_val = (mean(west_mat(feeder_idx(1):feeder_idx(2),:),1) - mean(west_mat(turn_idx(1):turn_idx(2),:),1))/2;


figure(1010)
subplot(1,30,1:4)
imagesc(tvec.C1, 1:length(all_mat_order), north_mat(:,sort_idx)')
title({'North' ;'banana x 3'})
% set(gca, 'ytick', 1:length(all_mat_order), 'yTickLabel', sessions(all_mat_order))
set(gca, 'ytick', [])

set(gca, 'xtick', -5:2.5:5)
xlabel('feeder time (s)')
ylabel('session')
vl = vline(0, '--k');
vl.LineWidth = 2;
caxis([Z_min Z_max]);
for ii = 1:length(all_mat_order)
   text(-5.8, ii, sessions(all_mat_order(ii)), 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left');
   if north_mat_sig_pre(ii)
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif north_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
   if N_deval_mat(ii)
       text(6.25, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
end
hl = hline(find(diff(all_nCells_c_ord))+.5, {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end


subplot(1,30, 6)
% imagesc(1, 1:length(all_mat_order), N_feed_turn_val(:,sort_idx)')
scatter(N_feed_turn_val(:,sort_idx), (1:length(all_mat_order))-.5,125, N_feed_turn_val(:,sort_idx), 'filled')
set(gca, 'YDir', 'reverse')
title('North Rew - turn')
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []); 


subplot(1,30, 7:10)
imagesc(tvec.C1, 1:length(all_mat_order), west_mat(:,sort_idx)')
set(gca, 'ytick', [])
title({'West'; 'grain x 3'})
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
   if west_mat_sig_pre(ii)
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif west_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
   
   if W_deval_mat(ii)
       text(6.25, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
end




subplot(1,30, 12)
% imagesc(1, 1:length(all_mat_order), W_feed_turn_val(:,sort_idx)');
scatter(W_feed_turn_val(:,sort_idx), (1:length(all_mat_order))-.5,125, W_feed_turn_val(:,sort_idx), 'filled')
set(gca, 'YDir', 'reverse')
title('West Rew - turn')
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []); 


subplot(1,30, 13:16)
imagesc(tvec.C1, 1:length(all_mat_order), south_mat(:,sort_idx)')
set(gca, 'ytick', [])
title({'South'; 'banana x 1'})
% set(gca, 'ytick', 1:size(mean_S_mat,1), 'yTickLabel', sess_label)
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
   if south_mat_sig_pre(ii)
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif south_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
   if S_deval_mat(ii)
       text(6.25, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
end




subplot(1,30, 18)
% imagesc(1, 1:length(all_mat_order), S_feed_turn_val(:,sort_idx)')
scatter(S_feed_turn_val(:,sort_idx), (1:length(all_mat_order))-.5,125, S_feed_turn_val(:,sort_idx), 'filled')
set(gca, 'YDir', 'reverse')
title('South Rew - turn')
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []); 


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
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif east_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
   if E_deval_mat(ii)
       text(6.25, ii, '\diamondsuit', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
end

subplot(1,30,24)
% imagesc(tvec.C1, 1:length(all_mat_order), E_feed_turn_val(:,sort_idx)')
scatter(E_feed_turn_val(:,sort_idx),(1:length(all_mat_order))-.5,125, E_feed_turn_val(:,sort_idx), 'filled')
set(gca, 'YDir', 'reverse')
title('East Rew - turn')
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []); 


subplot(1,30,25:28)
imagesc(tvec.C1, 1:length(all_mat_order), all_mat(:,sort_idx)')
title('All')
set(gca, 'ytick', [])
title('All feeders')
% set(gca, 'ytick', 1:size(mean_all_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
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
      text(5.25, ii, '\downarrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   elseif all_mat_sig_post(ii)
      text(5.75, ii, '\uparrow', 'fontweight', 'bold','fontsize', 10, 'color', c_ord(all_nCells_c_ord(ii),:), 'HorizontalAlignment', 'left','Interpreter','tex');
   end
end


subplot(1,30,30)
% imagesc(1, 1:length(all_mat_order), all_feed_turn_val(:,sort_idx)')
scatter(all_feed_turn_val(:,sort_idx), (1:length(all_mat_order))-.5,125, all_feed_turn_val(:,sort_idx), 'filled')
set(gca, 'YDir', 'reverse')
title('All Rew - turn')
hl = hline(find(diff(all_nCells_c_ord)), {'k', 'k', 'k'});
for ii = 1:length(hl)
    hl(ii).LineWidth = 3;
    hl(ii).Color = c_ord(ii+1,:);
end
set(gca, 'xtick', [- 2 0 2], 'xticklabel', {'T','0', 'R'}, 'ytick', []); 


SetFigure([], gcf)


% save the output data as a table/csv
for ii = 1:length(all_mat_order)
sess_id(ii) = sessions(all_mat_order(ii)); 
end

index_out = array2table([N_feed_turn_val(:,sort_idx)', W_feed_turn_val(:,sort_idx)', S_feed_turn_val(:,sort_idx)', E_feed_turn_val(:,sort_idx)'], 'VariableNames',{'North_index',  'West_index', 'South_index', 'East_index'}) ;

index_out = [cell2table(sess_id'), index_out]; 

writetable( index_out,'Rew_turn_idx.csv')

%% save all the figures
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


figure(103)
maximize

saveas(gcf, [inter_dir filesep 'All_cells.png'])
saveas(gcf, [inter_dir filesep 'All_cells.fig'])
print(gcf,[inter_dir filesep  'All_cells'],'-depsc')


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



%% plot stuff avec error bars
days = unique(CeArecordtrain1.day);
for ii = 1:length(days)
    
    perc_mean(ii) = mean(CeArecordtrain1.perc(CeArecordtrain1.day == days(ii)));
    perc_SEM(ii) = std(CeArecordtrain1.perc(CeArecordtrain1.day == days(ii))) / sqrt(length(CeArecordtrain1.perc(CeArecordtrain1.day == days(ii)))); 
    
    
end

h = shadedErrorBar(days, perc_mean*100, perc_SEM*100)
h.mainLine.Color = c_ord(1,:);
h.mainLine.LineWidth = 3;
h.edge(1).Color = h.patch.FaceColor; 
h.edge(2).Color = h.patch.FaceColor; 

hline(60)

% bar(CeArecordtrain1.day, CeArecordtrain1.perc')%, {@mean,@std})


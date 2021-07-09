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


%% allocate some empty arrays to fill in. 
% for ii = 1:15 % Acqusition stage could talke many days. Will trim empty ones later. 
%     all_out.(['A' num2str(ii)]) = [];
%     all_sig.(['A' num2str(ii)]) = [];
% 
% end

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
%% cycle and collect

cd(inter_dir); 
file_names = FindFiles('*TT*');  % get all the files containing 'TT' which should be all the outputs from KA_master. 


for iF = 1:length(file_names)
    load(file_names{iF}); % load the data as 'This_cell'
    
    % put additional criteria for inclusion here. 
    
    if str2double(file_names{iF}(strfind(file_names{iF}, '.mat')-1)) > 3
        fprintf('Quality is %d, skipping %s\n', str2double(file_names{iF}(strfind(file_names{iF}, '.mat')-1)) , file_names{iF})
        continue
    end
    
    this_sess = This_cell.session; % get the session type. 
    
    all_out.(this_sess) = [all_out.(this_sess) This_cell.Z{5}]; 
    north_out.(this_sess) = [north_out.(this_sess) This_cell.Z{1}];
    west_out.(this_sess) = [west_out.(this_sess) This_cell.Z{2}];
    south_out.(this_sess) = [south_out.(this_sess) This_cell.Z{3}]; 
    east_out.(this_sess) = [east_out.(this_sess) This_cell.Z{4}]; 

    all_sig.(this_sess) =  [all_sig.(this_sess) This_cell.H{5}]; 
    
    north_sig.(this_sess) =  [north_sig.(this_sess) This_cell.H{1}]; 
    west_sig.(this_sess) =  [west_sig.(this_sess) This_cell.H{2}]; 
    south_sig.(this_sess) =  [south_sig.(this_sess) This_cell.H{3}]; 
    east_sig.(this_sess) =  [east_sig.(this_sess) This_cell.H{4}]; 

    tvec.(this_sess) = This_cell.outputIT{1}; 
    clear This_cell
end % iF files

%% compile into sessions and get the mean Z score. 
sessions = fieldnames(all_out); 

% remove R1 due to oddness. 
R1_idx = find(contains(sessions,'R1'));
sessions(R1_idx) = []; 

for iS = length(sessions):-1:1
    if isempty(all_out.(sessions{iS}))
        continue
    end
    mean_out.(sessions{iS}) = nanmean(all_out.(sessions{iS}),2);
    mean_N_out.(sessions{iS}) = nanmean(north_out.(sessions{iS}),2);
    mean_W_out.(sessions{iS}) = nanmean(west_out.(sessions{iS}),2);
    mean_S_out.(sessions{iS}) = nanmean(south_out.(sessions{iS}),2);
    mean_E_out.(sessions{iS}) = nanmean(east_out.(sessions{iS}),2);
    
mean_all_mat(iS,:) = mean_out.(sessions{iS}) ;
mean_N_mat(iS,:) = mean_N_out.(sessions{iS}) ;
mean_W_mat(iS,:) = mean_W_out.(sessions{iS}) ;
mean_S_mat(iS,:) = mean_S_out.(sessions{iS}) ;
mean_E_mat(iS,:) = mean_E_out.(sessions{iS}) ;

sess_label{iS} = sessions{iS}; 
nCells(iS) = size(all_out.(sessions{iS}),2);
end

Z_min = min(min([mean_all_mat; mean_N_mat; mean_W_mat; mean_S_mat; mean_W_mat])); 
Z_max = max(max([mean_all_mat; mean_N_mat; mean_W_mat; mean_S_mat; mean_W_mat])); 


%% make some simple image imagesc

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
   text(5.5, ii, num2str(nCells(ii)), 'fontsize', 12);
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
caxis([Z_min Z_max]);


subplot(1,5,3)
imagesc(tvec.C1, 1:size(mean_S_mat,1), mean_S_mat)
title({'South'; 'banana x 1'})
% set(gca, 'ytick', 1:size(mean_S_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
caxis([Z_min Z_max]);



subplot(1,5,4)
imagesc(tvec.C1, 1:size(mean_E_mat,1), mean_E_mat)
title({'East'; 'grain x 1'})
% set(gca, 'ytick', 1:size(mean_E_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
caxis([Z_min Z_max]);


subplot(1,5,5)
imagesc(tvec.C1, 1:size(mean_all_mat,1), mean_all_mat)
title('All feeders')
% set(gca, 'ytick', 1:size(mean_all_mat,1), 'yTickLabel', sess_label)
set(gca, 'ytick', [])
set(gca, 'xtick', -5:2.5:5)
vline(0, 'k')
caxis([Z_min Z_max]);
cb=colorbar;
cb.Position(1) = cb.Position(1) + .075;
ylabel(cb, 'mean zscore','Rotation',90)


% add a main title
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,'Mean Zscore per session')
% title(currentFigure.Children(end), 'Mean Zscore per session');

SetFigure([], gcf)

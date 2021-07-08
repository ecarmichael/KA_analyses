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


% allocate some empty arrays to fill in. 
for ii = 1:15 % Acqusition stage could talke many days. Will trim empty ones later. 
    all_out.(['A' num2str(ii)]) = [];
    all_sig.(['A' num2str(ii)]) = [];

end

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

    
    clear This_cell
end % iF files



%% KA screener master script. 

% load data
addpath(genpath('/Users/jericcarmichael/Documents/Github/EC_State'));
addpath(genpath('/Users/jericcarmichael/Documents/Github/KA_analyses'));
addpath(genpath('/Users/jericcarmichael/Documents/Github/vandermeerlab/code-matlab/shared'))
data_dir = '/Users/jericcarmichael/Dropbox/KA_Data/Raw_data'; % where all the NLX data is. 
inter_dir = '/Users/jericcarmichael/Dropbox/KA_Data/inter';  % where to save the outputs. 

cd(data_dir); % move to the data dir. 

% make an intermediate directory if it doesn't exist. 
if ~exist(inter_dir,'dir')
    mkdir(inter_dir)
end
    
%% loop over sessions / cells

% get all the sessions
this_dir = dir(data_dir);
for ii = 1:length(this_dir)
    if strcmp(this_dir(ii).name(1), '.') % check for hidden dirs 
        continue
    else
        sess_list{ii} = this_dir(ii).name;
    end
end
sess_list =   sess_list(~cellfun('isempty',sess_list));


% loop over sessions in the data dir.
for iS = 1:length(sess_list)
    
    cd([data_dir filesep sess_list{iS}])
    
    % example
    
    cells_to_process = FindFiles('*.t');
    
    % check if there are any .t files.  if not continue. 
    if isempty(cells_to_process)
        continue
    end
    
    for iT = 1:length(cells_to_process)
        
        parts = strsplit(cells_to_process{iT}, filesep);
        this_file = parts{end};
        This_cell = KA_screener(this_file);
        
        % if there were too few spikes in the .t then skip this file. 
        if isempty(This_cell)
            fprintf('<strong>%s</strong>:  Cell %s had too few spikes to be included (<1200)', mfilename, this_file)
            continue
        end
        
        parts = strsplit(sess_list{iS}, '_');
        This_cell.subject = [parts{1} '_' parts{2}]; % get the subject ID
        This_cell.session = parts{3}; % get the session type and number Acquisition, Criteria, Overtrain, Extinction, Reacquisition. nSessions: A ? C 3, O 7, E 1, R 3-5
        This_cell.date = parts{end}; 
        
        save([inter_dir filesep sess_list{iS} '_' this_file(1:end-2) '.mat'], 'This_cell')
        
        close all
    end
    
end
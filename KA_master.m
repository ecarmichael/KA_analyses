%% KA screener master script. 

% load data
addpath(genpath('/Users/jericcarmichael/Documents/Github/EC_State'));
addpath(genpath('/Users/jericcarmichael/Documents/Github/KA_analyses'));
addpath(genpath('/Users/jericcarmichael/Documents/Github/vandermeerlab/code-matlab/shared'))
data_dir = '/Users/jericcarmichael/Dropbox/KA_Data/Raw_data';
inter_dir = '/Users/jericcarmichael/Dropbox/KA_Data/inter'; 

cd(data_dir)

% make an intermediate directory
if ~exist(inter_dir,'dir')
    mkdir(inter_dir)
end
    
%% loop over sessions / cells

% get all the sessions
this_dir = dir(data_dir);
for ii = 3:length(this_dir)
    sess_list{ii} = this_dir(ii).name;
end
sess_list =   sess_list(~cellfun('isempty',sess_list));


% loop over sessions in the data dir.

for iS = 1:length(sess_list)
    
    cd([data_dir filesep sess_list{iS}])
    
    % example
    
    cells_to_process = FindFiles('*.t');
    
    for iT = 1:length(cells_to_process)
        
        parts = strsplit(cells_to_process{iT}, filesep);
        this_file = parts{end};
        This_cell = KA_screener(this_file);
        
        parts = strsplit(sess_list{iS}, '_');
        This_cell.subject = [parts{1} '_' parts{2}]; % get the subject ID
        This_cell.session = parts{3}; % get the session type and number Acquisition, Criteria, Overtrain, Extinction, Reacquisition
        This_cell.date = parts{end}; 
        
        save([inter_dir filesep sess_list{iS} '_' this_file], 
        
        close all
    end
    
end
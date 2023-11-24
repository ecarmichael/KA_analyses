%% KA screener master script. 

% load data
if ismac
    addpath(genpath('/Users/jericcarmichael/Documents/Github/vandermeerlab/code-matlab/shared'))
    addpath(genpath('/Users/jericcarmichael/Documents/Github/EC_State'));
    addpath(genpath('/Users/jericcarmichael/Documents/Github/KA_analyses'));
    data_dir = '/Users/jericcarmichael/Downloads/filtered_eric'; % where all the NLX data is.
    % inter_dir = '/Users/jericcarmichael/Dropbox/KA_Data/inter_';  % where to save the outputs.
    inter_dir = '/Users/jericcarmichael/Dropbox/KA_Data/inter_approach_2p5_new_sig';
elseif ispc
    % load data
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'))
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\EC_State'));
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\KA_analyses'));
    data_dir = 'C:\Users\ecarm\Desktop\for_eric_only'; % where all the NLX data is.
    inter_dir = 'J:\KA_Data\inter_reward_23';
    inter_dir_app = 'J:\KA_Data\inter_reward_23_approach';

end

cd(data_dir); % move to the data dir. 

% make an intermediate directory if it doesn't exist. 
if ~exist(inter_dir,'dir')
    mkdir(inter_dir)
end
    
if ~exist(inter_dir_app,'dir')
    mkdir(inter_dir_app)
end

% flagged sessions which contain some oddity like events outside of the
% recording

omit_list = {'C5_2_O7_2021-04-30_DONE',... % Feeders start way before the recording. 
    'C6_3_O1_2021-09-24_DONE',...
    'C6_3_O4_2021-09-27_DONE',...
    'C6_3_O5_2021-09-29_DONE',... 
    'C6_4_O6_2021-09-27_DONE',... % HS likly fell out. 
    };
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
    
    for iT = 1:length(cells_to_process)
        
        parts = strsplit(cells_to_process{iT}, filesep);
        this_file = parts{end};
%         This_cell = KA_screener_zscore(this_file);

%         This_cell = KA_screener_approach(this_file); 
%         This_cell = KA_screener(this_file);
%         This_cell = KA_screener_feeder(this_file);
        This_cell = KA_screener_v2(this_file);

%         This_cell = KA_screener_pseudo_baseline(this_file); 
        
        % if there were too few spikes in the .t then skip this file. 
        if ischar(This_cell) || isempty(This_cell)
            fprintf('<strong>%s</strong>: Minimum requirments not met: %s  -   <strong>%s</strong>\n', mfilename, this_file, This_cell)
            if strcmpi(This_cell, 'too short')
                success(iS) = -1;
                continue
            end
        elseif isnumeric(This_cell)
                success(iS) = -3;
                FR(length(success)) = This_cell; 
            continue
        end
        success(iS) = 1;
        
        FR(length(success)) = length(This_cell.S.t{1})/(This_cell.pos.tvec(end) - This_cell.pos.tvec(1)); 
        
        parts = strsplit(sess_list{iS}, '_');
        This_cell.subject = [parts{1} '_' parts{2}]; % get the subject ID
        This_cell.session = parts{3}; % get the session type and number Acquisition, Criteria, Overtrain, Extinction, Reacquisition. nSessions: A ? C 3, O 7, E 1, R 3-5
        This_cell.date = parts{end}; 
        
        save([inter_dir filesep sess_list{iS} '_' this_file(1:strfind(this_file, '.')-1) '_Feeder.mat'], 'This_cell')
        
        close all
    end
    
end

% summarize the files
fprintf('<strong>Total Sessions: %2.0f\nnSucess: %2.0f (%2.2f%%)\nToo short: %2.0f (%2.2f%%)\nFR too low: %2.0f (%2.2f%%)\n</strong>',...
    length(success),sum(success==1), ((sum(success==1))/length(success))*100,sum(success==2), ((sum(success==2))/length(success))*100,...
    sum(success==3),((sum(success==3))/length(success))*100)
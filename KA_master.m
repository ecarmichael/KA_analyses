%% KA screener master script. 

% load data
addpath(genpath('/Users/jericcarmichael/Documents/Github/EC_State'));
addpath(genpath('/Users/jericcarmichael/Documents/Github/KA_analyses'));
addpath(genpath('/Users/jericcarmichael/Documents/Github/vandermeerlab/code-matlab/shared'))
cd('/Users/jericcarmichael/Dropbox/C4_3_R1_2021-03-12 copy')


%% loop over sessions / cells


% example
cells_to_process = FindFiles('*.t');

for iT = 1:length(cells_to_process)
    
   parts = strsplit(cells_to_process{iT}, filesep);
   this_file = parts{end}; 
   This_cell = KA_screener(this_file);  
    
end
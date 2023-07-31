function KA_screen_peth(data_dir, out_dir, fname, type)
%% KA_screen_cells: loop over the output figures from KA_screener in a dir
% and collect them in some tile plots for easier screening. 
%
%
%  Inputs:
%       data_dir: [string]  path to the intermediate files from
%       KA_screener
%       fname: [string]  what to look for in the processed files. would be 'Feeder' or 'entry'...
%       type: [string]   suffix for the file type to copy. default 'fig'. 
%% plot basics setup. 

if nargin < 4
    type = 'fig';
end

%% get the paths for all the rate plots. 

f_list = dir([data_dir filesep '**' filesep '*' fname '.' type]);

for iF = length(f_list) :-1:1
   if contains(f_list(iF).name, fname) 
       keep_idx(iF) = 1;
   else
       keep_idx(iF) = 0;
   end
end

f_list(~keep_idx) = []; 

for iF = 1:length(f_list)
    sess_name = strsplit(f_list(iF).folder, filesep); 
    sess_name = sess_name{end}(1:strfind(sess_name{end}, 'DONE')-2); 
    
    copyfile([f_list(iF).folder filesep f_list(iF).name], [out_dir filesep   strrep(sess_name, '-', '_') '_'  f_list(iF).name(1:6) '.' type])
    
end
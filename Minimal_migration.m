function Minimal_migration

%% pick folders

% OG directory. 
fprintf('<strong>Pick original data folder...')
data_dir = uigetdir(cd, 'Original Data Directory'); 

fprintf('<strong>Original data directory: %s</strong>\n', data_dir)


% new directory. 
fprintf('<strong>Pick the folder to copy everything to...')
new_dir = uigetdir(cd, 'New Data Directory');

fprintf('\n<strong>Copying data to: %s</strong>\n', new_dir)


%% start looping over files to copy. 
this_dir = dir(data_dir); 

sess_list = [];
for ii = 1:length(this_dir)
    if strcmp(this_dir(ii).name(1), '.') % check for hidden dirs 
        continue
    else
        sess_list{ii} = this_dir(ii).name;
    end
end
sess_list =   sess_list(~cellfun('isempty',sess_list));

% loop over session folders
for iS = 1:length(sess_list)
    
    cd([data_dir filesep sess_list{iS}]); % move to the session folder. 
    
    this_new_dir = [new_dir filesep sess_list{iS}]; % copy the folder name.
    mkdir(this_new_dir); 
    
    % get the .t files in this session dir
    t_list = dir('*.t'); 
    
    for iT = 1:length(t_list)
        copyfile(t_list(iT).name, [this_new_dir filesep t_list(iT).name])
    end
    
    % get the .nev files in this session dir
    e_list = dir('*.nev'); 
    copyfile(e_list.name, [this_new_dir filesep e_list.name])
    
    
    % get the .nvt files in this session dir
    v_list = dir('*.nvt'); 
    copyfile(v_list.name, [this_new_dir filesep v_list.name])
    
    % copy the matlab trial control output
    PM_dir = dir('PM*.mat');
    copyfile(PM_dir.name, [this_new_dir filesep PM_dir.name])

    
end

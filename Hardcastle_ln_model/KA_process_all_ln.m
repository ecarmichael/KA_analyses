%% KA_process_all_ln

ln_dir = 'C:\Users\ecar\Williams Lab Dropbox\Eric Carmichael\KA_Data\Ln_model_out';


cd(ln_dir)

load('Cell_id.mat')

ln_list = dir(fullfile(ln_dir, '*.mat'));

%% loop over cell models and see what the results were

for iM = 1:length(ln_list)




end



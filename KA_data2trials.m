function KA_trl_mat(cfg_in, data)
%% takes the KA data and splits it into trials based on time from reward. 


cfg_def = []; 
cfg_def.trl_time = -4; 


cfg = ProcessConfig(cfg_def, cfg_in)

%% restrict the data to the common number of samples for the trial time; 

trl_idx = cfg.trl_time./data.pos.tvec



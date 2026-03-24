function resp_profile = KA_peth_profile(cfg_in, tvec, peth_in)
%% KA_resp_profile: takes a PETH and then slides the windows of interest to determine maximal responses






%% initialize
cfg_def = [];

cfg_def.win_s = [-1 1];

cfg = ProcessConfig(cfg_def, cfg_in);


%% making a sliding response window to find any minima/maxima
win = mode(diff(tvec))/abs(diff(cfg.win_))
for ii = 1:size(peth_in,2)
    
mov_resp = movmean(peth_in(:,ii), )

end




function [out] = KA_spd_mod(cfg_in, S, velo); 
%% KA_spd_mod: computes the speed modulation for a given cell using Kropff et al. 2015 methods
%
%
%  Inputs
%    - cfg_in: [struct]    contains configuration parameters 
%
%    - S [struct]       contains spike times in the TS format. Note will
%      only process the first cell in the structure. 
%
%    - velo: [struct] contains velocity and time data in the TSD format.
%
%% initialize

cfg_def = [];
cfg_def.binsize = 0.03;
cfg_def.gauss_window = .25./cfg_def.binsize;
cfg_def.gauss_sd = 0.025./cfg_def.binsize;

cfg = ProcessConfig(cfg_def, cfg_in); 
%% smooth data

% set up gau kernal
gauss_window = cfg.gauss_window./cfg.binsize; % 1 second window
gauss_SD = cfg.gauss_sd./cfg.binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); 

tbin_edges = velo.tvec(1):cfg.binsize:velo.tvec(end);

spk_count = histc(S.t{1},tbin_edges);
S_gau_sdf = conv2(spk_count(1:end-1),gk,'same'); % convolve with gaussian window


% spk_count = histcounts(out.S.t{1},tbin_edges);
% S_gau_sdf =smooth(spk_count, 10); 

FR_velo_int = interp1(tbin_edges(1:end-1), S_gau_sdf, velo.tvec); 

nan_idx = isnan(FR_velo_int); % points to exclude

move_idx = velo.data > 2; 


% get the pearson product
spd_corr  = corrcoef(FR_velo_int(~nan_idx & move_idx), velo.data(~nan_idx & move_idx)); 
out.spd_corr = spd_corr(1,2); 


%% compare to shufflen
nShuf = 1000; 
Shuff_shifts = randperm(length(velo.data(~nan_idx)),nShuf); 
for ii = nShuf:-1:1
   this_shuff = corrcoef(FR_velo_int(~nan_idx & move_idx), circshift(velo.data(~nan_idx & move_idx), Shuff_shifts(ii))); 
    
   spd_shuff(ii)  = this_shuff(1,2); 
end
out.spd_shuff = prctile(spd_shuff,99);
out.spd_shuff_99_abs = prctile(abs(spd_shuff),99); 

if abs(out.spd_corr) > abs(out.spd_shuff_99_abs)
            fprintf('<strong>%s</strong> Cell:  %s shows a significant speed modulation <strong>%0.2f</strong>\n',...
            mfilename, S.label{1}, out.spd_corr);
        out.spd_mod = 1; 
else
    out.spd_mod = 0; 
end

% zscore spd mod
out.z_mod = (out.spd_corr - mean(spd_shuff))./std(spd_shuff); 

out.p_val = sum(spd_shuff > out.spd_corr,2)/nShuf; 

out.FR_velo_int = FR_velo_int./cfg.binsize; 

%% for plotting later on

% plot(velo.tvec, velo.data(1,:), 'color', [.4 .4 .4], 'linewidth', 1.5)
% yyaxis('right')
% plot(velo.tvec, FR_velo_int./cfg.binsize, 'linewidth', 1.5)
% ylabel('Firing rate (Hz)')
% yyaxis('left')
% ylabel('Speed (cm/s)')
% 
% % simple acceleration
% plot(velo.tvec(1:end-1), diff(velo.data(1,:)), 'color', [.4 .4 .4], 'linewidth', 1.5)
% yyaxis('right')
% plot(velo.tvec(1:end-1), FR_velo_int(1:end-1)./cfg.binsize, 'linewidth', 1.5)
% ylabel('Firing rate (Hz)')
% yyaxis('left')
% ylabel('Speed (cm/s)')



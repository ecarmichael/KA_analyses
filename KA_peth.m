function [peth] = KA_peth(cfg_in, S, IV); 
%% KA_peth: creates a simple histogram of spikes across some time intervals. 
%
%
%  Follows methods from Fraser et al. 2023 (https://www.biorxiv.org/content/10.1101/2023.06.28.546936v1.full)
% 
% 'Peri-stimulus time histograms (PSTHs) were constructed around
% event-related responses using 0.01 ms bins. The spiking activity of each
% neuron across these bins of the PSTH was smoothed using a half-normal
% filter (σ = 6.6) that used activity in previous, but not upcoming, bins.
% To visualize the normalized activity of neurons, the mean activity within
% each of the smoothed bins of the PSTH was transformed to a z-score as
% follows: (Fi – Fmean)/ FSD, where Fi is the firing rate of the ith bin of
% the PSTH, and Fmean and FSD are the mean and SD of the firing rate of the
% 10 s baseline period. Color-coded maps and average traces of individual
% neurons’ activity were constructed based on these z-scores.'
%
%
%  Inputs:
%     -cfg_in: [struct]   contains user configuration parameters. will
%     override default fields. 
%
%     - S: [struct]    spike data in the TS format. Will only process one
%     cell in index 1. 
%
%     - IV: [struct]  contains intervals in the IV format. 
%
%
%  Outputs:
%     -peth: [struct] contains the histogram of spikes along with the cfgs.
%     



%% initialize

cfg_def = [];
cfg_def.bin = 0.01; % assuming seconds not '0.01ms' as per Fraiser. 
cfg_def.b_line = -10; % 10s baseline. 
cfg_def.sd = 6.6; % as per Fraiser. 

cfg_def.plot = 0; % toggle for plots. 


cfg = ProcessConfig(cfg_def, cfg_in); 

%% create the histogram


%% apply filtering



%% plotting if needed. 

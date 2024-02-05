function [peth] = KA_spike_hist(cfg_in, S, t, ind);
%% KA_Spike_hist: creates a simple histogram of spikes across some time intervals.
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
cfg_def.binsize = 0.05; % assuming seconds not '0.01ms' as per Fraiser.
cfg_def.b_line = -10; % 10s baseline.
cfg_def.sd = 6.6; % as per Fraiser.
cfg_def.gauss_window = 1; 
cfg_def.gauss_sd = 0.02; 
cfg_def.plot = 0; % toggle for plots.



cfg = ProcessConfig(cfg_def, cfg_in);

%% convert times (t) to IV

% IV = iv(t + cfg.win(1), t+cfg.win(2));

%% Align spikes to trials.

% HN_k = makedist('HalfNormal','mu',0,'sigma',cfg.sd/2); 

gauss_window = cfg.gauss_window./cfg.binsize; % 1 second window
gauss_SD = cfg.gauss_sd./cfg.binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); 

%
hn_gk = gk;
[~, m_idx] = max(gk); 
hn_gk(1:m_idx) = 0; 
hn_gk = (hn_gk*2)./cfg.binsize; 
hn_gk(1:m_idx) = []; 
gk = gk./cfg.binsize; % normalize by binsize

S_out = []; 
for iT = length(t):-1:1
    S0 = restrict(S, t(iT)+cfg.window(1), t(iT)+cfg.window(2));


    tbin_edges = t(iT)+cfg.window(1):cfg.binsize:t(iT)+cfg.window(2);
    tbin_centers = tbin_edges(1:end-1)+cfg.binsize/2;


    % get the histogram around the events.

    [spk_count, B] = histcounts(S0.t{1},tbin_edges);
    spk_count =  spk_count(1:end-1);

        S_gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window

    S_gau_out(iT,:) = S_gau_sdf; 
    S_out(iT,:) = spk_count; % save the counts.
end


S_rate = sum(S_out,1)./(tbin_edges(end) - tbin_edges(1));

S_gau = sum(S_gau_out);


peth.S_gau = mean(S_gau_out)/cfg.binsize;
peth.tbins = tbin_edges(1:end-2); 
%% apply filtering


peth = [];
%% plotting if needed.
figure(999)
clf

subplot(2,1,1);
cla
plot(tbin_edges(1:end-2), mean(S_gau_out)/cfg.binsize)
hold on
% plot(tbin_edges(1:end-2), conv2(mean(S_out), gk, 'same'), 'b')

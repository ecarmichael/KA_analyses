%% KA_ln_modeling_CeA


% This is essentially a wrapper for the hardcastle et al. 2017 pipeline. 


% This script is segmented into several parts. First, the data (an
% example cell) is loaded. Then, 15 LN models are fit to the
% cell's spike train. Each model uses information about 
% position, head direction, running speed, theta phase,
% or some combination thereof, to predict a section of the
% spike train. Model fitting and model performance is computed through
% 10-fold cross-validation, and the minimization procedure is carried out
% through fminunc. Next, a forward-search procedure is
% implemented to find the simplest 'best' model describing this spike
% train. Following this, the firing rate tuning curves are computed, and
% these - along with the model-derived response profiles and the model
% performance and results of the selection procedure are plotted.

% Code as implemented in Hardcastle, Maheswaranthan, Ganguli, Giocomo,
% Neuron 2017
% V1: Kiah Hardcastle, March 16, 2017


%% Clear the workspace and load the data

clear all; close all; clc

% load the data
fprintf('(1/5) Loading data from example cell \n')
load data_for_cell77

load('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\KA_Data\inter_reward_23\C3_4_C3_2020-08-31_DONE_maze_data.mat')

iC = 1;
        % isolate the cell of interest in the session (if there are
        this_S = KA_isolate_S(data.S, data.S.label{iC});
        % speed modulation ( add spd and acc MI later)
        % data.spd_mod = KA_spd_mod([], this_S, data.velo_smooth);
        % data.FR{iC} = data.spd_mod.FR_velo_int; 

% firing rate to match data samples
dt = mode(diff(data.pos.tvec)); 
tbin_edges = data.pos.tvec(1):dt:data.pos.tvec(end);

spk_count = histc(this_S.t{1},tbin_edges);
FR_velo_int = interp1(tbin_edges(1:end-1), spk_count(1:end-1), data.pos.tvec); 


% compute a filter, which will be used to smooth the firing rate
filter = gaussmf(-4:4,[2 0]); filter = filter/sum(filter); 
 fr = FR_velo_int/dt;
smooth_fr = conv(fr,filter,'same');


% convert the time into time from events (rewards and approach time)
t_minus_r = NaN(size(data.pos.tvec)); 

for ii = length(data.rew.t):-1:1
    if ii == 1
        this_idx = nearest_idx3(data.rew.t(ii), data.pos.tvec);
        prior_idx = 1;
    else
        this_idx = nearest_idx3(data.rew.t(ii), data.pos.tvec);
        prior_idx = nearest_idx3(data.rew.t(ii-1), data.pos.tvec);
    end

    t_minus_r(prior_idx:this_idx) = data.pos.tvec(prior_idx:this_idx) - data.rew.t(ii); 
end

post = data.pos.tvec; 
spiketrain = FR_velo_int; 
posx_c = data.pos.data(1,:)'; 
posy_c = data.pos.data(2,:)'; 
boxSize = 150-20; 

mov_idx = data.velo_smooth.data > 5; % keep data when moving. 

% description of variables included:
% boxSize = length (in cm) of one side of the square box
% post = vector of time (seconds) at every 20 ms time bin
% spiketrain = vector of the # of spikes in each 20 ms time bin
% posx = x-position of left LED every 20 ms
% posx2 = x-position of right LED every 20 ms
% posx_c = x-position in middle of LEDs
% posy = y-position of left LED every 20 ms
% posy2 = y-posiiton of right LED every 20 ms
% posy_c = y-position in middle of LEDs
% filt_eeg = local field potential, filtered for theta frequency (4-12 Hz)
% eeg_sample_rate = sample rate of filt_eeg (250 Hz)
% sampleRate = sampling rate of neural data and behavioral variable (50Hz)

%% fit the model
fprintf('(2/5) Fitting all linear-nonlinear (LN) models\n')
fit_all_ln_models

%% find the simplest model that best describes the spike train
fprintf('(3/5) Performing forward model selection\n')
select_best_model

%% Compute the firing-rate tuning curves
fprintf('(4/5) Computing tuning curves\n')
compute_all_tuning_curves

%% plot the results
fprintf('(5/5) Plotting performance and parameters\n') 
plot_performance_and_parameters

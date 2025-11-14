function KA_run_ln(data, iC, fname, save_dir)


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


%% load the data
sname = strsplit(fname, '_DONE'); 

fprintf('(1/5) Loading data from %s cell #%d \n', sname{1},iC)
% load data_for_cell77

% isolate the cell of interest in the session (if there are
this_S = KA_isolate_S(data.S, data.S.label{iC});


% firing rate to match data samples
dt = mode(diff(data.pos.tvec)); 
tbin_edges = data.pos.tvec;

spk_count = histc(this_S.t{1},tbin_edges);
FR_velo_int = interp1(tbin_edges(1:end-1), spk_count(1:end-1), data.pos.tvec); 
 

% % compute a filter, which will be used to smooth the firing rate
filter = gaussmf(-4:4,[2 0]); filter = filter/sum(filter); 
 fr = FR_velo_int/dt;
smooth_fr = conv(fr,filter,'same');


% convert the time into time from events (rewards and approach time)
t_minus_r = NaN(size(data.pos.tvec)); 

for ii = length(data.rew.t):-1:1
    if ii == 1
        this_idx = nearest_idx(data.rew.t(ii), data.pos.tvec);
        prior_idx = 1;
    else
        this_idx = nearest_idx(data.rew.t(ii), data.pos.tvec);
        prior_idx = nearest_idx(data.rew.t(ii-1), data.pos.tvec);
    end

    t_minus_r(prior_idx:this_idx) = data.pos.tvec(prior_idx:this_idx) - data.rew.t(ii); 
end

t_minus_r = -t_minus_r; %make positive

hd_data = data.pos.data(3,:); 
velo_smooth = data.velo_smooth.data(1,:); 
post = data.pos.tvec; 
spiketrain = FR_velo_int'; 
posx_c = data.pos.data(1,:)'; 
posy_c = data.pos.data(2,:)'; 
boxSize = 150-20; 

spiketrain = fillmissing(spiketrain, "nearest"); 

% quick plot for sanity
figure(1000+iC)
hold on
plot(post, posx_c);
plot(post, t_minus_r);
plot(post, velo_smooth);
yyaxis right 
plot(post, spiketrain);

close(gcf)

% remove nans from the signal. 
% nan_idx = isnan(spiketrain); 
% 
% spiketrain(nan_idx) = []; 
% post(nan_idx) = []; 
% velo_smooth(nan_idx) = []; 
% posx_c(nan_idx) = []; 
% posy_c(nan_idx) = []; 
% hd_data(nan_idx) = []; 
% t_minus_r(nan_idx) = []; 

% remove periods that were too far removed from a reward. 
mov_idx = data.velo_smooth.data > 5; % keep data when moving. 
ITI_idx = t_minus_r >= 10;

rm_idx = ITI_idx | ~mov_idx; % remove periods of immobility too far from next reward. 

clear data
%% fit the model
fprintf('(2/5) Fitting all linear-nonlinear (LN) models\n')
fit_all_ln_models_KA

%% find the simplest model that best describes the spike train
fprintf('(3/5) Performing forward model selection\n')
select_best_model_KA

%% Compute the firing-rate tuning curves
fprintf('(4/5) Computing tuning curves\n')

compute_all_tuning_curves_KA

%% save the workspace for now. 

save([save_dir filesep sname{1} '_' num2str(iC) '.mat'])
%% plot the results
fprintf('(5/5) Plotting performance and parameters\n') 
% plot_performance_and_parameters

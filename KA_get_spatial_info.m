function [spd_metrics, t_metrics, p_metrics] = KA_get_spatial_info(data, iC)


spd_metrics = []; 
t_metrics = []; 
p_metrics = []; 

this_S = KA_isolate_S(data.S, data.S.label{iC});


dt = mode(diff(data.pos.tvec)); 

trials = [nearest_idx3(data.rew.t-5, data.pos.tvec), nearest_idx3(data.rew.t+2.5, data.pos.tvec)]; 

data.velo_smooth.data(data.velo_smooth.data>50) = NaN; 
% figure
[~, ~, spd_metrics ] =spatial_information_KA(data.velo_smooth.data, data.velo_smooth.tvec, this_S.t, trials, 25, 8, 1000);


% time to reward
t_minus_r = NaN(size(data.pos.tvec)); 
for ii = length(data.rew.t):-1:1
    if ii == 1
        this_idx = nearest_idx2(data.rew.t(ii)+2.5, data.pos.tvec);
        prior_idx = 1;
    else
        this_idx = nearest_idx2(data.rew.t(ii)+2.5, data.pos.tvec);
        prior_idx = nearest_idx2(data.rew.t(ii-1), data.pos.tvec);
    end
    t_minus_r(prior_idx:this_idx) = data.pos.tvec(prior_idx:this_idx) - data.rew.t(ii); 
end

t_minus_r = -t_minus_r; %make positive
t_minus_r(t_minus_r > 5) = NaN; 
% figure
[~, ~, t_metrics ] =spatial_information_KA(t_minus_r, data.velo_smooth.tvec, this_S.t, trials, 8, 8, 1000);


% figure

[~, ~, p_metrics ] =spatial_information_KA(rand(size(data.velo_smooth.tvec)), data.velo_smooth.tvec, this_S.t, trials, 8, 8, 1000);

%% 1D vectors:

% % speed
% nbins = 20;
% smax = 50; 
% speedVec = 0:smax/nbins:smax;
% % speed_grid = zeros(numel(data.velo_smooth.data),numel(speedVec));
% 
% [spd.MI, spd.posterior, spd.occupancy_vector, spd.p_active, spd.likelihood] = MS_get_spatial_information(binary_in, position_in, bin_vec)





% KA_SI(A_vec, B_vec, bin_size)


%% 

%% place
% 
% cfg.X_bins = 0:cfg.p_bin_size:ceil(max(behav_aligned.position(:,1)));
% X_bin_centers = cfg.X_bins +  cfg.p_bin_size/2;
% X_bin_centers = X_bin_centers(1:end-1);
% % same for Y bins
% cfg.Y_bins = 0:cfg.p_bin_size:ceil(max(behav_aligned.position(:,2)));
% Y_bin_centers = cfg.Y_bins +  cfg.p_bin_size/2;
% Y_bin_centers = Y_bin_centers(1:end-1);
% 
% %% Speed
% 
% % make speed bins
% Speed_bin_centers = 0 + cfg.s_bin_size/2;
% Speed_bin_centers = Speed_bin_centers(1:end-1);
% 
% % acceleration bins
% Acc_bin_centers = cfg.a_bins + cfg.a_bin_size/2;
% Acc_bin_centers = Acc_bin_centers(1:end-1);
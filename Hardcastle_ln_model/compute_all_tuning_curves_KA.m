% function [pos_curve, hd_curve, speed_curve, theta_curve]  = compute_all_tuning_curves(smooth_fr, posx_c, posy_c, direction, speed, phase, too_fast)

%% Description
% This will compute the firing rate tuning curves for position, head
% direction, running speed, and theta phase.

% take out times when the animal ran >= 50 cm/s
posx_c(rm_idx) = []; posy_c(rm_idx) = []; 
hd_data(rm_idx) = [];
velo_smooth(rm_idx) = [];
t_minus_r(rm_idx) = [];

% compute tuning curves for position, head direction, speed, and theta phase
[pos_curve] = compute_2d_tuning_curve(posx_c,posy_c,smooth_fr,n_pos_bins,0,boxSize);
[hd_curve] = compute_1d_tuning_curve(hd_data,smooth_fr,n_dir_bins,0,2*pi);
[speed_curve] = compute_1d_tuning_curve(velo_smooth,smooth_fr,n_speed_bins,0,50);
[t_minus_curve] = compute_1d_tuning_curve(abs(t_minus_r),fr,n_tminus_bins,0,20);
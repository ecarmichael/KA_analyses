function [MI, posterior, occ_vec, p_active, likelihood] = KA_SI(S, vec, bin_size)
%%KA_SI: computes the mutual information between spike times and a vector. 

% based on the original by GE from: https://github.com/etterguillaume/CaImDecoding/blob/master/extract_1D_information.m
%MSPLACE_CELL_OF_POSITIONS Analyse and plot spatial information
%   This function plots in space the position related to calcium

% Copyright (C) 2017-2019 by Guillaume Etter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or any
% later version.
% Contact: etterguillaume@gmail.com

% binary_trace: logical vector representing neurons active/inactive periods
%
% interp_behav_vec: m x 2  matrix representing behavioral state in two dimensions (eg position in a
% maze). Has to be the same size as binary_trace
%
% ca_time: time vector for both binary_trace and interp_behav_trace. If time correspondance
% is different between these two variables, you can use the fonction
% interp_behav to interpolate the behavioral trace
%
% OPTIONAL vectors:
% inclusion_vec: logical vector including (1) or excluding (0) corresponding timestamps. Has to be the same size as binary_trace.


%% set up the spike times as a vector. 

disp(S); 
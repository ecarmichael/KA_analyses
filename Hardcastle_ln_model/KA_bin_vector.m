[pos_mat, spd_mat, tminus_mat] = KA_bin_vector_data(data); 



post = data.pos.tvec; 
posx_c = data.pos.data(1,:)'; 
posy_c = data.pos.data(2,:)'; 


% extract the event timing and make a vector
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


n_pos_bins = 20;
% n_dir_bins = 18;
n_speed_bins = 10;
% n_theta_bins = 18;
n_tminus_bins = 20; 

% compute position matrix
[posgrid, posVec] = pos_map([data.pos.data(1,:)' data.pos.data(2,:)'], n_pos_bins, boxSize);

% compute head direction matrix
% [hdgrid,hdVec,direction] = hd_map(posx,posx2,posy,posy2,n_dir_bins);

% compute speed matrix
[speedgrid,speedVec] = velo_map(data.velo_smooth.data ,n_speed_bins);

% compute theta matrix
[tminus_grid,tminusVec] = tminus_map(-t_minus_r,n_tminus_bins);

% % remove times when the animal ran > 50 cm/s (these data points may contain artifacts)
% ITI_idx = find(tminus_grid >= 20);
% 
% posgrid(ITI_idx,:) = []; 
% speedgrid(ITI_idx,:) = []; 
% tminus_grid(ITI_idx,:) = [];
% spiketrain(ITI_idx) = [];
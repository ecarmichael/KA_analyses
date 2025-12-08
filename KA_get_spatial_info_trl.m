function [spd_metrics, t_metrics, p_metrics] = KA_get_spatial_info_trl(data, iC)


spd_metrics = []; 
t_metrics = []; 
p_metrics = []; 

this_S = KA_isolate_S(data.S, data.S.label{iC});

dt = mode(diff(data.pos.tvec)); 
data.velo_smooth.data(data.velo_smooth.data>50) = NaN; 
data.velo_smooth.data(data.velo_smooth.data<2) = NaN; 


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

%% loop over trial types.
t_types = 1:4;
for tt = length(t_types):-1:1

    t_idx = data.rew.in == t_types(tt);

    trials = [nearest_idx3(data.rew.t(t_idx)-5, data.pos.tvec), nearest_idx3(data.rew.t(t_idx)+2.5, data.pos.tvec)];
    trials_iv = iv(data.rew.t(t_idx)-5, data.rew.t(t_idx)+2.5); 
    S_trl = restrict(this_S, trials_iv); 


    if isempty(trials) || isempty (S_trl.t{1})
        spd_metrics{tt}.MI = NaN; 
        spd_metrics{tt}.Isec = NaN; 
        spd_metrics{tt}.Ispike = NaN; 
        spd_metrics{tt}.NormMI = NaN; 
        spd_metrics{tt}.NormIsec = NaN; 
        spd_metrics{tt}.NormIspike = NaN; 
        t_metrics{tt} = spd_metrics{tt}; 
        p_metrics{tt} = spd_metrics{tt};
    else

        [~, ~, spd_metrics{tt}] =spatial_information_KA(data.velo_smooth.data, data.velo_smooth.tvec, this_S.t, trials, 25, 8, 1000);

        [~, ~, t_metrics{tt} ] =spatial_information_KA(t_minus_r, data.velo_smooth.tvec, this_S.t, trials, 8, 8, 1000);

        [~, ~, p_metrics{tt}] =spatial_information_KA(rand(size(data.velo_smooth.tvec)), data.velo_smooth.tvec, this_S.t, trials, 8, 8, 1000);

    end

end


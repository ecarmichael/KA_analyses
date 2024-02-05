function wave_props = KA_Cell_stats(data_in, pos, plot_flag)
%% KA_Cell_stats: gathers basic spiking stats for a cell.
%
%
%
%  Inputs
%    - data_in [struct] contains spike 'S', and tracking 'pos' data in the
%    TSD and TS formats as well as other .
%
%    - cell_id: [string] which S to use (must match the label field of
%    data_in.S.label
%
%
%
% Outputs
%    - spike_stats: [struct] contains basic spike measures ISI, mean FR,
%

%% init
if nargin < 3
    plot_flag = 0; 
end


%%

% if ~ismember(data_in.S.label, cell_id)
%     error('cell_id does not match data_in.S.label')
% end
%
% % isolate the cell from the sesson data.
% idx = (ismember(data_in.S.label, cell_id));
% this_S = data_in.S;
% this_S.t(~idx) = [];
% this_S.label(~idx) = [];
% this_S.waves(~idx) = [];
% this_S.Q(~idx) = [];

%%
% s_stats.FR = length(data_in.t{1})/(pos.tvec(end) - pos.tvec(1));
% s_stats.ISI = mode(diff(data_in.t{1}))*1000;
% [s_stats.ISI_hist.n, s_stats.ISI_hist.b] = hist(diff(data_in.t{1}*1000), 0:1:500);
% 
rng(123, "twister"); 
for ii = 1000:-1:1
    t0 = randsample(pos.tvec(100:end-100), 1, true); 

    S0 = restrict(data_in, t0-2.5, t0+2.5); 
    fr_d(ii) = length(S0.t{1})/((t0+2.5) - (t0-2.5)); 

end


 [~, mu, z_std] = zscore(fr_d); 



%% waveform classification
wave_time = (0:31)/32000; %NLX wave samples. 
wave_props = MS_get_wave_properties(data_in, [wave_time(1:32)', data_in.waves{1}.mWV]', pos.tvec, plot_flag ); 

wave_props.std = z_std; 




end
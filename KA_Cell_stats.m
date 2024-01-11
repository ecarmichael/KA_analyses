function stats = KA_Cell_stats(data_in); 
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
s_stats.FR = length(data_in.S.t{1})/(data_in.pos.tvec(end) - data_in.pos.tvec(1)); 
s_stats.ISI = mode(diff(data_in.S.t{1}))*1000; 
[s_stats.ISI_hist.n, s_stats.ISI_hist.b] = hist(diff(data_in.S.t{1}*1000), 0:1:500);







end
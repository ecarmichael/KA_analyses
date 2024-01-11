function [stats] = KA_process(cfg_in, inter_dir)
%% KA_process: take the output from KA master and extract cell response metrics
%
%
%
%    Inputs: 
%    - cfg [struct]   configuration see the defaults below. 
%
%    - 
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2023-11-10   initial version 
%
%
%
%% initialize
cfg_def = [];
cfg_def.min_resp = .3; 


cfg = ProcessConfig(cfg_def, cfg_in);
%%  loop over cells and determine which ones to keep. 
cd(inter_dir)
sess_list = dir('*.mat');

keep_idx = [];nR = []; nT = []; 
for ii =1:length(sess_list)%:-1:1
   
    load(sess_list(ii).name)
    
    % check the number of trials with a spike
    
    trials = This_cell.Zone_ids;
    
    nS = []; nT = []; 
    for iT = 1:4
        this_trial = This_cell.Zone_times(trials == iT); 
        
        % remove
        nT(iT) = length(this_trial); 
    end
    
    fprintf('%s %s %s N: %0.0f%% | W: %0.0f%% | S: %0.0f%% | E: %0.0f%% \n', This_cell.subject, This_cell.session, This_cell.S.label{1}, nT(1), nT(2), nT(3), nT(4))

    %end remove
        
%         for jj = length(this_trial):-1:1
%            s_t = restrict(This_cell.S, this_trial(jj) + This_cell.cfg_peth.window(1), this_trial(jj) + This_cell.cfg_peth.window(2));
%             
%            nS{iT}(jj) = length(s_t.t{1}); 
%           
%         end
%         nT(ii,iT) = length(this_trial); 
%         nR(ii,iT) = sum(nS{iT} >0) ./length(nS{iT}); 
%     end
%     
%     fprintf('%s %s %s N: %0.0f%% | W: %0.0f%% | S: %0.0f%% | E: %0.0f%% \n', This_cell.subject, This_cell.session, This_cell.S.label{1}, nR(ii,1)*100, nR(ii,2)*100, nR(ii,3)*100, nR(ii,4)*100)
%     
    
%     sum(sum(This_cell.outputS{end},2) == 0)/size(This_cell.outputGau{end},2)
    
% wave_prop{ii} = MS_get_wave_properties(This_cell.S, [(0:length(This_cell.cfg_peth.waves.mWV)-1)./1/32000; This_cell.cfg_peth.waves.mWV'], This_cell.pos.tvec, 1);
    


end

%% make a plot of the firing per 

%% waveform class

for ii = length(sess_list):-1:1
    
    
    fr(ii) = wave_prop{ii}.firing_rate; 
    isi(ii) = mean(wave_prop{ii}.ISI);
    pt(ii) = wave_prop{ii}.pt_ratio; 
    slp(ii) = wave_prop{ii}.slopes_ratio; 
    width(ii) = wave_prop{ii}.spike_width ; 
    risefall(ii) = wave_prop{ii}.rise_fall_inter; 
end

data_in = [fr', isi', width'];

[g_idx, n_idx, frq_idx]= MS_kmean_scatter(data_in, 3,[3,2,1], 50);
clus1 = frq_idx(1); % get the top two clusters. 
clus2 = frq_idx(2); 

close(gcf);

figure(221)
subplot(2,3,[1,2,4,5])
hold on
col_ord = linspecer(2); 
scatter3(width(g_idx==clus1), isi(g_idx==clus1), fr(g_idx==clus1), 50,col_ord(1,:), 'filled');
scatter3(width(g_idx==clus2), isi(g_idx==clus2), fr(g_idx==clus2), 50,col_ord(2,:), 'filled');

grid on
xlabel('width'); ylabel('isi'); zlabel('Firing rate (Hz)')
legend('Cluster 1', 'Cluster 2')



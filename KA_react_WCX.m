function [p, h, all_base_fr, all_trial_fr] = KA_react_WCX(cfg_in, S, t, id)
%% KA_react_WCX: Peri-event response using a t statistic from Fraser et al. 2023 (https://www.biorxiv.org/content/10.1101/2023.06.28.546936v1.full)
% Method:
%
% neurons were determined to be modulated by an event if the spike rate in a
% custom window (−0.5 to 0.5 s for port entries and port exits and 0 to
% 0.03 s for lick) following each event significantly differed from a 10 s
% baseline period according to a Wilcoxon signed-rank test (p < 0.05,
% two-tailed).
% 
% ...For electrophysiological data, statistical tests were performed on
% unsmoothed data. The specific tests performed are noted throughout the
% text and figure legends. For electrophysiological data we did not test for
% normality, but made use of nonparametric tests (two-sided Wilcoxon’s
% rank-sum and signed-rank tests)



%

%% initialize

cfg_def = [];

cfg_def.baseline = [-10 0];
cfg_def.win = [-1 1];
cfg_def.alpha = 0.05;
cfg_def.label = {'North', 'West', 'South', 'East', 'All'}; 

%ploting
cfg_def.bin = 0.01; % says "using 0.01ms bins" but this seems very small.

cfg = ProcessConfig(cfg_def, cfg_in);



%%

for ii = unique(id)
    
    
    F_idx = find(id == ii);
    base_fr = [];
    % baseline per trial.
    for jj = length(F_idx):-1:1
        base_S = restrict(S, t(F_idx(jj))+cfg.baseline(1), t(F_idx(jj))+cfg.baseline(2));
        
        if isempty(base_S.t{1})
            base_fr(jj) = 0;
        else
            base_fr(jj) = length(base_S.t{1})./ ((t(F_idx(jj))+cfg.baseline(2)) - (t(F_idx(jj))+cfg.baseline(1)));
        end
    end
    
    
    % reward
    % same but wilcoxen
    trial_fr = [];
    for jj = length(F_idx):-1:1
        trial_S = restrict(S, t(F_idx(jj)) + cfg_def.win(1), t(F_idx(jj)) + cfg_def.win(2));
        if isempty(trial_S.t{1})
            trial_fr(jj) = 0;
        else
            trial_fr(jj) = length(trial_S.t{1})./ abs(diff(cfg_def.win));
        end
    end
    
    [p(ii), h(ii)] = signrank(base_fr, trial_fr, 'tail', 'both', 'alpha', 0.05);
    if h(ii)
    fprintf('<strong>%s</strong> Cell:  %s has significantly activity change at %0.2f : %0.2f s (%0.2fHz) vs baseline  (%0.2fHz | %0.2f : %0.2fs) <strong>(p = %0.5f) at  %s </strong> arm \n',...
        mfilename, S.label{1}, cfg.win(1), cfg.win(2),mean(trial_fr), mean(base_fr),cfg.baseline(1), cfg.baseline(2), p(ii), cfg.label{id(ii)});
    end
    all_base_fr(ii) = mean(base_fr); 
    all_trial_fr(ii) = mean(trial_fr); 
end

%% same thing but all the trial types.
ii = 5;

F_idx = 1:length(t);
base_fr = [];
% baseline per trial.
for jj = length(F_idx):-1:1
    base_S = restrict(S, t(F_idx(jj))+cfg.baseline(1), t(F_idx(jj))+cfg.baseline(2));
    
    if isempty(base_S.t{1})
        base_fr(jj) = 0;
    else
        base_fr(jj) = length(base_S.t{1})./ ((t(F_idx(jj))+cfg.baseline(2)) - (t(F_idx(jj))+cfg.baseline(1)));
    end
end


% reward
% same but wilcoxen
trial_fr = [];
for jj = length(F_idx):-1:1
    trial_S = restrict(S, t(F_idx(jj)) + cfg_def.win(1), t(F_idx(jj)) + cfg_def.win(2));
    if isempty(trial_S.t{1})
        trial_fr(jj) = 0;
    else
        trial_fr(jj) = length(trial_S.t{1})./ abs(diff(cfg_def.win));
    end
end

[p(ii), h(ii)] = signrank(base_fr, trial_fr, 'tail', 'both', 'alpha', 0.05);
if h(ii)
    fprintf('<strong>%s</strong> Cell:  %s has significantly activity change at %0.2f : %0.2fs (%0.2fHz) vs baseline (%0.2fHz | %0.2f : %0.2fs) <strong>(p = %0.5f) overall </strong>\n',...
        mfilename, S.label{1}, cfg.win(1), cfg.win(2),mean(trial_fr), mean(base_fr),cfg.baseline(1), cfg.baseline(2), p(ii));
end

    all_base_fr(ii) = mean(base_fr); 
    all_trial_fr(ii) = mean(trial_fr); 


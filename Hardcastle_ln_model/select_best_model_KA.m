%% Description
% This code will implement forward feature selection in order to determine
% the simplest model that best describes neural spiking. First, the
% highest-performing single-variable model is identified. Then, the
% highest-perfmoring double-variable model that includes the
% single-variable model is identified. This continues until the full model
% is identified. Next, statistical tests are applied to see if including
% extra variables significantly improves model performance. The first time
% that including variable does NOT signficantly improve performance, the
% procedure is stopped and the model at that point is recorded as the
% selected model.

% the model indexing scheme:
% phst, phs, pht, pst, hst, ph, ps, pt, hs, ht, st, p,  h,  s,  t
% 1      2    3    4    5    6  7   8   9   10  11  12  13  14  15

% KA version
% pst, ps, pt, st, p, s, t
%  1    2   3  4   5  6  7


testFit_mat = cell2mat(testFit);
LLH_values = reshape(testFit_mat(:,3),numFolds,numModels);



% find the best single model
singleModels = 5:7;
[~,top1] = max(nanmean(LLH_values(:,singleModels))); top1 = top1 + singleModels(1)-1;

% find the best double model that includes the single model
if top1 == 5 % P -> PS, PT
    [~,top2] = max(nanmean(LLH_values(:,[2 3])));
    vec = [2 3]; top2 = vec(top2);
elseif top1 == 6 % S -> PS, PT
    [~,top2] = max(nanmean(LLH_values(:,[2 4])));
    vec = [2 4]; top2 = vec(top2);
elseif top1 == 7 % S -> PS, HS, ST
    [~,top2] = max(nanmean(LLH_values(:,[3 4])));
    vec = [3 4]; top2 = vec(top2);
% else % T -> PT, HT, ST
%     [~,top2] = max(nanmean(LLH_values(:,[8 10 11])));
%     vec = [8 10 11]; top2 = vec(top2);
end

% find the best triple model that includes the double model
if top2 == 2 % PS-> PST
    [~,top3] = max(nanmean(LLH_values(:,[1])));
    vec = [1]; top3 = vec(top3);
elseif top2 == 3 % PT -> PST
    [~,top3] = max(nanmean(LLH_values(:,[1])));
    vec = [1]; top3 = vec(top3);
elseif top2 == 4 % ST -> PST
    [~,top3] = max(nanmean(LLH_values(:,[1])));
    vec = [1]; top3 = vec(top3);
% elseif top2 == 9 % HS -> PHS, HST
%     [~,top3] = max(nanmean(LLH_values(:,[2 5])));
%     vec = [2 5]; top3 = vec(top3);
% elseif top2 == 10 % HT -> PHT, HST
%     [~,top3] = max(nanmean(LLH_values(:,[3 5])));
%     vec = [3 5]; top3 = vec(top3);
% elseif top2 == 11 % ST -> PST, HST
%     [~,top3] = max(nanmean(LLH_values(:,[4 5])));
%     vec = [4 5]; top3 = vec(top3);
end

top4 = 1;
LLH1 = LLH_values(:,top1); LLH2 = LLH_values(:,top2);
LLH3 = LLH_values(:,top3); LLH4 = LLH_values(:,top4);

[p_llh_12,~] = signrank(LLH2,LLH1,'tail','right');
[p_llh_23,~] = signrank(LLH3,LLH2,'tail','right');
[p_llh_34,~] = signrank(LLH4,LLH3,'tail','right');


if p_llh_12 < 0.05 % double model is sig. better
    if p_llh_23 < 0.05  % triple model is sig. better
        if p_llh_34 < 0.05 % full model is sig. better
            selected_model = 1; % full model
        else
            selected_model = top3; %triple model
        end
    else
        selected_model = top2; %double model
    end
else
    selected_model = top1; %single model
end

% re-set if selected model is not above baseline
pval_baseline = signrank(LLH_values(:,selected_model),[],'tail','right');

if pval_baseline > 0.05
    selected_model = NaN;
end

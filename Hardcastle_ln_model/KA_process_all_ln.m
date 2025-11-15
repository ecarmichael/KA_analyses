%% KA_process_all_ln

ln_dir = 'C:\Users\ecar\Williams Lab Dropbox\Eric Carmichael\KA_Data\Ln_model_out';


cd(ln_dir)

load('Cell_id.mat')

ln_list = dir(fullfile(ln_dir, '*.mat'));

for ii =length(ln_list):-1:1
    if contains(ln_list(ii).name, 'Cell_id.mat')
        rm_idx(ii) = true; 
    else
        rm_idx(ii) = false;
    end
end
ln_list(rm_idx) = []; 
%% loop over cell models and see what the results were
models = {'pst', 'ps', 'pt', 'st', 'p', 's', 't'}; 
for iM = length(ln_list):-1:1
    load(ln_list(iM).name);
    out{iM}.selected_model = selected_model; 
    out{iM}.pval_baseline = pval_baseline; 

    if ~isnan(selected_model)
        disp([num2str(iM) ': ' models{selected_model}]);
        figure(iM);
        plot_performance_and_parameters_KA;
        clearvars -except iM ln_list ln_dir out models;
    end
end


%% count

for iM = length(out):-1:1

    winning_model(iM) = out{iM}.selected_model; 

end


m1 = (sum(winning_model == 1) ./ length((winning_model)))*100; 
m2 = (sum(winning_model == 2) ./ length((winning_model)))*100; 
m3 = (sum(winning_model == 3) ./ length((winning_model)))*100; 
m4 = (sum(winning_model == 4) ./ length((winning_model)))*100; 
m5 = (sum(winning_model == 5) ./ length((winning_model)))*100; 
m6 = (sum(winning_model == 6) ./ length((winning_model)))*100; 
m7 = (sum(winning_model == 7) ./ length((winning_model)))*100; 


figure(1010)
modelPercentages = [m1, m2, m3, m4, m5, m6, m7];
bar(modelPercentages);
xlabel('Model');
ylabel('Percentage of Winning Models (%)');
title('Winning Model Percentages');
set(gca, 'XTickLabel', models);
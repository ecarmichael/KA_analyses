cd(inter_dir);
file_names = FindFiles('*Feeder*');  % get all the files containing 'TT' which should be all the outputs from KA_master.

Zone_names = {'North', 'West', 'South', 'East', 'All'};
subj_list = {'C1_1', 'C3_2','C3_3', 'C3_4', 'C4_3', 'C5_2', 'C5_3','C6_3', 'C6_4'};
deval_type = {'WE', 'NS', 'NS', ' ' , 'WE', 'WE', 'NS', 'WE', 'NS'};
for ii = 1:3  % Criteria should be 3
    all_out.(['C' num2str(ii)]) = cell(1,5);
end

for ii = 1:7 % Overtaining should b 7
    all_out.(['O' num2str(ii)]) = cell(1,5);
end

for ii = 1:1 % Extinction is only 1
    all_out.(['E' num2str(ii)]) = cell(1,5);
end

for ii = 1:5 % Reacqusition is up to 5
    all_out.(['R' num2str(ii)]) = cell(1,5);
end

all_rew_pre_z = all_out;
all_rew_act_z = all_out;
all_act_pre_z = all_out;


%%

all_rew_pre_z = all_out;
all_rew_act_z = all_out;
all_act_pre_z = all_out;

for iF = 1:length(file_names)
    load(file_names{iF}); % load the data as 'This_cell'
    
    this_sess = This_cell.session;
    
    for ii = 1:5
        all_rew_pre_z.(this_sess){ii} = [all_rew_pre_z.(this_sess){ii}, This_cell.rew_pre_z{ii}];
        all_rew_act_z.(this_sess){ii} = [all_rew_act_z.(this_sess){ii}, This_cell.rew_act_z{ii}];
        all_act_pre_z.(this_sess){ii} = [all_act_pre_z.(this_sess){ii}, This_cell.act_pre_z{ii}];
    end

end

%% collect across sess types
sessions = fieldnames(all_rew_pre_z);

all_rew_pre_z_mat = cell(1,5);
all_rew_act_z_mat = cell(1,5);
all_act_pre_z_mat = cell(1,5);


for iS = 1:length(sessions)
    if isempty(all_out.(sessions{iS}))
        continue
    end
    
    for ii = 1:5
        all_rew_pre_z_mat{ii} = [all_rew_pre_z_mat{ii}, all_rew_pre_z.(sessions{iS}){ii}];
        all_rew_act_z_mat{ii} = [all_rew_act_z_mat{ii}, all_rew_act_z.(sessions{iS}){ii}];
        all_act_pre_z_mat{ii} = [all_act_pre_z_mat{ii}, all_act_pre_z.(sessions{iS}){ii}];
    end
    
end

%% plot some histograms

figure(808)

n = 3; m = 5; 
c_ord = linspecer(5); 
c_ord(5,:) = [0 0 0]; 

for ii = 1:5
subplot(n, m, ii)
histogram(all_rew_pre_z_mat{ii}, 20, 'facecolor', c_ord(ii,:)); title(['Rew-pre ' Zone_names{ii}]); xlim([-3, 3]); vline([-1.96, 1.96]);

subplot(n, m, m+ii)
histogram(all_rew_act_z_mat{ii}, 20, 'facecolor', c_ord(ii,:)); title(['Rew-act ' Zone_names{ii}]);xlim([-3, 3]); vline([-1.96, 1.96]);

subplot(n, m, (m*2)+ii)
histogram(all_act_pre_z_mat{ii}, 20, 'facecolor', c_ord(ii,:)); title(['Act_pre ' Zone_names{ii}]);xlim([-3, 3]); vline([-1.96, 1.96]);

end

SetFigure([], gcf)

    
%% KA_run_all_ln
clearvars
close all


usr_name = char(java.lang.System.getProperty('user.name')); 

if ismac

    addpath(genpath(['/Users/' usr_name '/Documents/Github/vandermeerlab/code-matlab/shared']))
    addpath(genpath(['/Users/' usr_name '/Documents/Github/EC_State']));
    addpath(genpath(['/Users/' usr_name '/Documents/Github/KA_analyses']));
    inter_dir = ['/Users/' usr_name '/Williams Lab Dropbox/Eric Carmichael/KA_Data/inter_reward_23'];
    save_dir = ['/Users/' usr_name '/Williams Lab Dropbox/Eric Carmichael/KA_Data/Ln_model_out']; 
elseif ispc


    % load data
    addpath(genpath(['C:\Users\' usr_name '\Documents\GitHub\vandermeerlab\code-matlab\shared']))
    addpath(genpath(['C:\Users\' usr_name '\Documents\GitHub\EC_State']));
    addpath(genpath(['C:\Users\' usr_name '\Documents\GitHub\KA_analyses']));
    inter_dir = ['C:\Users\' usr_name '\Williams Lab Dropbox\Eric Carmichael\KA_Data\inter_reward_23'];
    save_dir = ['C:\Users\' usr_name '\Williams Lab Dropbox\Eric Carmichael\KA_Data\Ln_model_out']; 
end


omit_list = {'C4_3_C3_2021-02-25_DONE',...
    'C5_2_O7_2021-04-30_DONE',... % Feeders start way before the recording.
    'C6_3_O1_2021-09-24_DONE',...
    'C6_3_O4_2021-09-27_DONE',...
    'C6_3_O5_2021-09-29_DONE',...
    'C6_4_O6_2021-09-27_DONE',... % HS likly fell out.
    'C6_4_E1_2021-10-07_DONE',... % 3 recordings.
    'C3_2_O7_2020-09-03_DONE',... % spike is all noise.
    };

omit_cells =  {'C1_1 O6_2020-07-12_DONE_maze_data_TT1_01.t64',...
    'C3_3_O2_2020-09-03_DONE_maze_data_TT4_02.t64',...
    'C3_3_O4_2020-09-05_DONE_maze_data_TT4_01.t64',...
    'C3_3_O6_2020-09-07_DONE_maze_data_TT2_01.t64',...
    'C5_3_C3_2021-04-24_DONE_maze_data_TT2_01.t64',...
    'C5_3_O1_2021-04-25_DONE_maze_data_TT3_01.t64',...
    'C5_3_O2_2021-04-26_DONE_maze_data_TT3_01.t64',...
    'C5_3_O4_2021-04-28_DONE_maze_data_TT2_01.t64',...
    'C5_3_O4_2021-04-28_DONE_maze_data_TT4_01.t64',...
    'C6_3_O6_2021-09-30_DONE_maze_data_TT2_02.t64'};

%%
cd(inter_dir)
sess_list = dir([inter_dir filesep '*.mat']);
spd_metrics = []; 
t_metrics = []; 

k = 0;
for iS = 1:length(sess_list)
warning off

     data = load([inter_dir filesep sess_list(iS).name]);

     data = data.data; 

     data.pos = KA_refine_pos([], data.pos) ;

     data.pos = KA_grid_pos([] ,data.pos);
    % get the mean velocity when the animal is moving.
    % mVelo(iS) = mean(data.velo_smooth.data(data.velo_smooth.data>5));

    for iC = 1:length(data.S.t)

        if ismember([sess_list(iS).name(1:end-4) '_' data.S.label{iC}], omit_cells)
            continue
        end
        k = k+1;

        cell_id{k} = [sess_list(iS).name(1:end-4) '_' data.S.label{iC}];

        % run the ln models
        % KA_run_ln(data, iC, sess_list(iS).name(1:end-4), save_dir)

            % KA_get_spatial_info(data, iC)

    dt = mode(diff(data.pos.tvec)); 
    this_S = KA_isolate_S(data.S, data.S.label{iC});
    [spd_metrics{k}, t_metrics{k}, p_metrics{k}] = KA_get_spatial_info(data, iC);

    % if spd_metrics{k}.zMI < 1.




    end

end



%% get the disrtibution of tuning values

spd.MI = []; spd.Is = []; spd.Ispk = []; 
spd.zMI = []; spd.zIs = []; spd.zIspk = []; 

t.MI = []; t.Is = []; t.Ispk = []; 
t.zMI = []; t.zIs = []; t.zIspk = []; 

for ii = length(spd_metrics):-1:1

% raw
spd.MI(ii) = spd_metrics{ii}.MI; 
spd.Is(ii) = spd_metrics{ii}.Isec; 
spd.Ispk(ii) = spd_metrics{ii}.Ispike; 

%norm
spd.zMI(ii) = spd_metrics{ii}.NormMI; 
spd.zIs(ii) = spd_metrics{ii}.NormIsec; 
spd.zIspk(ii) = spd_metrics{ii}.NormIspike; 

% raw
t.MI(ii) = t_metrics{ii}.MI; 
t.Is(ii) = t_metrics{ii}.Isec; 
t.Ispk(ii) = t_metrics{ii}.Ispike; 

%norm
t.zMI(ii) = t_metrics{ii}.NormMI; 
t.zIs(ii) = t_metrics{ii}.NormIsec; 
t.zIspk(ii) = t_metrics{ii}.NormIspike; 



end


figure(1010)
clf
subplot(3,3,1)
histogram(spd.MI,50); 
xlabel('MI')
title('Speed')

subplot(3,3,4)
histogram(spd.zMI, 50); 
xlabel('zMI')

subplot(3,3,7)
histogram(spd.zIs, 50); 
xlabel('zIs')

subplot(3,3,2)
histogram(t.MI, 50); 
xlabel('MI')
title('t-minus')

subplot(3,3,5)
histogram(t.zMI, 50); 
xlabel('zMI')

subplot(3,3,8)
histogram(t.zIs, 50); 
xlabel('zIs')


subplot(3,3, [3 6])
spd_sig = spd.zMI > 2.5; 
t_sig = t.zMI > 2.5; 
st_sig = spd_sig & t_sig; 
cla
hold on
scatter(spd.MI(~t_sig | ~spd_sig), t.MI(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.MI(spd_sig), t.MI(spd_sig), 50, "yellow", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.MI(t_sig), t.MI(t_sig), 50, "blue", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.MI(t_sig & spd_sig), t.MI(t_sig & spd_sig), 50, "green", 'filled') %[spd_sig.*(t_sig+1)])

ylabel('t MI')
xlabel('spd MI')
legend()

subplot(3,3, 9)
spd_sig = spd.zMI > 2.5; 
t_sig = t.zMI > 2.5; 
st_sig = spd_sig & t_sig; 
cla
hold on
scatter(spd.zMI(~t_sig | ~spd_sig), t.zMI(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zMI(spd_sig), t.zMI(spd_sig), 50, "yellow", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zMI(t_sig), t.zMI(t_sig), 50, "blue", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zMI(t_sig & spd_sig), t.zMI(t_sig & spd_sig), 50, "green", 'filled') %[spd_sig.*(t_sig+1)])

ylabel('t zMI')
xlabel('spd zMI')


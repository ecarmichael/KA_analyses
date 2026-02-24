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
p_metrics = []; 

t_metrics_trl = []; 
spd_metrics_trl = []; 
p_metrics_trl = []; 
FR = []; 
phase =[] ; 
sub = []; 

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

    try
    [spd_metrics_trl{k}, t_metrics_trl{k}, p_metrics_trl{k}] = KA_get_spatial_info_trl(data, iC);
    catch
        spd_metrics_trl{k}.M = NaN; spd_metrics_trl{k}.Isec = NaN; spd_metrics_trl{k}.Ispike = NaN;
        t_metrics_trl{k}.M = NaN; t_metrics_trl{k}.Isec = NaN; t_metrics_trl{k}.Ispike = NaN;
        p_metrics_trl{k}.M = NaN; p_metrics_trl{k}.Isec = NaN; p_metrics_trl{k}.Ispike = NaN;

    end
    FR{k} = length(this_S.t{1})./ (data.pos.tvec(end) - data.pos.tvec(1));

    phase{k} = sess_list(iS).name(6:7);
    sub{k} = str2double(sess_list(iS).name([2 4]));


    end

end



%% get the disrtibution of tuning values

spd.MI = []; spd.Is = []; spd.Ispk = []; 
spd.zMI = []; spd.zIs = []; spd.zIspk = []; 

t.MI = []; t.Is = []; t.Ispk = []; 
t.zMI = []; t.zIs = []; t.zIspk = []; 

% p = random shuffles to see what the false discovery rate was. 
p.MI = []; p.Is = []; p.Ispk = []; 
p.zMI = []; p.zIs = []; p.zIspk = []; 

fr = cell2mat(FR); 

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

% raw
p.MI(ii) = p_metrics{ii}.MI; 
p.Is(ii) = p_metrics{ii}.Isec; 
p.Ispk(ii) = p_metrics{ii}.Ispike; 

%norm
p.zMI(ii) = p_metrics{ii}.NormMI; 
p.zIs(ii) = p_metrics{ii}.NormIsec; 
p.zIspk(ii) = p_metrics{ii}.NormIspike; 


end


c_ord = MS_linspecer(5);
bin_c = 50; 
m = 6; n= 5; 

figure(1010)
clf
subplot(m, n ,1)
histogram(spd.MI,bin_c, FaceColor= c_ord(1,:), BinLimits=[0 .5]); 
xlabel('MI')
title('Speed')
xlim([0 .5])
ylabel('count')

subplot(m, n ,6)
histogram(spd.zMI, bin_c, FaceColor= c_ord(1,:), BinLimits=[-5 60]);
xlabel('zMI')
xlim([-5 60])
ylabel('count')
xline([-3.29 3.29])

subplot(m, n,11)
histogram(spd.Is, bin_c, FaceColor= c_ord(1,:), BinLimits=[0 0.015]);
xlabel('Is')
ylabel('count')
xlim([0 0.015])

subplot(m, n,16)
histogram(spd.zIs, bin_c, FaceColor= c_ord(1,:), BinLimits=[-5 60]);
xlabel('zIs')
xlim([-5 60])
ylabel('count')
xline([-3.29 3.29])

subplot(m, n,21)
histogram(spd.Ispk, bin_c, FaceColor= c_ord(1,:), BinLimits=[0 1.25]); 
xlabel('Ispk')
xlim([0 1.25])
ylabel('count')

subplot(m, n,26)
histogram(spd.zIspk, bin_c, FaceColor= c_ord(1,:), BinLimits=[-5 60]); 
xlabel('zIspk')
xlim([-5 60])
ylabel('count')
xline([-3.29 3.29])

subplot(m, n,2)
histogram(t.MI, bin_c, FaceColor= c_ord(2,:), BinLimits=[0 .5]); 
xlabel('MI')
title('t-minus')
xlim([0 .5])

subplot(m, n,7)
histogram(t.zMI, bin_c, FaceColor= c_ord(2,:), BinLimits=[-5 60]);
xlabel('zMI')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,12)
histogram(t.Is, bin_c, FaceColor= c_ord(2,:), BinLimits=[0 0.015]); 
xlabel('zIs')
xlim([0 0.015])

subplot(m, n,17)
histogram(t.zIs, bin_c, FaceColor= c_ord(2,:), BinLimits=[-5 60]); 
xlabel('zIs')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,22)
histogram(t.Ispk, bin_c, FaceColor= c_ord(2,:), BinLimits=[0 1.25]); 
xlabel('Ispk')
xlim([0 1.25])

subplot(m, n,27)
histogram(t.zIspk, bin_c, FaceColor= c_ord(2,:), BinLimits=[-5 60]);
xlabel('zIspk')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,3)
histogram(p.MI, bin_c,FaceColor='k', BinLimits=[0 .5]); 
xlabel('MI')
title('Shuff')
xlim([0 .5])

subplot(m, n,8)
histogram(p.zMI, bin_c,FaceColor='k', BinLimits=[-5 60]);
xlabel('zMI')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,13)
histogram(p.Is, bin_c,FaceColor='k', BinLimits=[0 0.015]);
xlabel('Is')
xlim([0 0.015])

subplot(m, n,18)
histogram(p.zIs, bin_c,FaceColor='k', BinLimits=[-5 60]);
xlabel('zIs')
xlim([-5 60])
xline([-3.29 3.29])


subplot(m, n,23)
histogram(p.Ispk, bin_c,FaceColor='k', BinLimits=[0 1.25]); 
xlabel('Is')
xlim([0 1.25])

subplot(m, n,28)
histogram(p.zIspk, bin_c,FaceColor='k', BinLimits=[-5 60]);
xlabel('zIspk')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n, [4 9])
spd_sig = spd.zMI > 3.29; 
t_sig = t.zMI > 3.29; 
st_sig = spd_sig & t_sig; 

cla
hold on
scatter(spd.zMI(~t_sig | ~spd_sig), t.zMI(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zMI(spd_sig), t.zMI(spd_sig), 50, c_ord(1,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zMI(t_sig), t.zMI(t_sig), 50, c_ord(2,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zMI(t_sig & spd_sig), t.zMI(t_sig & spd_sig), 50, c_ord(5,:), 'filled') %[spd_sig.*(t_sig+1)])

ylabel('t zMI')
xlabel('spd MI')
xlim([min(spd.zMI) max(spd.zMI)])
ylim([min(t.zMI) max(t.zMI)])
set(gca, 'XScale', 'log', 'YScale', 'log')
xline(3.29,'HandleVisibility','off'); yline(3.29,'HandleVisibility','off')
legend({"none", 'Spd mod', 'T2R mod', 'Spd + T2R mod'}, 'box', 'off', 'location', 'northwest')
% subtitle('MI')

subplot(m, n, [5 10])
cla
hold on
scatter3(spd.zMI(~t_sig | ~spd_sig), t.zMI(~t_sig | ~spd_sig),fr(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter3(spd.zMI(spd_sig), t.zMI(spd_sig),fr(spd_sig), 50, c_ord(1,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter3(spd.zMI(t_sig), t.zMI(t_sig), fr(t_sig), 50, c_ord(2,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter3(spd.zMI(t_sig & spd_sig), t.zMI(t_sig & spd_sig),fr(t_sig & spd_sig), 50, c_ord(5,:), 'filled') %[spd_sig.*(t_sig+1)])

zlabel('firing rate')
ylabel('t zMI')
xlabel('spd MI')
xlim([min(spd.zMI) max(spd.zMI)])
ylim([min(t.zMI) max(t.zMI)])
set(gca, 'XScale', 'log', 'YScale', 'log')
% xline(3.29,'HandleVisibility','off'); yline(3.29,'HandleVisibility','off')
% legend({"none", 'Spd mod', 'T2R mod', 'Spd + T2R mod'}, 'box', 'off', 'location', 'northwest')
% subtitle('MI')
view(gca,[-45 -45, 45])

subplot(m, n, [14 19])
spd_sig = spd.zIs > 3.29; 
t_sig = t.zIs > 3.29; 
st_sig = spd_sig & t_sig; 
cla
hold on
scatter(spd.zIs(~t_sig | ~spd_sig), t.zIs(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIs(spd_sig), t.zIs(spd_sig), 50, c_ord(1,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIs(t_sig), t.zIs(t_sig), 50, c_ord(2,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIs(t_sig & spd_sig), t.zIs(t_sig & spd_sig), 50,c_ord(5,:), 'filled') %[spd_sig.*(t_sig+1)])

ylabel('t zIs')
xlabel('spd zIs')
set(gca, 'XScale', 'log', 'YScale', 'log')
xline(3.29,'HandleVisibility','off'); yline(3.29,'HandleVisibility','off')
% legend({"none", 'Spd mod', 'T2R mod', 'Spd + T2R mod'}, 'box', 'off', 'location', 'northwest')
ylim([10^-2 10^2])
xlim([10^-2 10^2])

% FR
subplot(m, n, [15 20])
cla
hold on
scatter3(spd.zIs(~t_sig | ~spd_sig), t.zIs(~t_sig | ~spd_sig),fr(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter3(spd.zIs(spd_sig), t.zIs(spd_sig),fr(spd_sig), 50, c_ord(1,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter3(spd.zIs(t_sig), t.zIs(t_sig), fr(t_sig), 50, c_ord(2,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter3(spd.zIs(t_sig & spd_sig), t.zIs(t_sig & spd_sig),fr(t_sig & spd_sig), 50, c_ord(5,:), 'filled') %[spd_sig.*(t_sig+1)])

zlabel('firing rate')
ylabel('t zIs')
xlabel('spd zIs')
% xlim([min(spd.zIs) max(spd.zIs)])
ylim([10^-2 10^2])
xlim([10^-2 10^2])

set(gca, 'XScale', 'log', 'YScale', 'log')
view(gca,[-45 -45, 45])

subplot(m, n, [24 29])
spd_sig = spd.zIspk > 3.29; 
t_sig = t.zIspk > 3.29; 
st_sig = spd_sig & t_sig; 
cla
hold on
scatter(spd.zIspk(~t_sig | ~spd_sig), t.zIspk(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIspk(spd_sig), t.zIspk(spd_sig), 50, c_ord(1,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIspk(t_sig), t.zIspk(t_sig), 50, c_ord(2,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIspk(t_sig & spd_sig), t.zIspk(t_sig & spd_sig), 50,c_ord(5,:), 'filled') %[spd_sig.*(t_sig+1)])

ylabel('t zIspk (log)')
xlabel('spd zIspk (log)')
xline(3.29,'HandleVisibility','off'); yline(3.29,'HandleVisibility','off')
% legend({"none", 'Spd mod', 'T2R mod', 'Spd + T2R mod'}, 'box', 'off', 'location', 'northwest')
ylim([10^-1 10^2])
xlim([10^-1 10^2])
set(gca, 'XScale', 'log', 'YScale', 'log')

subplot(m, n, [25 30])
cla
hold on
scatter3(spd.zIspk(~t_sig | ~spd_sig), t.zIspk(~t_sig | ~spd_sig),fr(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter3(spd.zIspk(spd_sig), t.zIspk(spd_sig),fr(spd_sig), 50, c_ord(1,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter3(spd.zIspk(t_sig), t.zIspk(t_sig), fr(t_sig), 50, c_ord(2,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter3(spd.zIspk(t_sig & spd_sig), t.zIspk(t_sig & spd_sig),fr(t_sig & spd_sig), 50, c_ord(5,:), 'filled') %[spd_sig.*(t_sig+1)])

zlabel('firing rate')
ylabel('t zIspk')
xlabel('spd zIspk')
ylim([10^-1 10^2])
xlim([10^-1 10^2])
set(gca, 'XScale', 'log', 'YScale', 'log')
view(gca,[-45 -45, 45])

%% trial type


spd_trl.MI = []; spd_trl.Is = []; spd_trl.Ispk = []; 
spd_trl.zMI = []; spd_trl.zIs = []; spd_trl.zIspk = []; 

t_trl.MI = []; t_trl.Is = []; t_trl.Ispk = []; 
t_trl.zMI = []; t_trl.zIs = []; t_trl.zIspk = []; 

% p = random shuffles to see what the false discovery rate was. 
p_trl.MI = []; p_trl.Is = []; p_trl.Ispk = []; 
p_trl.zMI = []; p_trl.zIs = []; p_trl.zIspk = []; 

for ii = length(spd_metrics):-1:1
 
    for tt = 1:4
% raw
spd_trl.MI(ii, tt) = spd_metrics_trl{ii}{tt}.MI; 
spd_trl.Is(ii, tt) = spd_metrics_trl{ii}{tt}.Isec; 
spd_trl.Ispk(ii, tt) = spd_metrics_trl{ii}{tt}.Ispike; 

%norm
spd_trl.zMI(ii, tt) = spd_metrics_trl{ii}{tt}.NormMI; 
spd_trl.zIs(ii, tt) = spd_metrics_trl{ii}{tt}.NormIsec; 
spd_trl.zIspk(ii, tt) = spd_metrics_trl{ii}{tt}.NormIspike; 


% raw
t_trl.MI(ii, tt) = t_metrics_trl{ii}{tt}.MI; 
t_trl.Is(ii, tt) = t_metrics_trl{ii}{tt}.Isec; 
t_trl.Ispk(ii, tt) = t_metrics_trl{ii}{tt}.Ispike; 

%norm
t_trl.zMI(ii, tt) = t_metrics_trl{ii}{tt}.NormMI; 
t_trl.zIs(ii, tt) = t_metrics_trl{ii}{tt}.NormIsec; 
t_trl.zIspk(ii, tt) = t_metrics_trl{ii}{tt}.NormIspike; 

% raw
p_trl.MI(ii, tt) = p_metrics_trl{ii}{tt}.MI; 
p_trl.Is(ii, tt) = p_metrics_trl{ii}{tt}.Isec; 
p_trl.Ispk(ii, tt) = p_metrics_trl{ii}{tt}.Ispike; 

%norm
p_trl.zMI(ii, tt) = p_metrics_trl{ii}{tt}.NormMI; 
p_trl.zIs(ii, tt) = p_metrics_trl{ii}{tt}.NormIsec; 
p_trl.zIspk(ii, tt) = p_metrics_trl{ii}{tt}.NormIspike;

    end
end

figure(1001)
clf
jitter = (rand(1, length(spd_trl.MI))*.5); 
x_vec = [ones(1, length(spd_trl.MI)).*jitter+.75, ones(1, length(spd_trl.MI)).*jitter+1.75, ones(1, length(spd_trl.MI)).*jitter+2.75, ones(1, length(spd_trl.MI)).*jitter+3.75]; 
m = 6; n= 3; 
sz = 10; 

subplot(m, n ,1)
scatter(x_vec, [spd_trl.MI(:,1); spd_trl.MI(:,2); spd_trl.MI(:,3); spd_trl.MI(:,4)],sz, c_ord(1,:), 'filled')
set(gca, 'XTick', [1:4]);
ylim([0 3]); 
ylabel('MI')
set(gca, 'XTick', []);
title('Speed modulation')

subplot(m, n ,2)
scatter(x_vec, [t_trl.MI(:,1); t_trl.MI(:,2); t_trl.MI(:,3); t_trl.MI(:,4)],sz, c_ord(2,:), 'filled')
set(gca, 'XTick', [1:4]);
ylim([0 3]); 
set(gca, 'XTick', [], 'ytick', []);
title('time to reward modulation')

subplot(m, n ,3)
scatter(x_vec, [p_trl.MI(:,1); p_trl.MI(:,2); p_trl.MI(:,3); p_trl.MI(:,4)],sz, 'k', 'filled')
set(gca, 'XTick', [1:4]);
ylim([0 3]); 
set(gca, 'XTick', [], 'ytick', []);
title('Shuffles')

subplot(m, n ,4)
scatter(x_vec, [spd_trl.zMI(:,1); spd_trl.zMI(:,2); spd_trl.zMI(:,3); spd_trl.zMI(:,4)],sz, c_ord(1,:), 'filled')
set(gca, 'XTick', [1:4]);
ylim([-3.29 30]); 
ylabel('zscore MI')
set(gca, 'XTick', []);

subplot(m, n ,5)
scatter(x_vec, [t_trl.zMI(:,1); t_trl.zMI(:,2); t_trl.zMI(:,3); t_trl.zMI(:,4)],sz, c_ord(2,:), 'filled')
set(gca, 'XTick', [1:4]);
ylim([-3.29 30]); 
set(gca, 'XTick', [], 'ytick', []);

subplot(m, n ,6)
scatter(x_vec, [p_trl.zMI(:,1); p_trl.zMI(:,2); p_trl.zMI(:,3); p_trl.zMI(:,4)],sz, 'k', 'filled')
set(gca, 'XTick', 1:4);
ylim([-3.29 30]); 
set(gca, 'XTick', [], 'ytick', []);

% iS
subplot(m, n ,7)
scatter(x_vec, [spd_trl.Is(:,1); spd_trl.Is(:,2); spd_trl.Is(:,3); spd_trl.Is(:,4)],sz, c_ord(1,:), 'filled')
set(gca, 'XTick', 1:4);
ylim([0 0.1]); 
ylabel('info/sec')
set(gca, 'XTick', []);

subplot(m, n ,8)
scatter(x_vec, [t_trl.Is(:,1); t_trl.Is(:,2); t_trl.Is(:,3); t_trl.Is(:,4)],sz, c_ord(2,:), 'filled')
ylim([0 0.1]); 
set(gca, 'XTick', [], 'ytick', []);

subplot(m, n ,9)
scatter(x_vec, [p_trl.Is(:,1); p_trl.Is(:,2); p_trl.Is(:,3); p_trl.Is(:,4)],sz, 'k', 'filled')
ylim([0 0.1]); 
set(gca, 'XTick', [], 'ytick', []);

subplot(m, n ,10)
scatter(x_vec, [spd_trl.zIs(:,1); spd_trl.zIs(:,2); spd_trl.zIs(:,3); spd_trl.zIs(:,4)],sz, c_ord(1,:), 'filled')
ylim([-3.29 30]); 
ylabel('zscore info/sec')
set(gca, 'XTick', []);

subplot(m, n ,11)
scatter(x_vec, [t_trl.zIs(:,1); t_trl.zIs(:,2); t_trl.zIs(:,3); t_trl.zIs(:,4)],sz, c_ord(2,:), 'filled')
ylim([-3.29 30]); 
set(gca, 'XTick', [], 'ytick', []);

subplot(m, n ,12)
scatter(x_vec, [p_trl.zIs(:,1); p_trl.zIs(:,2); p_trl.zIs(:,3); p_trl.zIs(:,4)],sz, 'k', 'filled')
ylim([-3.29 30]); 
set(gca, 'XTick', [], 'ytick', []);

% iSpike
subplot(m, n ,13)
scatter(x_vec, [spd_trl.Ispk(:,1); spd_trl.Ispk(:,2); spd_trl.Ispk(:,3); spd_trl.Ispk(:,4)],sz, c_ord(1,:), 'filled')
ylim([0 6]); 
ylabel('zscore info/spike')
set(gca, 'XTick', []);

subplot(m, n ,14)
scatter(x_vec, [t_trl.Ispk(:,1); t_trl.Ispk(:,2); t_trl.Ispk(:,3); t_trl.Ispk(:,4)],sz, c_ord(2,:), 'filled')
ylim([0 6]); 
set(gca, 'XTick', [], 'ytick', []);

subplot(m, n ,15)
scatter(x_vec, [p_trl.Ispk(:,1); p_trl.Ispk(:,2); p_trl.Ispk(:,3); p_trl.Ispk(:,4)],sz, 'k', 'filled')
ylim([0 6]); 
set(gca, 'XTick', [], 'ytick', []);

subplot(m, n ,16)
scatter(x_vec, [spd_trl.zIspk(:,1); spd_trl.zIspk(:,2); spd_trl.zIspk(:,3); spd_trl.zIspk(:,4)],sz, c_ord(1,:), 'filled')
ylim([-3.29 30]); 
ylabel('zscore info/spike')
set(gca, 'XTick', 1:4);

subplot(m, n ,17)
scatter(x_vec, [t_trl.zIspk(:,1); t_trl.zIspk(:,2); t_trl.zIspk(:,3); t_trl.zIspk(:,4)],sz, c_ord(2,:), 'filled')
ylim([-3.29 30]); 
set(gca, 'XTick', 1:4, 'ytick', []);

subplot(m, n ,18)
scatter(x_vec, [p_trl.zIspk(:,1); p_trl.zIspk(:,2); p_trl.zIspk(:,3); p_trl.zIspk(:,4)],sz, 'k', 'filled')
ylim([-3.29 30]); 
set(gca, 'XTick', 1:4, 'ytick', []);


%% convert it to a table;
%simple
info_tbl = table(sub', phase', fr',spd.MI', t.MI', p.MI', spd.zMI', t.zMI', p.zMI',...
spd.Is', t.Is', p.Is', spd.zIs', t.zIs', p.zIs',...
spd.Ispk', t.Ispk', p.Ispk', spd.zIspk', t.zIspk', p.zIspk',...
    'VariableNames', {'Sub ID','Phase', 'Firing_rate', 'spd_MI', 't_MI', 'shuff_MI', 'spd_zMI', 't_zMI', 'shuff_zMI',...
    'spd_Is', 't_Is', 'shuff_Is', 'spd_zIs', 't_zIs', 'shuff_zIs',...
    'spd_Ispk', 't_Ispk', 'shuff_Ispk', 'spd_zIspk', 't_zIspk', 'shuff_zIspk'}); 

writetable(info_tbl,['info_metrics_tbl.csv'] )


% reshape things
c_id = repmat(1:length(spd_trl.MI),1, 4); 

f_id = [ones(1,length(spd_trl.MI)), ones(1,length(spd_trl.MI))*2, ones(1,length(spd_trl.MI))*3, ones(1,length(spd_trl.MI))*4];
f_names = cell(size(f_id)); 
f_mag = f_names; 
for ii =1 :4 
    f_names(f_id==ii) = data.PM.Feeder_type(ii); 
    f_mag(f_id==ii) = {num2str(data.PM.Feeder_mag(ii))}; 
end

spd_MI = reshape(spd_trl.MI,1, numel(spd_trl.MI)); 
spd_zMI = reshape(spd_trl.zMI,1, numel(spd_trl.zMI)); 
spd_Is = reshape(spd_trl.Is,1, numel(spd_trl.Is)); 
spd_zIs = reshape(spd_trl.zIs,1, numel(spd_trl.zIs)); 
spd_Ispk = reshape(spd_trl.Ispk,1, numel(spd_trl.Ispk)); 
spd_zIspk = reshape(spd_trl.zIspk,1, numel(spd_trl.zIspk)); 

t_MI = reshape(t_trl.MI,1, numel(t_trl.MI)); 
t_zMI = reshape(t_trl.zMI,1, numel(t_trl.zMI)); 
t_Is = reshape(t_trl.Is,1, numel(t_trl.Is)); 
t_zIs = reshape(t_trl.zIs,1, numel(t_trl.zIs)); 
t_Ispk = reshape(t_trl.Ispk,1, numel(t_trl.Ispk)); 
t_zIspk = reshape(t_trl.zIspk,1, numel(t_trl.zIspk)); 

p_MI = reshape(p_trl.MI,1, numel(p_trl.MI)); 
p_zMI = reshape(p_trl.zMI,1, numel(p_trl.zMI)); 
p_Is = reshape(p_trl.Is,1, numel(p_trl.Is)); 
p_zIs = reshape(p_trl.zIs,1, numel(p_trl.zIs)); 
p_Ispk = reshape(p_trl.Ispk,1, numel(p_trl.Ispk)); 
p_zIspk = reshape(p_trl.zIspk,1, numel(p_trl.zIspk)); 


info_trl_tbl = table(repmat(cell2mat(sub),1, 4)', repmat(phase,1, 4)',c_id',  f_id', f_names', f_mag', repmat(fr, 1,4)',spd_MI', t_MI', p_MI', spd_zMI', t_zMI', p_zMI',...
spd_Is', t_Is', p_Is', spd_zIs', t_zIs', p_zIs',...
spd_Ispk', t_Ispk', p_Ispk', spd_zIspk', t_zIspk', p_zIspk',...
    'VariableNames', {'Sub ID','Phase' 'Cell ID', 'Feeder_id', 'Feeder_type', 'Feeder_mag', 'Firing_rate', 'spd_MI', 't_MI', 'shuff_MI', 'spd_zMI', 't_zMI', 'shuff_zMI',...
    'spd_Is', 't_Is', 'shuff_Is', 'spd_zIs', 't_zIs', 'shuff_zIs',...
    'spd_Ispk', 't_Ispk', 'shuff_Ispk', 'spd_zIspk', 't_zIspk', 'shuff_zIspk'}); 

writetable(info_trl_tbl,'info_metrics_trial_tbl.csv')


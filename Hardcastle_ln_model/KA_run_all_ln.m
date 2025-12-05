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

    [spd_metrics_trl{k}, t_metrics_trl{k}, p_metrics_trl{k}] = KA_get_spatial_info_trl(data, iC);

    % if spd_metrics{k}.zMI < 1.




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
m = 6; n= 4; 

figure(1010)
clf
subplot(m, n ,1)
histogram(spd.MI,bin_c, FaceColor= c_ord(1,:), BinLimits=[0 .5]); 
xlabel('MI')
title('Speed')
xlim([0 .5])
ylabel('count')

subplot(m, n ,5)
histogram(spd.zMI, bin_c, FaceColor= c_ord(1,:), BinLimits=[-5 60]);
xlabel('zMI')
xlim([-5 60])
ylabel('count')
xline([-3.29 3.29])

subplot(m, n,9)
histogram(spd.Is, bin_c, FaceColor= c_ord(1,:), BinLimits=[0 0.015]);
xlabel('Is')
ylabel('count')
xlim([0 0.015])

subplot(m, n,13)
histogram(spd.zIs, bin_c, FaceColor= c_ord(1,:), BinLimits=[-5 60]);
xlabel('zIs')
xlim([-5 60])
ylabel('count')
xline([-3.29 3.29])

subplot(m, n,17)
histogram(spd.Ispk, bin_c, FaceColor= c_ord(1,:), BinLimits=[0 1.25]); 
xlabel('Ispk')
xlim([0 1.25])
ylabel('count')

subplot(m, n,21)
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

subplot(m, n,6)
histogram(t.zMI, bin_c, FaceColor= c_ord(2,:), BinLimits=[-5 60]);
xlabel('zMI')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,10)
histogram(t.Is, bin_c, FaceColor= c_ord(2,:), BinLimits=[0 0.015]); 
xlabel('zIs')
xlim([0 0.015])

subplot(m, n,14)
histogram(t.zIs, bin_c, FaceColor= c_ord(2,:), BinLimits=[-5 60]); 
xlabel('zIs')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,18)
histogram(t.Ispk, bin_c, FaceColor= c_ord(2,:), BinLimits=[0 1.25]); 
xlabel('Ispk')
xlim([0 1.25])

subplot(m, n,22)
histogram(t.zIspk, bin_c, FaceColor= c_ord(2,:), BinLimits=[-5 60]);
xlabel('zIspk')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,3)
histogram(p.MI, bin_c,FaceColor='k', BinLimits=[0 .5]); 
xlabel('MI')
title('Shuff')
xlim([0 .5])

subplot(m, n,7)
histogram(p.zMI, bin_c,FaceColor='k', BinLimits=[-5 60]);
xlabel('zMI')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,11)
histogram(p.Is, bin_c,FaceColor='k', BinLimits=[0 0.015]);
xlabel('Is')
xlim([0 0.015])

subplot(m, n,15)
histogram(p.zIs, bin_c,FaceColor='k', BinLimits=[-5 60]);
xlabel('zIs')
xlim([-5 60])
xline([-3.29 3.29])


subplot(m, n,19)
histogram(p.Ispk, bin_c,FaceColor='k', BinLimits=[0 1.25]); 
xlabel('Is')
xlim([0 1.25])

subplot(m, n,23)
histogram(p.zIspk, bin_c,FaceColor='k', BinLimits=[-5 60]);
xlabel('zIspk')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n, [4 8])
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

subplot(m, n, [12 16])
cla
hold on
scatter(spd.zIs(~t_sig | ~spd_sig), t.zIs(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIs(spd_sig), t.zIs(spd_sig), 50, c_ord(1,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIs(t_sig), t.zIs(t_sig), 50, c_ord(2,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIs(t_sig & spd_sig), t.zIs(t_sig & spd_sig), 50,c_ord(5,:), 'filled') %[spd_sig.*(t_sig+1)])

ylabel('t zIs')
xlabel('spd zIs')
xlim([min(spd.zIs) max(spd.zIs)])
ylim([min(t.zIs) max(t.zIs)])
set(gca, 'XScale', 'log', 'YScale', 'log')
xline(3.29,'HandleVisibility','off'); yline(3.29,'HandleVisibility','off')
legend({"none", 'Spd mod', 'T2R mod', 'Spd + T2R mod'}, 'box', 'off', 'location', 'northwest')


subplot(m, n, [20 24])
spd_sig = spd.zIs > 3.29; 
t_sig = t.zIs > 3.29; 
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
legend({"none", 'Spd mod', 'T2R mod', 'Spd + T2R mod'}, 'box', 'off', 'location', 'northwest')
xlim([min(spd.zIspk) max(spd.zIspk)])
ylim([min(t.zIspk) max(t.zIspk)])
set(gca, 'XScale', 'log', 'YScale', 'log')


%% trial type


spd_trL.MI = []; spd_trL.Is = []; spd_trL.Ispk = []; 
spd_trL.zMI = []; spd_trL.zIs = []; spd_trL.zIspk = []; 

t_trl.MI = []; t_trl.Is = []; t_trl.Ispk = []; 
t_trl.zMI = []; t_trl.zIs = []; t_trl.zIspk = []; 

% p = random shuffles to see what the false discovery rate was. 
p_trl.MI = []; p_trl.Is = []; p_trl.Ispk = []; 
p_trl.zMI = []; p_trl.zIs = []; p_trl.zIspk = []; 

for ii = length(spd_metrics):-1:1
 
    for tt = 1:4
% raw
spd_trL.MI(ii, tt) = spd_metrics_trl{ii}{tt}.MI; 
spd_trL.Is(ii, tt) = spd_metrics_trl{ii}{tt}.Isec; 
spd_trL.Ispk(ii, tt) = spd_metrics_trl{ii}{tt}.Ispike; 

%norm
spd_trL.zMI(ii, tt) = spd_metrics_trl{ii}{tt}.NormMI; 
spd_trL.zIs(ii, tt) = spd_metrics_trl{ii}{tt}.NormIsec; 
spd_trL.zIspk(ii, tt) = spd_metrics_trl{ii}{tt}.NormIspike; 


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


c_ord = MS_linspecer(5);
bin_c = 50; 
m = 6; n= 4; 

figure(1010)
clf
subplot(m, n ,1)
histogram(spd.MI,bin_c, FaceColor= c_ord(1,:), BinLimits=[0 .5]); 
xlabel('MI')
title('Speed')
xlim([0 .5])
ylabel('count')

subplot(m, n ,5)
histogram(spd.zMI, bin_c, FaceColor= c_ord(1,:), BinLimits=[-5 60]);
xlabel('zMI')
xlim([-5 60])
ylabel('count')
xline([-3.29 3.29])

subplot(m, n,9)
histogram(spd.Is, bin_c, FaceColor= c_ord(1,:), BinLimits=[0 0.015]);
xlabel('Is')
ylabel('count')
xlim([0 0.015])

subplot(m, n,13)
histogram(spd.zIs, bin_c, FaceColor= c_ord(1,:), BinLimits=[-5 60]);
xlabel('zIs')
xlim([-5 60])
ylabel('count')
xline([-3.29 3.29])

subplot(m, n,17)
histogram(spd.Ispk, bin_c, FaceColor= c_ord(1,:), BinLimits=[0 1.25]); 
xlabel('Ispk')
xlim([0 1.25])
ylabel('count')

subplot(m, n,21)
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

subplot(m, n,6)
histogram(t.zMI, bin_c, FaceColor= c_ord(2,:), BinLimits=[-5 60]);
xlabel('zMI')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,10)
histogram(t.Is, bin_c, FaceColor= c_ord(2,:), BinLimits=[0 0.015]); 
xlabel('zIs')
xlim([0 0.015])

subplot(m, n,14)
histogram(t.zIs, bin_c, FaceColor= c_ord(2,:), BinLimits=[-5 60]); 
xlabel('zIs')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,18)
histogram(t.Ispk, bin_c, FaceColor= c_ord(2,:), BinLimits=[0 1.25]); 
xlabel('Ispk')
xlim([0 1.25])

subplot(m, n,22)
histogram(t.zIspk, bin_c, FaceColor= c_ord(2,:), BinLimits=[-5 60]);
xlabel('zIspk')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,3)
histogram(p.MI, bin_c,FaceColor='k', BinLimits=[0 .5]); 
xlabel('MI')
title('Shuff')
xlim([0 .5])

subplot(m, n,7)
histogram(p.zMI, bin_c,FaceColor='k', BinLimits=[-5 60]);
xlabel('zMI')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n,11)
histogram(p.Is, bin_c,FaceColor='k', BinLimits=[0 0.015]);
xlabel('Is')
xlim([0 0.015])

subplot(m, n,15)
histogram(p.zIs, bin_c,FaceColor='k', BinLimits=[-5 60]);
xlabel('zIs')
xlim([-5 60])
xline([-3.29 3.29])


subplot(m, n,19)
histogram(p.Ispk, bin_c,FaceColor='k', BinLimits=[0 1.25]); 
xlabel('Is')
xlim([0 1.25])

subplot(m, n,23)
histogram(p.zIspk, bin_c,FaceColor='k', BinLimits=[-5 60]);
xlabel('zIspk')
xlim([-5 60])
xline([-3.29 3.29])

subplot(m, n, [4 8])
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

subplot(m, n, [12 16])
cla
hold on
scatter(spd.zIs(~t_sig | ~spd_sig), t.zIs(~t_sig | ~spd_sig), 50, "black", 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIs(spd_sig), t.zIs(spd_sig), 50, c_ord(1,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIs(t_sig), t.zIs(t_sig), 50, c_ord(2,:), 'filled') %[spd_sig.*(t_sig+1)])
scatter(spd.zIs(t_sig & spd_sig), t.zIs(t_sig & spd_sig), 50,c_ord(5,:), 'filled') %[spd_sig.*(t_sig+1)])

ylabel('t zIs')
xlabel('spd zIs')
xlim([min(spd.zIs) max(spd.zIs)])
ylim([min(t.zIs) max(t.zIs)])
set(gca, 'XScale', 'log', 'YScale', 'log')
xline(3.29,'HandleVisibility','off'); yline(3.29,'HandleVisibility','off')
legend({"none", 'Spd mod', 'T2R mod', 'Spd + T2R mod'}, 'box', 'off', 'location', 'northwest')


subplot(m, n, [20 24])
spd_sig = spd.zIs > 3.29; 
t_sig = t.zIs > 3.29; 
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
legend({"none", 'Spd mod', 'T2R mod', 'Spd + T2R mod'}, 'box', 'off', 'location', 'northwest')
xlim([min(spd.zIspk) max(spd.zIspk)])
ylim([min(t.zIspk) max(t.zIspk)])
set(gca, 'XScale', 'log', 'YScale', 'log')

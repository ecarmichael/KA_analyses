%% try some plots of the events data


%% convert jittery signals to first ts

enames = {'t_state1', 't_tone1', 't_valve1','t_state2','t_tone2', 't_valve2', 't_laser_on'}; 
for iE = 1:length(enames)
    
    this_diff = diff(events.(enames{iE})); 
    jump_idx = find(this_diff > 1); 
    
    events_clean.(enames{iE}) = events.(enames{iE})([1 jump_idx+1]);

end

%% generate some rates
bins = 0.05;
t = 0:bins:100;

tbin_edges = t:bins:t(end);
tbin_centers = tbin_edges(1:end-1)+bins/2;

gauss_window = 1./bins; % 1 second window
gauss_SD = 0.02./bins; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./bins; % normalize by binsize


% licks 1
lk_count = histc(events.licks_1,tbin_edges);
lk_count = lk_count(1:end-1);

S_gau_sdf = conv2(lk_count,gk,'same'); % convolve with gaussian window

lk1 = [S_gau_sdf S_gau_sdf(end)];
        
% licks 2
lk_count = histc(events.licks_2,tbin_edges);
lk_count = lk_count(1:end-1);

S_gau_sdf = conv2(lk_count,gk,'same'); % convolve with gaussian window

lk2 = [S_gau_sdf S_gau_sdf(end)];


% interp to a finer resolution

bins = 0.0001;
ti = t(1):bins:t(end);

lk1 = interp1(t, lk1, ti);
lk2 = interp1(t, lk2, ti);

%%
c_ord = linspecer(9);
% c_ord = c_ord([1:3, 5:end],:)

figure(101)
clf
ax(2) = subplot(4,1,3);
hold on
plot([events.licks_1; events.licks_1], [ones(length(events.licks_1),1)-.5 ones(length(events.licks_1),1)+.5]', 'color', c_ord(1,:));
plot([events.licks_2; events.licks_2], [ones(length(events.licks_2),1)+.5 ones(length(events.licks_2),1)+1.5]',  'color', c_ord(2,:));

xlim([t(1) t(end)]);
ylim([.5 2.5])
set(gca, 'ytick', 1:2, 'Yticklabels', {'licks 1', 'lick 2'})

ax(1) = subplot(4,1,1:2);
hold on

enames = {'t_state1', 't_tone1', 't_valve1','t_state2','t_tone2', 't_valve2', 't_laser_on'}; 
for iE = 1:length(enames)
plot([events.(enames{iE}); events.(enames{iE})], [ones(length(events.(enames{iE})),1)+(iE-1)-.5 ones(length(events.(enames{iE})),1)+(iE-1)+.5]', 'color', c_ord(iE+2,:), 'LineWidth', 3);
end
set(gca, 'ytick', 1:length(enames), 'Yticklabels', strrep(strrep(enames, 't_', ''), '_', ' '))
% legend(strrep(enames, 't_', ''), 'Orientation', 'horizontal')
xlim([t(1) t(end)])
ylim([.5 length(enames)+.5])


ax(3) = subplot(4,1,4);
hold on
plot(ti, lk1,  'color', c_ord(1,:))

plot(ti, lk2,  'color', c_ord(2,:))
h = bar(tbin_centers,lk_count);
set(h, 'BarWidth', 1)

ylabel('licks/s')
xlabel('time (s)')

linkaxes(ax, 'x')


%% state based per-event

% get the average lick rats for states 1 and 2

win = [-5 5];

% valve 1 
lk1_v1 = []; lk2_v1 = [];
for ii = length(events_clean.t_valve1):-1:1
   lk1_v1 = [lk1_v1; lk1(nearest_idx3(events_clean.t_valve1(ii)+win(1), ti):nearest_idx3(events_clean.t_valve1(ii)+win(2), ti))];
   lk2_v1 = [lk2_v1;lk2(nearest_idx3(events_clean.t_valve1(ii)+win(1), ti):nearest_idx3(events_clean.t_valve1(ii)+win(2), ti))];
end
    

% valve 1 
lk1_v2 = []; lk2_v2 = [];
for ii = 1:length(events_clean.t_valve2)
   lk1_v2 = [lk1_v2; lk1(nearest_idx3(events_clean.t_valve2(ii)+win(1), ti):nearest_idx3(events_clean.t_valve2(ii)+win(2), ti))];
   lk2_v2 = [lk2_v2;lk2(nearest_idx3(events_clean.t_valve2(ii)+win(1), ti):nearest_idx3(events_clean.t_valve2(ii)+win(2), ti))];
end


figure(601)
subplot(2,1,1)
title('Valve 1')
hold on
plot(-5:bins:5, mean(lk1_v1), 'color', c_ord(1,:), 'LineWidth', 2.5);
plot(-5:bins:5, mean(lk2_v1), 'color', c_ord(2,:), 'LineWidth', 2.5);
legend({'lick 1', 'lick 2'})

% plot(-5:bins:5, mean(lk1_v1) + std(lk1_v1)/sqrt(length(lk1_v1)), '--', 'color', [c_ord(1,:) .3]);
% plot(-5:bins:5, mean(lk1_v1) - std(lk1_v1)/sqrt(length(lk1_v1)), '--', 'color', [c_ord(1,:) .3]);

% plot(-5:bins:5, mean(lk2_v1) + std(lk2_v1)/sqrt(length(lk2_v1)), '--', 'color', [c_ord(2,:) .3]);
% plot(-5:bins:5, mean(lk2_v1) - std(lk2_v1)/sqrt(length(lk2_v1)), '--', 'color', [c_ord(2,:), .3]);
vline(0, 'k')

subplot(2,1,2)
title('Valve 2')
hold on
plot(-5:bins:5, mean(lk1_v2), 'color', c_ord(1,:), 'LineWidth', 2.5);
plot(-5:bins:5, mean(lk2_v2), 'color', c_ord(2,:), 'LineWidth', 2.5);

legend({'lick 1', 'lick 2'})

% plot(-5:bins:5, mean(lk1_v2) + std(lk1_v2)/sqrt(length(lk1_v2)), '--', 'color', [c_ord(1,:) .3]);
% plot(-5:bins:5, mean(lk1_v2) - std(lk1_v2)/sqrt(length(lk1_v2)), '--', 'color', [c_ord(1,:) .3]);

% plot(-5:bins:5, mean(lk2_v2) + std(lk2_v2)/sqrt(length(lk2_v2)), '--', 'color', [c_ord(2,:) .3]);
% plot(-5:bins:5, mean(lk2_v2) - std(lk2_v2)/sqrt(length(lk2_v2)), '--', 'color', [c_ord(2,:), .3]);
vline(0, 'k')
    
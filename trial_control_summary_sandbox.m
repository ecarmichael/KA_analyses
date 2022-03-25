%% try some plots of the events data

%% generate some rates
bins = 0.05;
t = 0:bins:90;

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
plot(t, lk1,  'color', c_ord(1,:))

plot(t, lk2,  'color', c_ord(2,:))
h = bar(tbin_centers,lk_count);
set(h, 'BarWidth', 1)

ylabel('licks/s')
xlabel('time (s)')

linkaxes(ax, 'x')


%% state based per-event

% get the average lick rats for states 1 and 2

% for ii = 1:length(events.t_valve1)
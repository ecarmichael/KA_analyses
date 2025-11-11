function [z_err, R2_out, S_R2_out, plt_mat, plt_s_mat] = KA_lin_decode(cfg_in, x, y, plt_flag)
%% KA_lin_decode: uses the built in fitrlinear function to use linear regression to predict y given x. Uses a 10 fld cross validation and compares it to a shuffle/
if nargin < 4
    plt_flag = 0;
end

cfg_def=[];
cfg_def.hold_out = .3;
cfg_def.xfold = 10;
cfg_def.nShuff = 500;
cfg_def.hold_out = .2;

cfg = ProcessConfig(cfg_def, cfg_in);


%remove NaNs
nan_idx = isnan(x)  | isnan(y);

x(nan_idx) = [];
y(nan_idx) = [];

%% 10 fold xval decoding using a 80/20 split
rng(100, 'twister');
f_err = [];

s_idx = datasample(1:length(x), cfg.nShuff, 'Replace', false);
s_err = [];

plt_mat = NaN(length(y), 10);
plt_s_mat = NaN(length(y), 10);
    cvp = cvpartition(length(y), 'kFold', 10); %HoldOut = cfg.hold_out);

% loop for x val
for ii = 1:10

    %
    t_idx  = training(cvp, ii);
    test_idx = test(cvp, ii);
    idx_n = find(t_idx);

    % disp(idx_n(1:10)')
    % t_idx = 1:floor(length(x)*.7);
    % test_idx = floor(length(x)*.7)+1:length(x);

    [p, S, mu] = polyfit(x(t_idx)', y(t_idx)', 3); % predict y given x

    R2(ii) = 1 - (S.normr/norm(y(t_idx) - mean(y(t_idx))))^2; 

    pred = polyval(p, x(test_idx)', [], mu);

    f_err(ii) = mean((y(test_idx)' - pred).^2);

plt_mat(test_idx, ii) = pred; 

    for kk = 1:length(s_idx) % shuffles

        x_s = circshift(x, s_idx(kk));
        [p, S, mu] = polyfit(x_s(t_idx)', y(t_idx)', 3); % predict y given x

        s_R2(ii,kk) = 1 - (S.normr/norm(y(t_idx) - mean(y(t_idx))))^2; 

        pred_s = polyval(p, x_s(test_idx)', [], mu);

        s_err(ii, kk) = mean((y(test_idx)' - pred_s).^2);
    end
    plt_s_mat(test_idx, ii) = pred_s; 

end

f_mean_err = mean(f_err);
s_mean_err = mean(s_err, 'all');

z_err = (f_mean_err - s_mean_err) / std(s_err,[], 'all');

R2_out = mean(R2); 

S_R2_out = mean(s_R2, 2); 

%% plot if needed

if plt_flag
    fprintf('Mean fold error: %.2f  Vs shuffle err: %.2f\n', mean(f_err), mean(s_err, 'all'))
    figure(909)
    clf
    ax(1) = subplot(4,1,1);
    plot(1:length(y), x, 'b')
    ylabel('FR')

    ax(2) = subplot(4,1,2);
    hold on
    plot(1:length(y), mean(plt_mat, 2, 'omitnan'), '-', 'color', 'r')

    ylabel('predicted spd')
    ylim([0 inf])

    ax(3) = subplot(4,1,3);
    % yyaxis right
    % hold on
    plot(1:length(y), y)
    ylabel('Spd')


    ax(4) = subplot(4,1,4);
    plot(1:length(y), y)
    ylabel('Spd')
    y_lim = ylim;

    yyaxis right
    hold on
    plot(1:length(plt_mat),mean(plt_mat, 2, 'omitnan'), '--')
    plot(1:length(plt_mat), mean(plt_s_mat, 2, 'omitnan'), '--', 'Color', [.7 .7 .7])

    ylabel('predicted spd')

    ylim(y_lim)
    linkaxes(ax, 'x')

end






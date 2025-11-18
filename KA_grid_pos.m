function [pos_grid] = KA_grid_pos(cfg_in, pos)



%% initialize
cfg_def = []; 
cfg_def.plot = 0; 
cfg_def.grid = [86 140; 86 1; 19 68; 147 68]; 
cfg_def.n_points = 100; 


cfg = ProcessConfig2(cfg_def, cfg_in);


%% use GUi to make a grid

if isempty(cfg.grid)

    figure(10011)
    cla
    plot(pos.data(1,:), pos.data(2,:), '.k')

    disp('Draw a line from the Northern feeder to the Southern')
    NS_line_roi = MS_drawline_wait();

    disp('Draw a line from the Western feeder to the Eastern')
    WE_line_roi = MS_drawline_wait();


    cfg.grid = [NS_line_roi.Position; WE_line_roi.Position]; 
end



x = pos.data(1,:); 
y = pos.data(2,:); 

keep_idx = ~isnan(x) & ~isnan(y); 


% make lines between the ends

NS_line = [linspace(cfg.grid(1,1), cfg.grid(2,1), cfg.n_points); linspace(cfg.grid(1,2), cfg.grid(2,2), cfg.n_points)];
WE_line = [linspace(cfg.grid(3,1), cfg.grid(4,1), cfg.n_points); linspace(cfg.grid(3,2), cfg.grid(4,2), cfg.n_points)];

%% snap to the new grid
line_grid = [NS_line, WE_line];
[knn_idx, ~] = knnsearch(line_grid', [x(keep_idx); y(keep_idx)]'); 



x_out = nan(size(x)); 
x_out(keep_idx) = line_grid(1,knn_idx); 

y_out = nan(size(y)); 
y_out(keep_idx) = line_grid(2,knn_idx); 

out = fillmissing2([x_out; y_out], 'nearest'); 


pos_grid = pos; 
pos_grid.data(1,:) = out(1,:); 
pos_grid.data(2,:) = out(2,:); 


if cfg.plot

figure(10010)
cla
hold on
plot(x, y, '.k')
plot(NS_line(1,:), NS_line(2,:), '.r')
plot(line_grid(1,knn_idx), line_grid(2,knn_idx), 'ob')
plot(pos_grid.data(1,:), pos_grid.data(2,:),'xm')
end

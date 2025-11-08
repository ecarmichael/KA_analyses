function [speed_grid,speedVec] = velo_map(velo,nbins)

%compute velocity

maxSpeed = max(velo); 

speedVec = maxSpeed/nbins/2:maxSpeed/nbins:maxSpeed-maxSpeed/nbins/2;
speed_grid = zeros(numel(velo),numel(speedVec));

for i = 1:numel(velo)

    % figure out the speed index
    [~, idx] = min(abs(velo(i)-speedVec));
    speed_grid(i,idx) = 1;
    

end

return
function [t_minus_grid,tminusVec] = tminus_map(t_minus_in,nbins, maxT)


t_minus_in(t_minus_in>maxT) = maxT; %send everything over 50 cm/s to 50 cm/s


tminusVec = maxT/nbins/2:maxT/nbins:maxT-maxT/nbins/2;
t_minus_grid = zeros(numel(t_minus_in),numel(tminusVec));

for i = 1:numel(t_minus_in)

    % figure out the tminus index
    [~, idx] = min(abs(t_minus_in(i)-tminusVec));
    t_minus_grid(i,idx) = 1;
    

end

return
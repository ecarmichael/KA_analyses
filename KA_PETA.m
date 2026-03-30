function [mat_out, tvec] = KA_PETA(S_vec, evt, win)
%% KA_PETA: peri-event-time-average: takes a vector and generates a matrix across trials. 





%% loop over events and take a snippet. 
% remove events that fall outside of the data with the window
s_keep = evt+win(1) > S_vec.tvec(1); 
e_keep = evt+win(2) < S_vec.tvec(end); 

evt(~s_keep | ~e_keep) = []; 

s_idx = nearest_idx3(evt+win(1),S_vec.tvec); 
e_idx = nearest_idx3(evt+win(2),S_vec.tvec); 

mat_out = nan(length(evt), mode(e_idx - s_idx + 1)); % Prep output

% loop over events and get the rate. 
for ii = length(evt):-1:1
    mat_out(ii, :) = S_vec.data(s_idx(ii):e_idx(ii)); % Extract snippet for each event
end

tvec = S_vec.tvec(s_idx(1):e_idx(1)) - S_vec.tvec(s_idx(1)) + win(1); 
function mat_out = KA_PETA(S_vec, evt, win)
%% KA_PETA: peri-event-time-average: takes a vector and generates a matrix across trials. 





%% loop over events and take a snippet. 

s_idx = nearest_idx(evt+win(1),S_vec.tvec); 
e_idx = nearest_idx(evt+win(2),S_vec.tvec); 
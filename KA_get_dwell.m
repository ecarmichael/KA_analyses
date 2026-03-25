function dwell_iv = KA_get_dwell(dwell_tbl, sess_id)
%% KA_get_dwell: looks up the dwell start and stop times in the dwell table from JS. 
%
%
%
%  Rat conversion code:
sub_code= {'R001', 'C1_1'; 'R002', 'C3_2'; 'R003', 'C3_3';  'R004', 'C3_4'; 'R005', 'C4_3'; 'R006', 'C5_2'; 'R007', 'C5_3'; 'R008', 'C6_3'; 'R009', 'C6_4'}; 

% find the rat idx
this_rat = sub_code{contains(sess_id(1:4), sub_code(:,2)),1}; 
this_sess = [this_rat '-' sess_id(9:12) '-' sess_id(14:15) '-' sess_id(17:18)]; 

keep_idx = contains(dwell_tbl.SSN, this_sess); 

% get all the start and stop times as an iv format.

dwell_iv = iv(dwell_tbl.tStart(keep_idx), dwell_tbl.tEnd(keep_idx)); 
dwell_iv.usr = dwell_tbl.lapZone(keep_idx);
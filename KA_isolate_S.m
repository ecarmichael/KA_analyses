function S_out = KA_isolate_S(S, id)
%% KA_isolate_S: isolates one cell of interest from a spike structure 'S'. 
%
%   Inputs:
%     - S: [struct] contains the spike times, labels, ...
%     - id: [string OR index]   label string or index to use. 
%
%
%  Outputs: 
%     - S_out: [struct]   Spike structure with all but the cell of
%     interest removed. 
%
%% determine the id format and remove all others accordingly. 

if isnumeric(id) && length(S.t) <= id
    fprintf('<strong>%s:</strong> input was numeric (%0.0d) = %s\n', mfilename, id, S.label{id})
    idx = id;
elseif ischar(id)
    fprintf('<strong>%s:</strong> input was a string: %s\n', mfilename, id)
    idx = (ismember(S.label, id));
end

%% remove extra cells from fields from data and label fields. 
S_out = S; 

f_list = fieldnames(S); 

for ii = length(f_list):-1:1
    
    if length(S_out.(f_list{ii})) == length(S.t)
        S_out.(f_list{ii})(~idx) = [];
    end
    
end



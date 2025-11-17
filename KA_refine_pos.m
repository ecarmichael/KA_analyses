function [pos_out] = KA_refine_pos(cfg_in, pos_in)
%%  KA_refine_pos: takes the position data and applies spatial limits (ie regions where the animal could not have left) and intepolates to fill the missing pieces

% inputs
%   cfg_in: [struct] contains parameters
%         - limits: cell array with x -y pairs; 
%
%
%    pos_in: [strcut] data in the tsd position format (pos.tvec, pos.data,...)
%
%
% Outputs:
%
%   pos_out: [struct] refined position datat structure in the same format
%   as pos_in


%% initialize


cfg_def = []; 
cfg_def.limits = [70 105; 40 90]; 
cfg_def.method = 'nearest';

cfg = ProcessConfig2(cfg_def, cfg_in); 

% S_pix_x = [70 100]; S_pix_y = 130;
% W_pix_x = 25;     W_pix_y = [55, 85];
% N_pix_x = [70 100]; N_pix_y = 10;
% E_pix_x = 145;      E_pix_y = [55, 85];
%% apply limits to block out data
rm_idx = zeros(size(pos_in.tvec)); 


if ~isempty(cfg.limits)
   for ii = length(pos_in.data):-1:1
        if (pos_in.data(1,ii) < cfg.limits(1,1))  && (pos_in.data(2,ii) < cfg.limits(2,1))
            rm_idx(ii) = 1; 
        elseif (pos_in.data(1,ii) < cfg.limits(1,1))  && (pos_in.data(2,ii) > cfg.limits(2,2))
                rm_idx(ii) = 1;
        elseif (pos_in.data(1,ii) > cfg.limits(1,2))  && (pos_in.data(2,ii) < cfg.limits(2,1))
                    rm_idx(ii) = 1;
        elseif (pos_in.data(1,ii) > cfg.limits(1,2))  && (pos_in.data(2,ii) > cfg.limits(2,2))
                        rm_idx(ii) = 1;
        end

   end
end

rm_idx = logical(rm_idx);

%% interpolate over the removed data

pos_out = pos_in; 

pos_out.data(:, rm_idx) = NaN;

pos_out.data = fillmissing(pos_out.data, cfg.method, 'EndValues', cfg.method);






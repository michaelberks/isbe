function [v_spatial, v_temporal] = varnames()
% Return a list of variable names

v_spatial  = {'Q', 'uu', 'vv', 'F', 'theta'};
% v_spatial  = {'Q', 'F', 'theta'};
% v_spatial  = {'Q', 'uu', 'vv'};

v_temporal = {'f','displacements'};

% Concatenate variable names if only one return value requested
if (nargout == 1)
    v_spatial = [v_spatial v_temporal];
end

function [lb, ub] = varlimits(varname)
% Return the upper and lower bounds of a variable

switch varname
    case {'Q'},
        % Vessel occupancy.
        lb = 0; 
        ub = 1;
        
    case {'F','uu','vv'}, 
        % Spatial flow is defined as having a maximum of 1.
        lb = -1; 
        ub = 1;
        
    case {'theta','f'},
        % The optimizer is unlikely to cycle through multiple periods of
        % theta.
        % No limits on temporal flow.
        lb = -inf;
        ub = inf;
        
    otherwise,
        lb = -inf;
        ub = inf;
end

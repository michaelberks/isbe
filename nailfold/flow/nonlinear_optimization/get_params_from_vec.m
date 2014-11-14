function [p, wts, stackSize, uu, vv] = ...
    get_params_from_vec(observations, variables)

p = merge_structs(observations, variables);
p = reshape_variables(p);

wts = p.Q;

stackSize = size(p.imgStack);

if isfield(p, 'uu')
    uu = p.uu;
    vv = p.vv;
else
    uu = p.F .* cos(p.theta);
    vv = p.F .* sin(p.theta);
end

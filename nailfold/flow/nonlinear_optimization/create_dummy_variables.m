function [variables, observations, var_vec, sz_vec] = ...
    create_dummy_variables(f_use_uv)

if ~exist('f_use_uv','var'), f_use_uv = true; end

imgStack = rand(3,3,8);
imsz = size(imgStack);

nObservations = numel(imgStack)+1;
observationMask = create_observation_mask(size(imgStack), nObservations);

gt.mask = ones(size(imgStack(:,:,1)));
gt.presence = [];
gt.orientation = [];

displacements = rand(size(imgStack,3)-1,2);

uu = zeros(size(gt.mask));
vv = uu;
FF = zeros(size(gt.mask));
tt = zeros(size(gt.mask));



if f_use_uv
    [variables, observations] = define_variables_and_observations(...
        imgStack, imgStack, gt, observationMask, ...
        f_use_uv, uu, vv, displacements);
else
    [variables, observations] = define_variables_and_observations(...
        imgStack, imgStack, gt, observationMask, ...
        f_use_uv, FF, tt, displacements);
end
    
[var_vec, sz_vec, lbound, ubound] = ncm_pack(variables);
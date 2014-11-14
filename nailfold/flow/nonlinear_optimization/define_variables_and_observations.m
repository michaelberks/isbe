function [variables, observations] = ...
    define_variables_and_observations(imgStack, imgStackWarped, ...
                                      gt, observationMask, ...
                                      f_use_uv, flow1, flow2, displacements)
% Create two structures that will hold the parameters over which we'll
% optimize (variables) and the parameters that will remain fixed
% (observations).
% In practice, everything is present in observations but the variables take
% precedent during the optimization.
                                  
if ~exist('displacements','var'), displacements = []; end
    
% Displacements with respect to first image.
nImages = size(imgStack,3);
if isempty(displacements), displacements = zeros(nImages-1,2); end

% These will always be observations
observations.imgStack = imgStack;
observations.imgStackWarped = imgStackWarped;
observations.penalty = 0;
observations.obsMask = observationMask;

% These can be either observations or variables
if f_use_uv
    observations.uu = flow1;
    observations.vv = flow2;
else
    observations.F = flow1;
    observations.theta = flow2;
end
observations.Q = 0.5 * ones(size(gt.mask));
observations.f = ones([nImages,1]);
observations.displacements = displacements;
observations.orientation = gt.orientation;
observations.presence = gt.presence;

% Comment/uncomment these lines to define the variables (uncommented)

% Overwrite observations with these variables when the parameter vector 
% is 'ncm_pack'ed
if f_use_uv
    variables.uu = observations.uu;
    variables.vv = observations.vv;
else
    variables.F = observations.F;
    variables.theta = observations.theta;
end
% variables.Q = observations.Q;
% variables.f = observations.f;
% variables.displacements = observations.displacements;
% variables.orientation = observations.orientation;
% variables.presence = observations.presence;

function [target_pos] = select_corresponding_position(source_breast, target_breast, source_pos, flip)
%SELECT_CORRESPONDING_POSITION given a source and target breast shape and a
%source position, finds a corresponding position in the target breast
%   [target_pos] = select_corresponding_position(source_breast, target_breast, source_pos, flip)
%
% Inputs:
%      source_breast - N*2 array of xy coordinates specifying shape of
%      source breast
%
%      target_breast - N*2 array of xy coordinates specifying shape of
%      target breast
%
%      source_pos - xy coordinates of position
%
%      flip - logical value specifying whether or not to horizontally flip
%      the source breast prior to finding correspondence. Set to 1 if
%      matching left/right pairs or 0 if match temporal pairs of the same
%      breast
%
%
% Outputs:
%      target_pos - xy coordinates of target position
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 22-Jun-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

%Check source and target shapes have same number of points
if size(source_breast,1) ~= size(target_breast,1)
    cum_dist = [0; cumsum(sqrt(sum(diff(source_breast).^2,2)))];
    source_breast = interp1(cum_dist, source_breast, linspace(0, cum_dist(end), size(target_breast,1)), 'linear');
    
end

%Check whether we need to flip the source breast shape and position
if flip
    source_breast(:,1) = -source_breast(:,1);
    source_pos(:,1) = -source_pos(:,1);
end

%Procrustes align the source shape to the target shape
[dummy dummy t] = mb_procrustes(target_breast, source_breast);

%Map source position into the space of the target breast shape
target_pos = t.b*source_pos*t.T + repmat(t.c(1,:), size(source_pos,1),1);
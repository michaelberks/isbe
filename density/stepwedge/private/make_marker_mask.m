function [marker_mask] = make_marker_mask(IMAGE, x_marker, y_marker, r_marker, resize_factor)
%MAKE_MARKER_MASK *Insert a one line summary here*
%   [breast_border,breast_air,errorcheck] = make_marker_mask(IMAGE,segmentation)
%
% Inputs:
%      IMAGE- *Insert description of input variable here*
%
%      x_marker- X-coordinates of markers
%      y_marker- Y-coordinates of markers
%      r_marker- radii of markers
%
%
% Outputs:
%      marker_mask- BW mask, same size as IMAGE, that masks the markers so
%      these won't influence the density calculations
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Sep-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
dimensions = size(IMAGE);


% make a mask with ones inside marker and zeros outside - we'll use
% this later to mask the breast so the bright white markers don't
% influence the density calculations
marker_mask = false(dimensions);
xx = repmat(1:dimensions(:,2), dimensions(1), 1);
yy = repmat((1:dimensions(:,1))', 1, dimensions(2));

%Resize the marker co-ordinates and radii to new image size
x_marker = (x_marker-1)*resize_factor + 1;
y_marker = (y_marker-1)*resize_factor + 1;
r_marker = r_marker * 1.1 * resize_factor;

for ii = 1:numel(x_marker)
    marker_mask = marker_mask | ...
        ((xx-x_marker(ii)).^2 + (yy-y_marker(ii)).^2 < (r_marker(ii)+5)^2);
end
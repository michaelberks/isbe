function [discard_pts, vessel_centre] = ...
    discard_edge_preds(vessel_centre, fov_mask, border_width)
%DISCARD_EDGE_PREDS Discard centre lines near the edge of the image
%   [] = extract_vessel_centres()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('border_width', 'var') || isempty(border_width)
    border_width = 32;
end
[nrows ncols] = size(fov_mask);

%Make a mask of centre pixels from the vessel centre struc
vessel_centre_idx = sub2ind([nrows ncols], vessel_centre.y, vessel_centre.x);

%Make mask of the border and the orientation of the its edges
border_mask = fov_mask & ~imerode(fov_mask, strel('disk', border_width));
discard_pts = border_mask(vessel_centre_idx);
    
if nargout == 2
    f = fieldnames(vessel_centre);
    for i_f = 1:length(f)
        vessel_centre.(f{i_f})(discard_pts) = [];
    end
end
    

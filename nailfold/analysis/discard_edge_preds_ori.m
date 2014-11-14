function [discard_pts] = ...
    discard_edge_preds_ori(vessel_centre, fov_mask, border_width, ori_thresh)
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
if ~exist('ori_thresh', 'var') || isempty(ori_thresh)
    ori_thresh = pi*10/180;
end
[nrows ncols] = size(fov_mask);

%Make a mask of centre pixels from the vessel centre struc
vessel_centre_mask = false(nrows, ncols);
vessel_centre_idx = sub2ind([nrows ncols], vessel_centre.y, vessel_centre.x);
vessel_centre_mask(vessel_centre_idx) = 1;

%Make mask of the border and the orientation of the its edges
border_mask = fov_mask & ~imerode(fov_mask, strel('disk', border_width));
border_centre = border_mask(vessel_centre_mask);
[~, border_ori] = gaussian_1st_derivative_gradient(fov_mask, 16);   
border_ori = exp(2i*border_ori(vessel_centre_mask));

%Get labelling of connected components in the vessel centre mask
centre_labels = bwlabel(vessel_centre_mask, 8);
centre_labels = centre_labels(vessel_centre_mask);

%Go through each component and discard it if aligns with an edge
discard_pts = false(size(vessel_centre.x));
for i_lab = 1:max(centre_labels(:))
    label_mask = centre_labels == i_lab;

    label_ori = vessel_centre.ori(label_mask);        
    label_border = border_centre(label_mask);
    label_border_ori = border_ori(label_mask);

    %Pixels near the border, with matching orientation - probably picke
    %dup edge sof the image frames
    if any(label_border) && ...
        median(abs(ori_error(label_border_ori(label_border), label_ori(label_border)))) < ori_thresh
        discard_pts(label_mask) = 1;
    end
end

% discard_pts = border_centre &...
%     abs(ori_error(border_ori, vessel_centre.ori)) < ori_thresh;
    
    

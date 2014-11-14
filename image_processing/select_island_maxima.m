function [maxima_xy, maxima_vals, discarded_xy, discarded_vals] = select_island_maxima(maxima_xy, maxima_vals, mask)
%LOCAL_IMAGE_MAXIMA compute the positions of local maxima within an image
%   [maxima_pos, maxima_vals] = local_image_maxima(image_in, exclusion_zone, mask, threshold)
%
% Inputs:
%      image_in - image (or transformation of image) in whihc to find local
%      maxima
%
%      exclusion_zone - scalar value defining radius about each maxima that
%      will exclude other maxima (the largest local maxima will be included)
%
%      mask - (optional) mask of pixels to consider in the image (if not
%      set all pixels will be considered)
%
%      threshold - (optional) scalar value that if sets allows only those
%      pixels with greater than threshold to count as maxima
%
%
% Outputs:
%      maxima_pos - Nx2 array giving the x/y coordinates of each local
%      maxima
%
%      maxima_vals - Nx1 vector of image values at the local maxima
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 11-Oct-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
connected_label_mask = bwlabel(mask);

%Sort the maxima by value
[maxima_vals sorted_idx] = sort(maxima_vals, 'descend');
maxima_xy = maxima_xy(sorted_idx,:);

%Get the labels for each maxima
island_labels = interp2(connected_label_mask, maxima_xy(:,1), maxima_xy(:,2), 'nearest');
 
%1) Get the maximum val and position of each island
num_maxima = size(maxima_xy,1);
discard_idx = false(num_maxima,1);
for i_m = 1:num_maxima

    if ~discard_idx(i_m); 
        %Throw away all the other maxima on this island
        discard_idx(island_labels==island_labels(i_m)) = 1;
        discard_idx(i_m) = 0; %Don't throw yourself out!
    end
end

discarded_xy = maxima_xy(discard_idx,:);
discarded_vals = maxima_vals(discard_idx,:);
        
maxima_vals(discard_idx,:) = [];
maxima_xy(discard_idx,:) = [];
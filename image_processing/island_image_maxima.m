function [maxima_pos, maxima_vals] = island_image_maxima(image_in, island_threshold, exclusion_zone, mask, threshold, do_plot)
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
 [r c] = size(image_in);
 if ~exist('exclusion_zone', 'var') || isempty(exclusion_zone)
     exclusion_zone = 0;
 end
 if ~exist('mask', 'var') || isempty(mask)
     mask = true(r, c);
 end
 if ~exist('threshold', 'var')
     threshold = -inf;
 end
 if ~exist('do_plot', 'var')
     do_plot = 0;
 end
 mask = mask & (image_in > threshold);
 
 connected_label_mask = bwlabel(image_in > island_threshold);
 connected_label_mask(~mask) = 0;
 
 
 %Check we've got at least some points to work with, otherwise bale out now
 if ~any(connected_label_mask(:))
     maxima_pos = []; 
     maxima_vals = [];
     return;
 end
 
 %1) Get the maximum val and position of each island
 num_maxima = max(connected_label_mask(:));
 discard_idx = false(num_maxima,1);
 maxima_x = zeros(num_maxima,1);
 maxima_y = zeros(num_maxima,1);
 maxima_vals = zeros(num_maxima,1);
 for i_m = 1:num_maxima
     
     island_mask = connected_label_mask == i_m;
     
     if  ~any(island_mask(:))
         %This label got wiped out by the inclusion mask
         discard_idx(i_m) = 1;
     else
         [island_y island_x] = find(island_mask);
         island_vals = image_in(island_mask);
         [maxima_vals(i_m) max_idx] = max(island_vals);
         maxima_x(i_m) = island_x(max_idx);
         maxima_y(i_m) = island_y(max_idx);
     end
 end
 maxima_vals(discard_idx,:) = [];
 maxima_x(discard_idx,:) = [];
 maxima_y(discard_idx,:) = [];
 
 %Sort the maxima by value
 [sorted_maxima sorted_idx] = sort(maxima_vals, 'descend');
 sorted_x = maxima_x(sorted_idx);
 sorted_y = maxima_y(sorted_idx);
 
 %Go through each maxima and discard any other that lie within its
 %exclusion zone
 num_maxima = length(sorted_maxima);
 keep_idx = true(num_maxima, 1);
 for ii = 1:num_maxima
     
     %Don't process pts we've already exlcuded
     if keep_idx(ii)
         
         %get x/y coordinates of current maxima
         x = sorted_x(ii);
         y = sorted_y(ii);
         
         %Only keep pts which are further away than the exlcusion zone radius
         keep_idx = keep_idx & ((sorted_x-x).^2 + (sorted_y-y).^2 > exclusion_zone^2);
         
         %Add current keep idx back in (it will be excluded above)
         keep_idx(ii) = true;
     end
 end
 
 %Return the x/y coordinates of the maxima
 maxima_pos = [sorted_x(keep_idx) sorted_y(keep_idx)];
 
 %Return the image values at the maxima if necessary
 if nargout > 1
     maxima_vals = sorted_maxima(keep_idx);
 end
 
 if do_plot
     figure; imgray(image_in);
     plot(maxima_pos(:,1), maxima_pos(:,2), 'x');
 end
     
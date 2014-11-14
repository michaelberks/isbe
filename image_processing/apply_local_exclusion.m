function [maxima_pos, maxima_vals, keep_idx] = apply_local_exclusion(maxima_pos, maxima_vals, exclusion_zone)
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
 if nargin < 2
     exclusion_zone = 1;
 end
 
 %Sort the maxima by value in the original image
 [sorted_maxima sorted_idx] = sort(maxima_vals, 'descend');
 sorted_x = maxima_pos(sorted_idx,1);
 sorted_y = maxima_pos(sorted_idx,2);
 
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
     
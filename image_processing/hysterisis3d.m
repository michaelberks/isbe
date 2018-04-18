function [BW] = hysterisis3d(magnitudes,thresh,percent_not_edges,thresh_ratio,apply_skel)
%HYSTERISIS *Insert a one line summary here*
%   [BW] = hysterisis(magnitudes,thresh)
%
% Inputs:
%      magnitudes - *Insert description of input variable here*
%
%      thresh - *Insert description of input variable here*
%
%       percent_not_edges -
%
%       thresh_ratio -
%
% Outputs:
%      BW- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 07-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

if nargin < 2
    thresh = [];
end
if nargin < 3
    percent_not_edges = .99;
end
if nargin < 4
    thresh_ratio = .7;
end   
if nargin < 5
    apply_skel = false;
end
magnitudes = magnitudes / max(magnitudes(:));

% Select the thresholds
    if isempty(thresh) 
        high_thresh = prctile(magnitudes(:), 100*percent_not_edges);
        low_thresh = thresh_ratio*high_thresh;

    elseif length(thresh)==1
        high_thresh = thresh;
        if thresh >= 1
            eid = sprintf('Images:%s:thresholdMustBeLessThanOne', mfilename);
            msg = 'The threshold must be less than 1.'; 
            error(eid,'%s',msg);
        end
        low_thresh = thresh_ratio*thresh;

    elseif length(thresh) == 2
        low_thresh = thresh(1);
        high_thresh = thresh(2);
        
        if (low_thresh >= high_thresh) || (high_thresh >= 1)
            eid = sprintf('Images:%s:thresholdOutOfRange', mfilename);
            msg = 'Thresh must be [low high], where low < high < 1.'; 
            error(eid,'%s',msg);
        end
    end
    
    %Now do the edge hysterisis
    %Compute the strong edges and the weak edges
    weak_edges = magnitudes > low_thresh;
    strong_edges = magnitudes > high_thresh;
    
    %Now find weak edges that are 8-connected to a strong edge
    [n_y, n_x, n_z] = size(strong_edges);
    BW = false(n_y, n_x, n_z);
    for i_z = 1:n_z
        strong_edges_i = strong_edges(:,:,i_z);
        if i_z > 1
            strong_edges_i = strong_edges_i | strong_edges(:,:,i_z-1);
        end
        if i_z < n_z
            strong_edges_i = strong_edges_i | strong_edges(:,:,i_z+1);
        end
        
        if any(strong_edges(:))
            [rstrong, cstrong] = find(strong_edges_i);
            combined_edges = ...
                bwselect(weak_edges(:,:,i_z), cstrong, rstrong, 8);
        else
            combined_edges = weak_edges(:,:,i_z);
        end
        %finally, thin the combined edges
        if apply_skel
            BW(:,:,i_z) = bwmorph(combined_edges, 'skel', inf);
        else
            BW(:,:,i_z) = combined_edges;
        end
    end
    
        
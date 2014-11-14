function [BW] = hysterisis(magnitudes,thresh,percent_not_edges,thresh_ratio)
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
[m n] = size(magnitudes);
magnitudes = magnitudes / max(magnitudes(:));

% Select the thresholds
    if isempty(thresh) 
        counts = imhist(magnitudes, 64);
        high_thresh = ...
            find(cumsum(counts) > percent_not_edges*m*n, 1,'first') / 64;
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
    if any(strong_edges(:))
        [rstrong cstrong] = find(strong_edges);
        combined_edges = bwselect(weak_edges, cstrong, rstrong, 8);
    else
        combined_edges = weak_edges;
    end
    
    %finally, thin the combined edges
    BW = bwmorph(combined_edges, 'skel', inf);
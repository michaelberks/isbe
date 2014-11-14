function [BW] = ilp_hysterisis(ilp_band, mag_thresh, phase_thresh)
%HYSTERISIS *Insert a one line summary here*
%   [BW] = hysterisis(magnitudes,thresh)
%
% Inputs:
%      magnitudes - *Insert description of input variable here*
%
%      mag_thresh - *Insert description of input variable here*
%
%      phase_thresh -
%
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

magnitudes = abs(ilp_band);
phase = angle(ilp_band);
magnitudes = magnitudes / max(magnitudes(:));

%Now do the edge hysterisis
%Compute the strong edges and the weak edges
weak_lines = (magnitudes > mag_thresh(1)) & (phase > phase_thresh(1));
strong_lines = (magnitudes > mag_thresh(2)) & (phase > phase_thresh(2));

%Now find weak edges that are 8-connected to a strong edge
if any(strong_lines(:))
    [rstrong cstrong] = find(strong_lines);
    combined_lines = bwselect(weak_lines, cstrong, rstrong, 8);
else
    combined_lines = weak_lines;
end

%finally, thin the combined edges
%BW = bwmorph(combined_lines, 'skel', inf);
%BW = combined_lines;
BW = bwmorph(combined_lines, 'thin', 1);
function [mass_centre] = sample_mass_location(location_model,segmentation, thresh, debug_mode)
%SAMPLE_MASS_LOCATION *Insert a one line summary here*
%   [mass_centre] = sample_mass_location(location_model,segmentation,debug_mode)
%
% Inputs:
%      location_model- *Insert description of input variable here*
%
%      segmentation- *Insert description of input variable here*
%
%      debug_mode- *Insert description of input variable here*
%
%
% Outputs:
%      mass_centre- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

if nargin < 4
    debug_mode = 1;
end

%Get a sorted list of the density distribution
[D_sort D_idx] = sort(location_model.D_a);

%Take cumulative sum
cdf = cumsum(D_sort);

%Sample a rand threshold from which to select a point
if nargin < 3 || isempty(thresh)
    thresh = rand;
end

%Work out the points with density higher than thresh    
selected = D_idx(cdf > thresh);

%Randomly sample one of these points as the mass centre
samp_idx = selected(ceil(rand * length(selected)));
mass_centre = [location_model.x(samp_idx) location_model.y(samp_idx)];

if debug_mode
    figure; axis ij equal; hold on;
    plot(location_model.mean_shape(:,1), location_model.mean_shape(:,2));
    plot(location_model.x(selected), location_model.y(selected), 'r.');
    plot(mass_centre(1), mass_centre(2), 'y*', 'MarkerSize', 10);
end

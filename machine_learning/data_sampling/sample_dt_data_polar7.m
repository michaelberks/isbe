function [dt_samples] = sample_dt_polar7(dt, rows, cols, rotate)
%SAMPLE_DT_POLAR7 *Insert a one line summary here*
%   [dt_samples] = sample_dt_polar7(dt, rows, cols, rotate)
%
% Inputs:
%      dt - Dual-tree decomposition from which to sample
%
%      rows - row subscripts of sample locations
%
%      cols - column subscripts of sample locations
%
%      rotate - if true (default false) circular shift the columns of the
%      final representation so maximal orientation is in first row, thus
%      achieving approximate rotation invariance (+/- 15 degrees)
%
%
% Outputs:
%      dt_samples - *Insert description of input variable here*
%
%
% Example:
%
% Notes: Representation is based on that presented by Kingsbury in
% "Rotation-Invariant Local Feature Matching with Complex Wavelets"
%
% See also: DT_TO_PIXEL_SUBSET
%
% Created: 23-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Set default value for rotation
if nargin < 4
    rotate = false;
end


%Count number of samples
num_samples = length(rows);

%Create coordinate points at 12 clock positions about each sample location - 
%first point of clock is at 6 o'clock, subsequent points move anticlockwise
rr = [rows bsxfun(@plus, rows, sin(2*(0:-1:-5)*pi/6))];
cc = [cols bsxfun(@plus, cols, cos(2*(0:-1:-5)*pi/6))];

%Extract interpolated DT coefficients at each sample point and count number
%of levels
if isstruct(dt)
    %Use NAG interpolation
    samples = dt_interp_nag_out(...
        dt.knot_mag, dt.knot_im, dt.knot_re, dt.dt_dims, rr, cc);
    num_levels = length(dt.knot_mag);
else
    %Use matlab interpolation
    samples = dt_to_pixel_subset(dt, rr, cc);
    num_levels = length(dt) - 1;
end

%Pre-allocate output matrix
dt_samples = zeros(num_samples,6,7,num_levels);

%Copy DT coeffs from the 12 clock positions at each location, circular
%shifting to achieve rotation invariant representation
for ii = 1:6;
    dt_samples(:,ii,2:7,:) = circshift(samples(:,1+ii,:,:), [0 0 2*(1-ii) 0]);
end

%Swap 2nd and 3rd dimensions, so 6 subbands are across 3rd dim
dt_samples = permute(dt_samples, [1 3 2 4]);

%Copy DT coeffs from clock centre of each sample location
dt_samples(:,1,:,:) = samples(:,1,:,:);

%Check if we need to rotate by orientation
if rotate
    %Get magnitudes of each band at centre of clock
    centre_mags = reshape(abs(dt_samples(:,1,:,:)), num_samples, []);
    
    %Find the sub-band that corresponds to the maximum magnitude in
    %each row
    [dummy max_ori] = max(centre_mags, [], 2);
    max_ori = rem(max_ori-1,6)+1;
    
    %Circular shift the data according to the maximum subband
    for ori = 1:6
        shift_idx = max_ori == ori;
        dt_samples(shift_idx,:,:,:) = circshift(dt_samples(shift_idx,:,:,:), [0 0 1-ori 0]);
    end 
end
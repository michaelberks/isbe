function [dt_samples] = sample_dt_data_rectangle(dt, rows, cols, varargin)
%SAMPLE_DT_RECTANGLE *Insert a one line summary here*
%   [dt_samples] = sample_dt_polar13(dt, rows, cols, rotate)
%
% Inputs:
%      dt - Dual-tree decomposition from which to sample
%
%      rows - row subscripts of sample locations
%
%      cols - column subscripts of sample locations
%
%      win_size - dimensions of window to sample
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
args = u_packargs(varargin, '0',...
    'levels', 1:5,...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 3,...
    'interpmethod', 'cubic',...
    'band_fequencies', 1i*[-3 -1; -sqrt(5) -sqrt(5); -1 -3; 1 -3; sqrt(5) -sqrt(5); 3 -1]*pi/2.15,...
    'unwrap_phase', 1,...
    'interp_mag_phase', 1,...
    'correct_phase', 1);

clear varargin;

levels = args.levels;
win_size = args.win_size;

%Make copies of sample rows and cols at positions of local window patch
pad_w = floor(win_size/2);
win_idx = -pad_w:pad_w;
num_samples = length(rows);
rr = repmat(rows*ones(1,win_size) + ones(num_samples,1)*win_idx, 1, win_size);
cc = kron(cols*ones(1,win_size) + ones(num_samples,1)*win_idx, ones(1,win_size));

%Extract interpolated DT coefficients at each sample point
if isstruct(dt)
    %Use NAG interpolation
    dt_samples = dt_interp_nag_out(...
        dt.knot_mag, dt.knot_im, dt.knot_re, dt.dt_dims, rr, cc, levels);
else
    %Use matlab interpolation
    interp_args.levels = args.levels;
    interp_args.interpmethod = args.interpmethod;
    interp_args.band_fequencies = args.band_fequencies;
    interp_args.unwrap_phase = args.unwrap_phase;
    interp_args.interp_mag_phase = args.interp_mag_phase;
    interp_args.correct_phase = args.correct_phase;

    dt_samples = dt_to_pixel_subset(dt, rr, cc, interp_args);
end
clear dt rr cc;

%Check if we need to rotate by orientation
if (args.rotate == 2) || (args.do_max == 2)
    %***********OLD CODE, may remove? *********************************
    %Get magnitudes at each orientation/level of centre pixel in window
    centre_mags = reshape(abs(dt_samples(:,ceil(win_size^2 / 2),:,:)), num_samples, []);
    
    %Find the sub-band that corresponds to the maximum magnitude in
    %each row
    [dummy max_ori] = max(centre_mags, [], 2); %#ok
    
    %Circular shift the data according to the maximum subband
    for ori = 0:5
        shift_idx = rem(max_ori,6) == ori;
        dt_samples(shift_idx,:,:,:) = circshift(dt_samples(shift_idx,:,:,:), [0 0 1-ori 0]);
    end
    
    %If we taking the maximum over orientation we can now discard all but
    %the first band
    if args.do_max
        dt_samples(shift_idx,:,2:end,:) = [];
    end
elseif (args.rotate == 1) || (args.do_max == 1)
    %Need to reorder the dimensions so the orientation bands are in the last column,
    % then swap them back after the reordering function has been called
    dt_samples = permute(dt_samples, [1 2 4 3]);
    dt_samples = reorder_bands(dt_samples, args.do_max);
    dt_samples = permute(dt_samples, [1 2 4 3]);
end
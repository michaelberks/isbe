function [linop_samples] = sample_linop_data(image_in, rows, cols, varargin)
%SAMPLE_LINOP_DATA *Insert a one line summary here*
%   [dt_samples] = sample_dt_polar13(dt, rows, cols, rotate)
%
% Inputs:
%      image_in - image to decompose and sample
%
%      num_levels - num_levels to include
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
%      linop_samples - *Insert description of input variable here*
%
%
% Example:
%
% Notes: Representation is based on that presented by Kingsbury in
% "Rotation-Invariant Local Feature Matching with Complex Wavelets"
%
% See also:
%
% Created: 23-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    'num_levels', 5,...
    'num_angles', 8,...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 3,...
    'pca', []);

win_size = args.win_size;
num_angles = args.num_angles;
num_levels = args.num_levels;

[rr, cc, num_samples] = get_indices_from_winsize(rows, cols, win_size);

%Now get linop samples from the image
linop_samples = line_operator_octave_subset(image_in, num_angles, ...
                                            num_levels, rr, cc); 
clear image_in rr cc;

%Check if we need to rotate by orientation
if args.rotate
    %Get magnitudes at each orientation/level of centre pixel in window
    centre_mags = reshape(linop_samples(:,ceil(win_size^2 / 2),:,:), num_samples, []);
    
    %Find the sub-band that corresponds to the maximum magnitude in
    %each row
    [dummy max_ori] = max(centre_mags, [], 2);
    max_ori = rem(max_ori-1,num_angles)+1;
    
    %Circular shift the data according to the maximum subband
    for ori = 1:num_angles
        shift_idx = max_ori == ori;
        linop_samples(shift_idx,:,:,:) = circshift(linop_samples(shift_idx,:,:,:), [0 0 1-ori 0]);
    end 
elseif args.do_max
    %get the maximum response across orientations
    linop_samples = squeeze(max(linop_samples, [], 3));
end

linop_samples = reshape(linop_samples, num_samples, []);
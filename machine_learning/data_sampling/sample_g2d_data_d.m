function [g2d_samples] = sample_g2d_data_d(g2d_responses, rows, cols, varargin)
%SAMPLE_G2D_DATA_D *Insert a one line summary here*
%   [dt_samples] = sample_g2d_data_d(g2d_responses, rows, cols, rotate)
%
% Inputs:
%      g2d_responses - Pyramid of Gaussian 2nd derivative responses from which to sample
%
%      rows - row subscripts of sample locations
%
%      cols - column subscripts of sample locations
%
%      win_size - dimensions of window to sample
%
%
% Outputs:
%      g2d_samples - *Insert description of input variable here*
%
%
% Example:
%
% Notes: Representation is based on that presented by Kingsbury in
% "Rotation-Invariant Local Feature Matching with Complex Wavelets"
%
% See also: SAMPLE_G2D_DATA
%
% Created: 23-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    'interpmethod', 'cubic',...
    'win_size', 3,...
    'pca', []);

win_size = args.win_size; 

%Make copies of sample rows and cols at positions of local window patch
pad_w = floor(win_size/2);
win_idx = -pad_w:pad_w;
num_samples = length(rows);
yi = repmat(rows*ones(1,win_size) + ones(num_samples,1)*win_idx, 1, win_size);
xi = kron(cols*ones(1,win_size) + ones(num_samples,1)*win_idx, ones(1,win_size));

%Get number of levels to put in full tree
num_levels = length(g2d_responses);
[r1 c1 num_bands] = size(g2d_responses{1});

%Pre-allocate space for g2d_samples - can be an array since each level has
%same number of coeffs
g2d_samples = zeros(num_samples, win_size^2, num_bands, num_levels);

%Copy coefficients from level 1 - don't need to interpolate
band_idx = sub2ind([r1 c1], yi, xi);
for band = 1:num_bands
    offset = (band-1)*r1*c1;
    g2d_samples(:,:,band,1) = g2d_responses{1}(band_idx+offset);
end

%for each sub-band from the 2nd level down, interpolate up the coefficients
%to the full size
for lev = 2:num_levels
    
    %Generate x,y coordinates of sub-sampled pixels - take account of
    %padding
    offset = 0.5 - 2^(lev-2);
    sub_x = (0:size(g2d_responses{lev},2)+1) * 2^(lev-1) + offset;
    sub_y = (0:size(g2d_responses{lev},1)+1)' * 2^(lev-1) + offset;
    
    %Interpolate responses of each sub-band
    for band = 1:num_bands
        band_image = padarray(g2d_responses{lev}(:,:,band), [1 1], 'replicate');
        g2d_samples(:,:,band,lev) =...
            interp2(sub_x, sub_y, band_image, xi, yi, args.interpmethod);
    end
end
if any(isnan(g2d_samples(:)))
    display('Stop!');
end
g2d_samples = g2d_samples(:,:);
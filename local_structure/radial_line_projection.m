function [angle_bands] = radial_line_projection(prob_map, orientation_map, num_angles, dist_filt)
%RADIAL_LINE_PROJECTION *Insert a one line summary here*
%   [] = radial_line_projection()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 26-May-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

if nargin < 3 || rem(num_angles(1), num_angles(2))
    num_angles = [36 12];
end

dist_flag = (nargin > 3) && ~isempty(dist_filt);

%Get size of regions
[row col] = size(prob_map);

%Workout padding needed to allow rotation
pad_size = ceil(0.5*(sqrt(2)-1)*max([row col]));

%Pre-allocation outputs
line_sum = zeros(row, col);
angle_bands = zeros(row, col, num_angles(2));

% if dist_flag
%     if nargin < 4
%         %Generate filter used to sum, whilst dividing by distance
%         dist_filt = ones(1, 200);
%     end
%     dist_sum = zeros(row, col);
% end

%Compute angular resolution
ang_res = 180 / num_angles(1);

%Compute number of angles in each band and initialise band counter
band_size = num_angles(1) / num_angles(2);
band = 1;

%For each angle
for ii = 1:num_angles

    %Compute theta
    theta = (ii - 0.5)*ang_res;
    
    %Get mask of pixels that have orientation within theta range
    theta_mask = (orientation_map > theta - 0.5*ang_res) & (orientation_map <= theta+ 0.5*ang_res);
    
    %Mask out all but selected pixels and pad image for rotation
    pad_im = padarray(prob_map .* theta_mask, [pad_size pad_size], 0);
    
    %Rotate image by -theta, so selected pixels project horizontally
    rot_im = imrotate(pad_im, -theta, 'bilinear', 'crop');
    
    if dist_flag
        %To weight by distance, convolve with dist filter across columns rather
        %than simply summing
        sum_im = conv2(rot_im, dist_filt, 'same');
    else
        %Sum horizontally to get projected lines
        sum_im = repmat(sum(rot_im, 2), 1, size(rot_im,2));
    end

    %Invert the rotation
    inv_im = imrotate(sum_im, theta, 'bilinear', 'crop');
    
    %Add the central part of inverted rotation to the current line sum
    line_sum = line_sum + inv_im(pad_size+1:pad_size+row, pad_size+1:pad_size+col);
    
    %Now check if we've got a complete set of angles for angle subband
    if ~rem(ii, band_size)
        %Save the current line sum in the angle sub-band
        angle_bands(:,:,band) = line_sum;
        
        %figure; imagesc(angle_bands(:,:,band)); axis image;
        
        %Reset the line sum to zeros
        line_sum = zeros(row, col);
        
        %Increment the band counter
        band = band + 1;
        
    end
end
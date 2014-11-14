function [f_i1 f_i2 n_i n_plus N_i p_i  K_i] =...
    karssemeijer_radial_projection(line_map, orientation_map, r_min, r_max, R, num_angles, min_thresh)
%KARSSEMEIJER_RADIAL_PROJECTION *Insert a one line summary here*
%   [] = radial_line_projection()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes: Based on the spicule measures presented in the paper "Detection of
% Stellate Distortions in Mammograms", Karssemeijer & Te Brake, IEE TMI,
% 1996
% 
% This implements a fixed version, see karssemeijer_radial_projection2 for
% an adaptable scale version
% See also:
%
% Created: 26-May-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 
if nargin < 7
    min_thresh = 1;
end
if nargin < 6 || rem(num_angles(1), num_angles(2))
    num_angles = [48 24];
end
if nargin < 5
    R = 5;
end
if nargin < 4
    r_max = 90;
end
if nargin < 3
    r_min = 10;
end


%Make filters:
max_height = ceil(r_max*sin(pi / (2*num_angles(2))));
xx = repmat(-r_max:r_max, 2*max_height+1, 1);
yy = repmat((-max_height:max_height)', 1, 2*r_max+1);
dd = sqrt(xx.^2 + yy.^2);

%Row and column filters to compute n_i - sum of all pixels of a given angle
%directed at a given pixel with tolerance radius R
f_ni_h = [ones(1,r_max-r_min) zeros(1,2*r_min+1) ones(1,r_max-r_min)];
f_ni_v = ones(2*R+1, 1);

%Sum of all pixels within a particular theta and distance range used to
%compute N_i
f_Ni = (abs(atan(yy ./ xx)) < (pi / (2*num_angles(2)))) & dd < r_max & dd > r_min;

%Average distance of pixels within the range used to compute N_i
f_dist = 2*R ./ (pi * dd);
f_dist(~f_Ni) = 0;

%Get size of regions
[row col] = size(line_map);

%Workout padding needed to allow rotation
pad_size = ceil(0.5*(sqrt(2)-1)*max([row col]));

%Pre-allocation outputs
N_i = zeros(row, col);
n_i = zeros(row, col);
p_i = zeros(row, col);
n_plus = zeros(row, col);
K_i = zeros(row, col);

%Compute angular resolution
ang_res = 180 / num_angles(1);

%Compute number of angles in each band and initialise band counter
band_size = num_angles(1) / num_angles(2);

n_ik = zeros(row, col);

%For each angle
for ii = 1:num_angles(1)

    %Compute theta
    theta = (ii - 0.5)*ang_res;
    
    %Get mask of pixels that have orientation within theta range
    theta_mask = (orientation_map > theta - 0.5*ang_res) & (orientation_map <= theta+ 0.5*ang_res);
    
    %Mask out all but selected pixels and pad image for rotation
    pad_im = padarray(line_map .* theta_mask, [pad_size pad_size], 0);
    
    %Rotate image by -theta, so selected pixels project horizontally
    rot_im = imrotate(pad_im, -theta, 'bilinear', 'crop');
    
    %Convolve rotated image with f to compute Ni 
    n_ik_rot = conv2(f_ni_v, f_ni_h, double(rot_im), 'same'); clear rot_im;
    
    %Invert the rotation and add the central portion to n_ik
    inv_im = imrotate(n_ik_rot, theta, 'bilinear', 'crop'); clear nik_rot;
    n_ik = n_ik + inv_im(pad_size+1:pad_size+row, pad_size+1:pad_size+col);
    
    if ~rem(ii, band_size) 
        %Add n_ik to n_i and reset n_ik to zero
        n_i = n_i + n_ik;
        
        %Compute N_ik and R_ik using a rotated copy of unmasked line_map
        pad_im = padarray(line_map, [pad_size pad_size], 0);
    
        %Rotate image by -theta, so selected pixels project horizontally
        rot_im = imrotate(pad_im, -theta, 'bilinear', 'crop'); clear pad_im;
    
        %Convolve rotated image with f_Ni to compute Nik 
        N_ik_rot = conv2(double(rot_im), double(f_Ni), 'same');
    
        %Invert the rotation and take the central portion
        inv_im = imrotate(N_ik_rot, theta, 'bilinear', 'crop'); clear Nik_rot;
        N_ik = inv_im(pad_size+1:pad_size+row, pad_size+1:pad_size+col);
    
        %Add N_ik to N_i
        N_i = N_i + N_ik;
    
        %Convolve rotated image f_dist to compute pik
        p_ik_rot = conv2(double(rot_im), f_dist, 'same'); clear rot_im;
    
        %Invert the rotation and take the central portion
        inv_im = imrotate(p_ik_rot, theta, 'bilinear', 'crop'); clear pik_rot;
        p_ik = inv_im(pad_size+1:pad_size+row, pad_size+1:pad_size+col);
    
        %Add to p_i (later p_i will be divided by Ni)
        p_i = p_i + p_ik;
    
        %Divide pik by Nik
        %p_ik = p_ik ./ N_ik;
        %p_ik(N_ik < 1) = 0; %deal with errors where Nik is unreliable
    
        %Now use Nik and pik to compute media values m_k
        %mik = interp2(N_looukp, p_lookup, m_lookup, Nik, pik, 'linear');
        %mik = Nik .* pik;
        m_k = floor(p_ik);
    
        %See if nik is greater than mik (using non-integer increments - see
        %paper referenced in notes)
        n_plus_k = min(max(n_ik - m_k, 0), 1);
    
        %display(['Number of points to worry about : ' num2str(sum(n_ik(:) < (m_k(:)+1) & n_ik(:) > m_k(:)))...
        %    ' out of ' num2str(length(n_ik(:)))]);
        
        %However, check we've enough weight in Nik to use this bin
        usable_idx = N_ik > min_thresh;
    
        %Update Ki map and ni_plus map
        n_plus(usable_idx) = n_plus(usable_idx) + n_plus_k(usable_idx);
        K_i(usable_idx) = K_i(usable_idx) + 1;
        
        %Finally, reset n_ik to zero
        n_ik = zeros(row, col);
    end
end

%Now we can divide p_i by N_i
p_i = p_i ./ N_i;
usable_idx = N_i >= 1;
p_i(~usable_idx) = 0; %deal with errors where N_i unreliable

%Finally compute the two features from n_i, N_i, p_i, n_plus and K_i
f_i1 = (n_i - p_i.*N_i) ./ sqrt(N_i.*p_i.*(1-p_i));
f_i1(~usable_idx) = min(f_i1(usable_idx));

f_i2 = (n_plus - K_i/2) ./ sqrt(K_i / 4);
f_i2(K_i < 1) = min(f_i2(K_i >= 1));
function [f_i1 f_i2 n_i n_plus N_i p_i  K_i] =...
    karssemeijer_radial_projection_pix(line_map, orientation_map, r_min, r_max, R, num_angles, min_thresh, spacing)
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
if nargin < 8,	spacing = 1; end
if nargin < 7,	min_thresh = 1; end
if nargin < 6,	num_angles = 24; end
if nargin < 5,	R = 5; end
if nargin < 4,	r_max = 90; end
if nargin < 3,	r_min = 10; end

%Make distance and angle array needed to compute n_i, N_i etc:
[xx,yy] = meshgrid(-r_max:r_max,r_max:-1:-r_max);

%theta is angle of any point in square region to centre pixel - angles are
%4 sector defined over [0, 2*pi)
theta = mod(atan2(yy, xx), 2*pi);

%rij is distance of any point in square region to centre pixel
rij = sqrt(xx.^2 + yy.^2);

%Phi is the maximum angular difference allowed between a the orientation
%assigned to a point and it's a associated theta, such that the point is
%directed within a circle of radius R about the centre pixel
phi = abs(asin(R./rij));

%pp, gives the probability of any point being directed to a circle of
%radius R assuming the point has a random orientation
pp = 2*R ./ (pi * rij); 
pp(r_max+1,r_max+1) = 0; %Central pixel is never used in final features but inf will cause nan's later

%Get size of regions
[row col] = size(line_map);

%Zero pad the maps so we can use pixels near the edge of the image - the
%zero padding is taken into in the normalisation factors of the features
%computed
line_map = padarray(line_map, [r_max r_max], 0);
orientation_map = padarray(orientation_map, [r_max r_max], 0);

%Pre-allocation outputs
N_i = zeros(row, col);
n_i = zeros(row, col);
p_i = zeros(row, col);
n_plus = zeros(row, col);
K_i = zeros(row, col);

%Compute angular resolution
ang_res = 2*pi / num_angles(1);

%Make mask for each angle now and store
Ni_mask = cell(1,num_angles);
for ii = 1:num_angles
    %Make filters:
    theta_min = ang_res*(ii-1);
    theta_max = ang_res*ii;

    %Sum of all pixels within a particular theta and distance range used to
    %compute N_i
    Ni_mask{ii} = logical(	(rij > r_min) & (rij < r_max) & ...
							(theta >= theta_min) & (theta < theta_max) );
end

% preset some indices rather than use find() later
inds0 = (1:numel(pp))';

%Run through each pixel in the image, computing N_i, n_i, p_i, K_i and
%n_plus as we go
for xi = 1:spacing:col
    for yi = 1:spacing:row
		%Extract region of line and orientation map - remember padding
		%on image
		line_roi = line_map(yi:yi+2*r_max,xi:xi+2*r_max);
		ori_roi = orientation_map(yi:yi+2*r_max,xi:xi+2*r_max);
		
		% precompute indices of line pixels
		line_roi_inds = inds0(line_roi);
		
		%For each angle
        for ii = 1:num_angles
			% get indices of line pixels for this bin
			% this is equivalent to (but faster than):
			%	inds = find(line_roi & Ni_mask(:,:,ii));
			inds = line_roi_inds(Ni_mask{ii}(line_roi_inds));

			if isempty(inds), continue;	end
			
            %Get all line points that lie within this sector
            N_ik = length(inds);
			
            %Get all line points that lie within the sector and are
            %directed towards the centre
			theta_diff = abs(theta(inds)-ori_roi(inds));
            n_ik_roi =	(theta_diff < phi(inds)) | ...
                (abs(theta_diff-pi) < phi(inds));
            n_ik = sum(n_ik_roi(:));
            
            %Compute sum of probabilities of line points lying within the
            %sector - note we won't divide by N_i until we've finished as
            p_k_N_ik = sum( pp(inds) );

			%Approximate median of n_ik if pixels were randomly oriented as
            %the floor of this sum
            m_k = floor(p_k_N_ik);
            
			%Increment sums across all angle sectors
            N_i(yi, xi) = N_i(yi, xi) + N_ik;
            n_i(yi, xi) = n_i(yi, xi) + n_ik;
			p_i(yi, xi) = p_i(yi, xi) + p_k_N_ik; % normalize later

			%If there are enough line points in this sector, update
            %n_plus and K
            if N_ik > min_thresh;
                n_plus(yi, xi) = n_plus(yi, xi) + min(max(n_ik - m_k, 0), 1);
                K_i(yi, xi) = K_i(yi, xi) + 1;
            end	
        end
    end
end

% Now we can divide p_i by N_i, dealing with errors where N_i unreliable
p_i = p_i ./ N_i;
usable_idx = (N_i >= 1);
p_i(~usable_idx) = 0; 

%Finally compute the two features from n_i, N_i, p_i, n_plus and K_i
f_i1 = (n_i - p_i.*N_i) ./ sqrt(N_i.*p_i.*(1-p_i));
f_i1(~usable_idx) = min(f_i1(usable_idx));

f_i2 = (n_plus - K_i/2) ./ sqrt(K_i / 4);
f_i2(K_i < 1) = min(f_i2(K_i >= 1));

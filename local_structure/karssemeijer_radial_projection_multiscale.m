function [f_i1 f_i2 mask] =...
    karssemeijer_radial_projection_multiscale(line_map, orientation_map, varargin)
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
% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'r_min', 10,...
    'r_max', 90,...
    'R', 5,...
    'num_angles', 24,...
    'min_thresh', 10,...
    'spacing', 1,...
    'offset_x', 0,...
    'offset_y', 0,...
    'mask', [],...
    'ignore_map',[]);
clear varargin;

r_min = args.r_min;
r_max = args.r_max;
R = args.R;
num_angles = args.num_angles;
min_thresh = args.min_thresh;
spacing = args.spacing;
offset_x = args.offset_x;
offset_y = args.offset_y;
mask = args.mask;
ignore_map = args.ignore_map;
clear args;

%Get number of scales
num_scales = length(r_max);

%Make distance and angle array needed to compute n_i, N_i etc:
[xx,yy] = meshgrid(-r_max(end):r_max(end),r_max(end):-1:-r_max(end));

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
pp(r_max(end)+1,r_max(end)+1) = 0; %Central pixel is never used in final features but inf will cause nan's later

%Get size of regions
[row col] = size(line_map);

%Make mask of points at which to compute features
if isempty(mask)
    mask = true(row, col);
end
spacing_mask = false(row,col);
spacing_mask(offset_y+(1:spacing:end),offset_x+(1:spacing:end)) = true;
mask = spacing_mask & mask; clear spacing_mask;
[y_pts x_pts] = find(mask);
num_pts = length(y_pts);

%Check if we have a map specifying sectors to ignore
if isempty(ignore_map)
    ignore_map = inf(row, col);
end

%Zero pad the maps so we can use pixels near the edge of the image - the
%zero padding is taken into in the normalisation factors of the features
%computed
line_map = padarray(line_map, [r_max(end) r_max(end)], 0);
orientation_map = mod(orientation_map, pi);
orientation_map = padarray(orientation_map, [r_max(end) r_max(end)], 0);

%Pre-allocation outputs
f_i1 = zeros(num_pts, num_scales);
f_i2 = zeros(num_pts, num_scales);
min_f1 = 0;
min_f2 = 0;

%Compute angular resolution
ang_res = 2*pi / num_angles(1);

%Make mask for each angle now and store
N_ik_mask = cell(num_angles, num_scales);
phi_k = cell(num_angles);
theta_k = cell(num_angles);
pp_k = cell(num_angles);
for kk = 1:num_angles
    %Make filters:
    theta_min = ang_res*(kk-1);
    theta_max = ang_res*kk;

    %Sum of all pixels within a particular theta and distance range used to
    %compute N_i
    N_ik_mask{kk,num_scales} = (rij(:) > r_min) & (rij(:) < r_max(num_scales)) & ...
        (theta(:) >= theta_min) & (theta(:) < theta_max) ;
    
    %Use this mask to pre-compute the thetas and phis and angle probabilities
    %used in each sector
    phi_k{kk} = phi(N_ik_mask{kk,num_scales});
    %theta_k{kk} = theta(N_ik_mask{kk,num_scales});
    theta_k{kk} = mod(theta(N_ik_mask{kk,num_scales}),pi);
    pp_k{kk} = pp(N_ik_mask{kk,num_scales});
    
    %Now compute the subset of points used for each smaller r_max relative
    %to the largest r_max
    for jj = 1:num_scales-1
        scale_idx = (rij(:) > r_min) & (rij(:) < r_max(jj)) & ...
            (theta(:) >= theta_min) & (theta(:) < theta_max) ; 
        N_ik_mask{kk,jj} = scale_idx( N_ik_mask{kk,num_scales} );
    end
end
            
%Run through each pixel in the image, computing N_i, n_i, p_i, K_i and
%n_plus as we go
for ii = 1:length(x_pts)

    xi = x_pts(ii);
    yi = y_pts(ii);
    
    %Extract region of line and orientation map - remember padding
    %on image
    line_roi = line_map(yi:yi+2*r_max(end),xi:xi+2*r_max(end));
    ori_roi = orientation_map(yi:yi+2*r_max(end),xi:xi+2*r_max(end));

    %Pre-allocate intermdeiate features
    N_i = zeros(1, num_scales);
    n_i = zeros(1, num_scales);
    p_i = zeros(1, num_scales);
    n_plus = zeros(1, num_scales);
    K_i = zeros(1, num_scales);

    %For each angle
    for kk = 1:num_angles
        
        if ignore_map(yi, xi) == kk || (ignore_map(yi, xi) + num_angles/2) == kk
            continue;
        end

        % get indices of line pixels for this bin
        % this is equivalent to (but faster than):
        %	inds = find(line_roi & N_ik_mask(:,:,kk));
        inds = N_ik_mask{kk,num_scales};

        %Get all line points that lie within this sector
        N_ik_roi = line_roi(inds);

        %Compute N_ik now for top level - if this is near zero we can
        %continue to the next sector with out making further
        %computations
        N_ik = sum(N_ik_roi(:));

        if N_ik < 1e-6, continue;	end

        %Get all line points that lie within the sector and are
        %directed towards the centre
        theta_diff = abs(theta_k{kk}-ori_roi(inds));
        theta_inds = (theta_diff < phi_k{kk}) | ...
            (abs(theta_diff-pi) < phi_k{kk});
        n_ik_roi = theta_inds .* N_ik_roi;

        %Compute sum of probabilities of line points lying within the
        %sector - note we won't divide by N_i until we've finished as
        p_ik_roi = pp_k{kk} .* N_ik_roi;

        for jj = num_scales:-1:1
            %Compute N_ik, n_ik, p_k by summing over the necessary pts
            %in the pre-computed regions
            if jj == num_scales
                %N_ik already computed for largest r_max above
                n_ik = sum(n_ik_roi(:));
                p_k_N_ik = sum(p_ik_roi(:));
            else
                scale_inds = N_ik_mask{kk,jj};
                N_ik = sum(N_ik_roi(scale_inds));
                n_ik = sum(n_ik_roi(scale_inds));
                p_k_N_ik = sum(p_ik_roi(scale_inds));
            end

            %Approximate median of n_ik if pixels were randomly oriented as
            %the floor of this sum
            m_k = floor(p_k_N_ik);

            %Increment sums across all angle sectors
            N_i(jj) = N_i(jj) + N_ik;
            n_i(jj) = n_i(jj) + n_ik;
            p_i(jj) = p_i(jj) + p_k_N_ik; % normalize later

            %If there are enough line points in this sector, update
            %n_plus and K
            if N_ik > min_thresh;
                n_plus(jj) = n_plus(jj) + min(max(n_ik - m_k, 0), 1);
                K_i(jj) = K_i(jj) + 1;
            end
        end
    end

    %Check where we have usable data
    usable_f1 = (N_i >= 1);
    usable_f2 = (K_i >= 1);

    % Now we can divide p_i by N_i, dealing with errors where N_i unreliable
    p_i = p_i ./ N_i;

    %Finally compute the two features from n_i, N_i, p_i, n_plus and K_i
    f_i1_scales = (n_i - p_i.*N_i) ./ sqrt(N_i.*p_i.*(1-p_i));
    f_i2_scales = (n_plus - K_i/2) ./ sqrt(K_i / 4);

    %Now deal with places in f1 and f2 where the counts make unreliable
    %features
    min_f1 = min([min_f1 f_i1_scales(usable_f1)]);
    min_f2 = min([min_f2 f_i2_scales(usable_f2)]);
    f_i1_scales(~usable_f1) = min_f1;
    f_i2_scales(~usable_f2) = min_f2;

    %Finally copy feature scores into main pre-allocated storage
    f_i1(ii,:) = f_i1_scales;
    f_i2(ii,:) = f_i2_scales;       
end

        

        

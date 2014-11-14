function [particles_x particles_y] = jetstream(I_mag, I_ori, I_c, init_pt, init_vec, varargin)
%JETSTREAM contour tracking algorithm
%   [contour, particles] = jetstream(varargin)
%
% JETSTREAM uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%      I_mag - Image magnitudes (orieginal paper assumed gradient
%      magnitude, as we're tracing centre lines we prob use soemthing
%      different - e.g. 2nd deriv. response)
%
%      I_ori - Map of orientations (original paper assumed gradient
%      direction here - see note above)
%
%      I_c - Binary map of corners, not implemented fully yet
%
%      init_pt - [x y] coord of start point
%
%      init_vec - [u v] initial step
%
% Optional Arguments:
%      M - (100) number of particles in stream
%
%      N - (100) max number of steps to take (process will also break if no
%      valid particles are found)
%
%      nu - (0.01) expected proportion of corners
%
%      sigma_theta - (0.05) var of Normal dist to sample new directions
%
%      sigma_psi - (1) estimated variance of relationship between image
%      magnitude and psi - the difference between contour direction and
%      image orientation at a given point on the contour
%
%      lambda - parameter of exponential distribution for image magnitude
%      over no-contour points (i.e. used to computed p_off)
%
% Outputs:
%
%      particles - *Insert description of input variable here*
%
%
% Example:
%
% Notes: see "JetStream: probabilistic contour extraction with particles", Computer Vision, 2001. ICCV 2001 
%
% See also: JETSTREAM_RF, JETSTREAM2, JETSTREAM3
%
% Created: 25-Oct-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'd', 1,...
    'nu', 0.01,...
    'M', 100,...
    'sigma_theta', 0.05,...
    'sigma_psi', 1,...
    'lambda', [],...
    'N', 100, ...
    'plot', 0);
clear varargin;

[rows cols] = size(I_mag);
if isempty(I_c);
    I_c = zeros(rows, cols);
end
if isempty(args.lambda)
    %Compute lambda from I
     args.lambda = mean(abs(I_mag(:)));
end

%set up initial particles and container for future particles
particles_x = zeros(args.M, args.N + 2);
particles_y = zeros(args.M, args.N + 2);
particles_x(:,1) = ones(args.M,1)*init_pt(1);
particles_y(:,1) = ones(args.M,1)*init_pt(2);
particles_x(:,2) = particles_x(:,1) + init_vec(1);
particles_y(:,2) = particles_y(:,1) + init_vec(2);

if args.plot
    figure; imagesc(I_mag); axis image; colormap(gray(256)); hold all;
end

%Now loop iterate particles trhough each step
for ii = 2:args.N+1
    
    %sample new thetas - NEED TO ADD SAMPLING BASED ON CORNER DETECTION
    theta_i = args.sigma_theta * randn(args.M, 1);
    
    %predict new particle positions
    xi = [particles_x(:,ii) particles_y(:,ii)];
    xim1 = [particles_x(:,ii-1) particles_y(:,ii-1)];
    for mm = 1:args.M
        R_theta = [cos(theta_i(mm)) sin(theta_i(mm)); -sin(theta_i(mm)) cos(theta_i(mm))];
        xip1 = xi(mm,:) + (xi(mm,:) - xim1(mm,:))*R_theta;
        particles_x(mm,ii+1) = xip1(1);
        particles_y(mm,ii+1) = xip1(2);
    end
    
    %compute weightings for new particle positions
    N_theta_i = pdf('norm', theta_i, 0, args.d * args.sigma_theta);
    q_theta_i = args.nu / pi + (1 - args.nu)*N_theta_i;
    
    psi_ip1 = mb_mod(...
        interp2(I_ori, particles_x(:,ii+1), particles_y(:,ii+1), 'bilinear') - ...
        atan((xim1(:,2)-xi(:,2))./ (xi(:,1)-xim1(:,1))), pi); %Negate  the y-component because of ij indexing
    I_mag_ip1 = interp2(I_mag, particles_x(:,ii+1), particles_y(:,ii+1), 'bilinear');
    p_on = pdf('norm', psi_ip1, 0, args.sigma_psi ./ sqrt(abs(I_mag_ip1)));
    p_off = exp(-I_mag_ip1 / args.lambda);
    p_off(p_off > 1) = 1;
    
    l_ip1 = p_on ./ p_off;
    
    pi_ip1 = q_theta_i .* l_ip1 .* N_theta_i;
    pi_ip1(isnan(pi_ip1)) = 0;
    
    if ~any(pi_ip1)
        break;
    end
    
    %Now sample from the predicted points using the weights as sampling
    %probabilities
    probs_ip1 = [0; cumsum(pi_ip1)] / sum(pi_ip1);
    am = zeros(args.M, 1);
    for mm = 1:args.M
        am(mm) = find(rand >= probs_ip1, 1, 'last');
    end
    particles_x = particles_x(am,:);
    particles_y = particles_y(am,:);
    
    if args.plot
        plot(particles_x(:,ii+1), particles_y(:,ii+1), '.');
    end
end
if args.plot
    plot(mean(particles_x), mean(particles_y), 'k');
end
    

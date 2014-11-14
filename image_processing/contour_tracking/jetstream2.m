function [particles_x particles_y] = jetstream2(I_mag, I_ori, I_c, init_pt, init_vec, init_m, varargin)
%JETSTREAM2 contour tracking algorithm - extension to track edges of linear
%structure rather than a single contour (i.e. 2 parallel streams)
%   [contour, particles] = jetstream(varargin)
%
% JETSTREAM2 uses the U_PACKARGS interface function
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
% See also: JETSTREAM, JETSTREAM_RF, JETSTREAM3
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
    'sigma_m', 0,...
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
particles_xm = zeros(args.M, args.N + 2);
particles_ym = zeros(args.M, args.N + 2);
particles_xp = zeros(args.M, args.N + 2);
particles_yp = zeros(args.M, args.N + 2);
particles_x(:,1) = ones(args.M,1)*init_pt(1);
particles_y(:,1) = ones(args.M,1)*init_pt(2);
particles_x(:,2) = particles_x(:,1) + init_vec(1);
particles_y(:,2) = particles_y(:,1) + init_vec(2);

m_i = zeros(args.M, args.N + 2);
m_i(:,1:2) = ones(args.M,2)*init_m;

if args.plot
    figure; imagesc(I_mag); axis image; colormap(gray(256)); hold all;
end

%Now loop iterate particles trhough each step
for ii = 2:args.N+1
    
    %sample new thetas - NEED TO ADD SAMPLING BASED ON CORNER DETECTION
    theta_i = args.sigma_theta * randn(args.M, 1);
    
    %sample new thetas
    m_i(:,ii+1) = m_i(:,ii) + args.sigma_m * randn(args.M, 1);
    
    %predict new particle positions
    xi = [particles_x(:,ii) particles_y(:,ii)];
    xim1 = [particles_x(:,ii-1) particles_y(:,ii-1)];
    for mm = 1:args.M
        R_theta = [cos(theta_i(mm)) sin(theta_i(mm)); -sin(theta_i(mm)) cos(theta_i(mm))];
        xip1 = xi(mm,:) + (xi(mm,:) - xim1(mm,:))*R_theta;
        particles_x(mm,ii+1) = xip1(1);
        particles_y(mm,ii+1) = xip1(2);
        
        particles_xm(mm,ii+1) = xip1(1) + m_i(mm)*(xi(mm,2) - xim1(mm,2));
        particles_ym(mm,ii+1) = xip1(2) - m_i(mm)*(xi(mm,1) - xim1(mm,1));
        
        particles_xp(mm,ii+1) = xip1(1) - m_i(mm)*(xi(mm,2) - xim1(mm,2));
        particles_yp(mm,ii+1) = xip1(2) + m_i(mm)*(xi(mm,1) - xim1(mm,1));
    end
    
    %compute weightings for new particle positions
    N_theta_i = pdf('norm', theta_i, 0, args.d * args.sigma_theta);
    q_theta_i = args.nu / pi + (1 - args.nu)*N_theta_i;
    
    psi_ip1p = mb_mod(...
        interp2(I_ori, particles_xp(:,ii+1), particles_yp(:,ii+1), 'bilinear') - ...
        atan((xim1(:,2)-xi(:,2))./ (xi(:,1)-xim1(:,1))), pi); %Negate  the y-component because of ij indexing
    
    psi_ip1m = mb_mod(...
        interp2(I_ori, particles_xm(:,ii+1), particles_ym(:,ii+1), 'bilinear') - ...
        atan((xim1(:,2)-xi(:,2))./ (xi(:,1)-xim1(:,1))), pi); %Negate  the y-component because of ij indexing
    
    
    I_mag_ip1 = interp2(I_mag, particles_x(:,ii+1), particles_y(:,ii+1), 'bilinear');
    I_mag_ip1p = interp2(I_mag, particles_xp(:,ii+1), particles_yp(:,ii+1), 'bilinear');
    I_mag_ip1m = interp2(I_mag, particles_xm(:,ii+1), particles_ym(:,ii+1), 'bilinear');
    
    p_onp = pdf('norm', psi_ip1p, 0, args.sigma_psi ./ sqrt(abs(I_mag_ip1p)));
    p_onm = pdf('norm', psi_ip1m, 0, args.sigma_psi ./ sqrt(abs(I_mag_ip1m)));
    
    p_off = exp(-I_mag_ip1 / args.lambda);
    p_off(p_off > 1) = 1;
    
    p_offp = exp(-I_mag_ip1p / args.lambda);
    p_offp(p_offp > 1) = 1;
    
    p_offm = exp(-I_mag_ip1m / args.lambda);
    p_offm(p_offm > 1) = 1;
    
    l_ip1 = (p_off .* p_onp .* p_onm) ./ (p_offp .* p_offm);
    
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
    particles_xp = particles_xp(am,:);
    particles_yp = particles_yp(am,:);
    particles_xm = particles_xm(am,:);
    particles_ym = particles_ym(am,:);
    m_i = m_i(am,:);
    
    if args.plot
        plot(particles_x(:,ii+1), particles_y(:,ii+1), '.');
        %plot(particles_xm(:,ii+1), particles_ym(:,ii+1), '.');
    end
end
if args.plot
    plot(mean(particles_xp), mean(particles_yp), 'r');
    plot(mean(particles_xm), mean(particles_ym), 'g');
end
    

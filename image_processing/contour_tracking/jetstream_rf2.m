function [particles_x particles_y] = jetstream_rf2(I_mag, I_ori, I_jun, init_pt, init_vec, varargin)%, roc_vals
%JETSTREAM *Insert a one line summary here*
%   [contour, particles] = jetstream(varargin)
%
% JETSTREAM uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%      contour - *Insert description of input variable here*
%
%      particles - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
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

if isempty(args.lambda)
    %Compute lambda from I
     args.lambda = mean(abs(I_mag(:)));
end
[rows cols] = size(I_mag);
if isempty(I_jun);
    I_jun = zeros(rows, cols);
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
    
    %Interpolate junction status
    I_jun_i = I_jun(sub2ind([rows cols], round(particles_y(:,ii)), round(particles_x(:,ii))));
    
    %sample new thetas - either uniform at junction or normal if not
    theta_i = I_jun_i.*(rand(args.M, 1)*pi - pi/2) + ... %junctions
        (1-I_jun_i).* args.sigma_theta .* randn(args.M, 1); %non-junctions
    
    %predict new particle positions
    xi = [particles_x(:,ii) particles_y(:,ii)];
    xim1 = [particles_x(:,ii-1) particles_y(:,ii-1)];
    for mm = 1:args.M
        R_theta = [cos(theta_i(mm)) sin(theta_i(mm)); -sin(theta_i(mm)) cos(theta_i(mm))];
        xip1 = xi(mm,:) + (xi(mm,:) - xim1(mm,:))*R_theta;
        particles_x(mm,ii+1) = xip1(1);
        particles_y(mm,ii+1) = xip1(2);
    end
    
    %Convert particle xy coords into indices
    part_idx = sub2ind([rows cols], round(particles_y(:,ii+1)), round(particles_x(:,ii+1)));
    
    %compute weightings for new particle positions
    N_theta_i = pdf('norm', theta_i, 0, args.d * args.sigma_theta);
    q_theta_i = args.nu / pi + (1 - args.nu)*N_theta_i;
    
    
    %Interpolate magnitudes at potential particles
    I_mag_i = I_mag(part_idx);
    
    %Interpolate orientations at potential particles
    I_ori_i = I_ori(part_idx);
    
    %IDEA: interpolate ROC curve to compute p_on and p_off given the vessel
    %prob map?
    %p_on_off = interp(args.roc_pts, roc_vals, I_mag_ip1);
    %l_ip1 = p_on_off(:,1) ./ p_on_off(:,2);
    
    %Current solution (HACK?) - use vessel map as positive weighting and
    %prob of angle psi as positive weighting
    
    %Compute probabilities of psi given the sampled image orientation
    % - double the angle of the current vector
    psi_i = 2*atan((xim1(:,2)-xi(:,2))./ (xi(:,1)-xim1(:,1)));
    
    %Standardise (with wrapping)
    psi_iz = mb_mod(psi_i - angle(I_ori_i), 2*pi) ./ abs(I_ori_i);
    
    %Compute probability of psi: if if junction (I_jun_i==1), then p_psi = 1/pi,
    %otherwise p_psi_z ~ N(0,1)
    p_psi = I_jun_i/pi + (1-I_jun_i).*pdf('norm', psi_iz, 0, 1);
    
    %Compute final probability of a particle
    p_i = q_theta_i .* p_psi .* I_mag_i ./ (I_jun_i/pi + (1-I_jun_i).* N_theta_i);
    p_i(isnan(p_i)) = 0;
    
    if ~any(p_i)
        break;
    end
    
    %Now sample from the predicted points using the weights as sampling
    %probabilities
    probs_ip1 = [0; cumsum(p_i)] / sum(p_i);
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
    

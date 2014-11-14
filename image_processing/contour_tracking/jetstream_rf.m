function [particles_x particles_y] = jetstream_rf(I_mag, I_ori, I_jun, init_pt, init_vec, varargin)
%JETSTREAM_RF contour tracking algorithm that hacks jetsream for use with our
%RF orientation and detection outputs
%   [contour, particles] = jetstream(varargin)
%
% JETSTREAM_RF uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%      I_mag - Image magnitudes (e.g. vessel probability map from RF classifier)
%
%      I_ori - Map of orientations in complex form (e.g. from RF regressor)
%
%      I_jun - Binary map of junction locations
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
%      double_angle - (1) flag set to 1 if orientation is in double angle
%      form (and therefore needs to be halved)
%
%      step_length
%
% Outputs:
%
%      particles - *Insert description of input variable here*
%
%
% Example:
%
% Notes: see "JetStream: probabilistic contour extraction with particles", Computer Vision, 2001. ICCV 2001. 
%
% See also: JETSTREAM, JETSTREAM2, JETSTREAM3
%
% Created: 25-Oct-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'd', 1,...
    'M', 100,...
    'step_length', 1,...
    'N', 100, ...
    'double_angle', 1,...
    'plot', 0);
clear varargin;

[rows cols] = size(I_mag);

if isempty(I_jun);
    I_jun = zeros(rows, cols);
end
if args.double_angle
    dbl_ang = 2;
else
    dbl_ang = 1;
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
    
    idx_i = sub2ind([rows cols], round(particles_y(:,ii)), round(particles_x(:,ii)));
    
    %Interpolate junction status
    I_jun_i = I_jun(idx_i);
    
    %Interpolate orientations at potential particles
    I_ori_i = I_ori(idx_i);
    
    %predict new particle positions
    xi = [particles_x(:,ii) particles_y(:,ii)];
    xim1 = [particles_x(:,ii-1) particles_y(:,ii-1)];
    for mm = 1:args.M
        %if at junction sample direction
        if I_jun_i
            R_theta = 2*pi*rand;
        else
        %otherwise sample from distrbution defined by direction and
        %dispersion at current position
            theta = angle(I_ori_i(mm))/dbl_ang;
            rho = abs(I_ori_i(mm));
            
            if (cos(theta)*(xi(mm,1) - xim1(mm,1)) - sin(theta)*(xi(mm,2) - xim1(mm,2))) < 0
                theta = theta - pi;
            end
            sigma = sqrt(-2*log(rho));
            R_theta = theta + sigma * randn / 2;
            %R_theta = theta + wrapped_normal_sample(0, rho, 1)/2;
        end
        
        %Convert direction to x,y steps and check we're not going back on
        %ourselves
        u = cos(R_theta);
        v = -sin(R_theta);
%         if (u*(xi(mm,1) - xim1(mm,1)) + v*(xi(mm,2) - xim1(mm,2))) < -0.5
%             u = -u;
%             v = -v;
%         end

        %Add steps to generate new particles
        particles_x(mm,ii+1) = xi(mm,1) + args.step_length*u;
        particles_y(mm,ii+1) = xi(mm,2) + args.step_length*v;
    end
    
    %Convert particle xy coords into indices
    idx_ip1 = sub2ind([rows cols], round(particles_y(:,ii+1)), round(particles_x(:,ii+1)));    
    
    %Interpolate magnitudes at potential particles
    I_mag_ip1 = I_mag(idx_ip1);
    
    %Compute final probability of a particle
    p_i = I_mag_ip1;
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
    plot(mean(particles_x), mean(particles_y), 'g', 'linewidth', 2);
end
    

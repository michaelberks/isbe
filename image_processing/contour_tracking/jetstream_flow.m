function [particles_x particles_y particles_p particles_d hit_sink] =...
    jetstream_flow(I_mag, I_flow, I_dist, I_sink, init_pt, varargin)
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
%      I_flow - Map of flow in complex form (e.g. from Phil's flow optimisation algorithm)
%
%      I_jun - Binary map of junction locations (now in options)
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
%      normal_max_step_length
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
    'normal_max_step_length', 4,...
    'junction_max_step_length', 25,...
    'junction_angle_spread', 2*pi/3,...
    'N', 100, ...
    'double_angle', 0,...
    'use_angle_weights', 0,...
    'I_jun', [],...
    'plot', 0);
clear varargin;

[rows cols] = size(I_mag);

if ~isempty(args.I_jun);
    check_junctions = true;
    I_jun = args.I_jun;
else
    check_junctions = false;
end
if args.double_angle
    dbl_ang = 2;
else
    dbl_ang = 1;
end

%set up initial particles and container for future particles
particles_x = zeros(args.M, args.N + 2);
particles_y = zeros(args.M, args.N + 2);
particles_p = zeros(args.M, args.N + 2);
particles_d = zeros(args.M, args.N + 2);

particles_x(:,1) = ones(args.M,1)*init_pt(1);
particles_y(:,1) = ones(args.M,1)*init_pt(2);
particles_p(:,1) = 1 / args.M;
particles_d(:,1) = interp2(I_dist, particles_x(1,1), particles_y(1,1), '*nearest');

if args.plot
    figure; imagesc(I_mag); axis image; colormap(gray(256)); hold all;
    plot(particles_x(:,1), particles_y(:,1), '.');
end

%Now loop iterate particles trhough each step
max_flow = max(abs(I_flow(:)));

last_step = args.N+1;
for ii = 1:last_step
    
    idx_i = sub2ind([rows cols], round(particles_y(:,ii)), round(particles_x(:,ii)));
    hit_sink = I_sink(idx_i);
    
    if all(hit_sink) %sum(I_sink(idx_i)) > args.M/2
        last_step = ii;
        break;
    end
    
    %Interpolate junction status
    if check_junctions
        I_jun_i = I_jun(idx_i);
    end
    
    %Interpolate orientations at potential particles
    I_ori_i = I_flow(idx_i) / max_flow;
    
    %predict new particle positions
    xi = [particles_x(:,ii) particles_y(:,ii)];
    
    if args.use_angle_weights && ii > 1
        xim1 = [particles_x(:,ii-1) particles_y(:,ii-1)];       
    end
    
    angle_weights = ones(args.M,1);
    dist_weights = zeros(args.M,1);
    vessel_weights = zeros(args.M,1);
    for mm = 1:args.M
        
        %If particle in the sink, don't bother streaming new values
        if hit_sink(mm)
            continue;
        end
        
        %if at junction sample direction from anywhere in half circle
        %centre on previous direction (i.e. anywhere 'forwards')
        if check_junctions && I_jun_i(mm) && ii > 1
            prev_theta = atan2(...
                particles_y(mm,ii) - particles_y(mm,ii-1),...
                particles_x(mm,ii) - particles_x(mm,ii-1));
            R_theta = prev_theta + args.junction_angle_spread*(rand-0.5);
        else
        %otherwise sample from distrbution defined by direction and
        %dispersion at current position
            theta = angle(I_ori_i(mm))/dbl_ang;
            rho = abs(I_ori_i(mm));
            
            sigma = sqrt(-2*log(rho));
            R_theta = theta + sigma * randn / 2;
            %R_theta = theta + wrapped_normal_sample(0, rho, 1)/2;
        end
        
        %Convert direction to x,y steps and check we're not going back on
        %ourselves
        u = cos(R_theta);
        v = sin(R_theta);
        
        if args.use_angle_weights && ii > 1
            u0 = xi(mm,1) - xim1(mm,1);
            v0 = xi(mm,2) - xim1(mm,2);
            m0 = sqrt(u0^2 + v0^2);

            angle_weights(mm) = sqrt(((u+u0/m0)/2)^2 + ((v+v0/m0)/2)^2);
        end


        %If at junction, look further afield to try and cross to the
        %otherside
        if check_junctions && I_jun_i(mm)
            step_length = args.junction_max_step_length*rand;
        else
            step_length = args.normal_max_step_length*rand;
        end
        particles_x(mm,ii+1) = min(max(xi(mm,1) + step_length*u,1),cols);
        particles_y(mm,ii+1) = min(max(xi(mm,2) + step_length*v,1),rows);

        % Quick approximation to interp2()
        dx = abs((1:cols) - particles_x(mm,ii+1));
        xinds = (dx < 1);
        dy = abs((1:rows)' - particles_y(mm,ii+1));
        yinds = (dy < 1);
        
        if any(xinds) || any(yinds)
            xwts = 1 - dx(xinds);
            ywts = 1 - dy(yinds);
            
            %Iterpolated vessel prob
            vessel_weights(mm) = (ywts' * I_mag(yinds,xinds) * xwts');

            %Interpolate distance
            particles_d(mm,ii+1) = (ywts' * I_dist(yinds,xinds) * xwts'); 
        end
        dist = particles_d(mm,ii+1) - particles_d(mm,ii);
        dist_weights(mm) = exp(-abs(dist - step_length)/step_length);
        
    end
    
    %Normalise weights
    angle_weights = angle_weights / sum(angle_weights(~hit_sink));
    dist_weights = dist_weights / sum(dist_weights(~hit_sink));
    vessel_weights = vessel_weights / sum(vessel_weights(~hit_sink));
    
    %Compute final probability of a particle
    p_i = vessel_weights .* dist_weights .*angle_weights;
    p_i(isnan(p_i)) = 0; %particles in the sink will return zero because their vessel and dist weights are intialised to zero and not changed
    
    if ~any(p_i)
        last_step = ii;
        break;
    end
    
    %For particles in the sink, just copy over their data
    particles_x(hit_sink,ii+1) = particles_x(hit_sink,ii);
    particles_y(hit_sink,ii+1) = particles_y(hit_sink,ii);
    particles_p(hit_sink,ii+1) = particles_p(hit_sink,ii);
    particles_d(hit_sink,ii+1) = particles_d(hit_sink,ii);
    
    %Now sample from the predicted points using the weights as sampling
    %probabilities
    particles_p(~hit_sink,ii+1) = p_i(~hit_sink) / sum(p_i);
    probs_ip1 = [0; cumsum(particles_p(~hit_sink,ii+1))];
    
    free_particles = sum(~hit_sink);
    free_particles_idx = find(~hit_sink);
    am = zeros(free_particles, 1);
    for mm = 1:free_particles
        am(mm) = free_particles_idx(find(rand >= probs_ip1, 1, 'last'));
    end
    particles_x(~hit_sink,ii:ii+1) = particles_x(am,ii:ii+1);
    particles_y(~hit_sink,ii:ii+1) = particles_y(am,ii:ii+1);
    particles_p(~hit_sink,ii:ii+1) = particles_p(am,ii:ii+1);
    particles_d(~hit_sink,ii:ii+1) = particles_d(am,ii:ii+1);   
    
    if args.plot
        plot(particles_x(:,ii+1), particles_y(:,ii+1), '.');
    end
end

particles_x(:,last_step+1:end) = [];
particles_y(:,last_step+1:end) = [];
particles_p(:,last_step+1:end) = [];
particles_d(:,last_step+1:end) = [];    

if args.plot
    [~,max_idx] = max(particles_p(:,end));
    plot(mean(particles_x), mean(particles_y), 'g', 'linewidth', 2);
    plot(particles_x(max_idx,:), particles_y(max_idx,:), 'r', 'linewidth', 2);
end
    

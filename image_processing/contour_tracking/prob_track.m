function [path_map num_steps] = prob_track(I_prob, I_ori_theta, I_ori_D, x1, y1, step_length, junction_map)
%PROB_TRACK *Insert a one line summary here*
%   [path_map] = prob_track(I_mag, I_ori_theta, I_ori_D, xs, xy)
%
% Inputs:
%      I_prob - Map of structure (e.g. vessel) probabilites, most likely
%       returned from a classifier
%
%      I_ori_theta - Map of structures orientations - direction (i.e.
%      angle) component (this should NOT be doubled)
%
%      I_ori_D - Map of structures orientations - dispersion (i.e magnitude) component
%
%      x1 - x coord of start point
%
%      y1 - y coord of start point
%
%      step_length - does what it says on the tin
%
%      junction_map - (optional) binary map indicating whether there is
%       structure junction (i.e. bifurcation or crossing) present
%
%
% Outputs:
%      path_map - Map traced by path
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 23-Feb-2012
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if nargin < 7 || isempty(junction_map)
    %assume no junctions present if no argument given
    junction_map = false(size(I_prob));
end

%preallocate path map
path_map = zeros(size(I_prob));

%Set up intitial points and driection
xi = x1;
yi = y1;
t1 = I_ori_theta(round(yi), round(xi));
rd = 1 - 2*(rand > .5);
xs0 = rd*cos(t1);
ys0 = -rd*sin(t1);

%Uncomment if we want to return the xy coords of the path 
%x = xi;
%y = yi;
go_on = true;
step = 1;
while go_on

    %Determine next direction to travel
    if junction_map(round(yi), round(xi)) > rand
        %if at junction, sample from uniform distribution
        new_theta = 2*pi*rand;
    else
        %otherwise sample from orientation distribution (defined by
        %direction and dispersion) at current location
        theta = I_ori_theta(round(yi), round(xi));
        rho = I_ori_D(round(yi), round(xi));
        new_theta = theta + wrapped_normal_sample(0, rho, 1)/2;
    end

    %Convert direction into x,y steps
    xs1 = cos(new_theta);
    ys1 = -sin(new_theta);

    %Check we're not going back on ourselves by more than 45 degrees
    if xs1*xs0 + ys1*ys0 < -0.5
        xs1 = -xs1;
        ys1 = -ys1;
    end

    %Step in the chosen direction
    xi = xi + step_length*xs1;
    yi = yi + step_length*ys1;
    xs0 = xs1;
    ys0 = ys1;
    
    %Save the new point into the path map
    path_map(round(yi), round(xi)) = path_map(round(yi), round(xi)) + 1;
        
    %Uncomment if we want to return the xy coords of the path 
    %x = [x; xi];
    %y = [y; yi];
    
    %Check whether we're going to continue from the new point
    p = I_prob(round(yi), round(xi));
    go_on = (p > 0.5) || (p > rand); %This is obviously a heuristic hack at the moment
    
    %Increment the step count
    step = step + 1;
end
num_steps = step;
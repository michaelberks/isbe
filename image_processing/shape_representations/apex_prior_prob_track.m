function [apex_prior_score path_map] = ...
    apex_prior_prob_track(I_prob, I_ori_theta, I_ori_sigma, x1, y1, step_length, num_streams, ...
        max_steps, lambda, apex_prior)
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

%preallocate path map
patch_sz = size(I_prob);
path_map = zeros(patch_sz);
t1 = I_ori_theta(y1, x1);
apex_prior_score = zeros(max_steps,2);


[prior_length, n_ori, ~] = size(apex_prior);
max_steps = min(max_steps, prior_length);

%get orientation at initisl point
switch_count = 0;
for i_stream = 1:num_streams
    %Set up intitial points and driection
    xi = x1;
    yi = y1;
    xii = round(xi);
    yii = round(yi);
    
    %Pick an initial random direction (split 50 between either orientation
    %hemisphere)
    rd = 1 - 2*(rand > .5);
    xs0 = rd*cos(t1);
    ys0 = -rd*sin(t1);
    
    if xs0 < 0
        prior_idx = 1;
    else
        prior_idx = 2;
    end

    go_on = true;
    step = 0;

    while go_on
        %Increment the step count
        step = step + 1;
        
        if step == max_steps;
            go_on = false;
            continue;
        end
            
        %sample from orientation distribution (defined by direction and dispersion) at current location
        theta = I_ori_theta(yii, xii);
        sigma = I_ori_sigma(yii, xii);
        
        %Switch theta to be in half circle we already travelling in
        if cos(theta)*xs0 - sin(theta)*ys0 < 0
            theta = theta - pi;
            switch_count = switch_count + 1;
        end
        
        new_theta = theta + sigma * randn / 2;

        %Convert direction into x,y steps
        xs1 = cos(new_theta);
        ys1 = -sin(new_theta);
         
        %Work orientation relative to t1
        ori_change = mb_mod(2*(new_theta-t1),2*pi);
        ori_bin = ceil(.5*n_ori*(ori_change + pi)/pi);
        prior_weight = apex_prior(step, ori_bin, prior_idx);

        %Step in the chosen direction
        xi = xi + step_length*xs1;
        yi = yi + step_length*ys1;
        xs0 = xs1;
        ys0 = ys1;
        
        %Check path hasn't gone off the edge of the image
        xii = round(xi);
        yii = round(yi);
        if (xii > 0) && (xi <= patch_sz(2)) && (yii > 0) && (yi <= patch_sz(1))            

            %Save the new point into the path map
            path_map(yii, xii) = path_map(yii, xii) + 1;%prior_weight;
            apex_prior_score(step, prior_idx) = apex_prior_score(step, prior_idx) +...
                prior_weight;

            %Check whether we're going to continue from the new point
            p_im = I_prob(yii, xii);
            p_rand = rand;
            go_on = (p_rand > lambda) || (p_im > p_rand); %This is obviously a heuristic hack at the moment
        else
            go_on = false;
        end        
    end
end
%display(['Total swicthes = ' num2str(switch_count)]);
            
            
    
    






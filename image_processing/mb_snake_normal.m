function [snake_pnts,e] = mb_snake_normal(pnts, alpha, beta, max_delta_x, resol_x, ...
                          normal_p, normal_x, normal_y)
                      
%MB_SNAKE_NORMAL Mike Berks adaption of snake algorithm, where image energy
% comes from normal profiles sampled from an image, whilst shape energy
% comes from the position of the points used to form the profiles in the
% original image (see MB's thesis for more details)
%
% Inputs:
%   pnts          Starting contour. Each row is a [x,y] coordinate.
%   alpha         Energy contributed by the distance between control points.
%                 Set to zero if length of slope doesn't matter.
%   beta          Energy contributed by the curvature of the snake.  Larger
%                 values of beta cause bends in the snake to have a high cost
%                 and lead to smoother snakes.
%   max_delta_y   Max number of pixels to move each contour point vertically
%   resol_y       Contour points will be moved by multiples of resol_y
%   max_delta_x   Analog to max_delta_y
%   resol_x       Analog to resol_y
%   normal_p      2D-Array of the feature responses in the image.  For example
%                 it can contain the magnitude of the image gradients
%   normal_x/y    x/y-coordinates of points in the original image (shape
%                 energy is derived from these points, NOT the coordinates
%                 of the points in the profile image
%
% Outputs:
%   snake_pnts    New contour points.
%   e             Energy value of these new contour points
%
%
% Example:  [snake_pnts,e] = mb_snake_normal(pnts, alpha, beta, ...
%                        max_delta_y, resol_y, max_delta_x, resol_x, normal_p, normal_x, normal_y)
%
% Notes: Adapted from SNAKE to which the copyright probably belongs, see
% MBs thesis for full details
%
% See also: SNAKE
%
% Created: 10-Jun-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5125 
% Copyright: (C) University of Manchester
% SNAKES
% By Chris Bregler and Malcolm Slaney, Interval Technical Report IRC 1995-017
% Copyright (c) 1995 Interval Research Corporation.

  %Minimum resolution to move points in x/y is 1 pixel  
  if resol_x < 1; resol_x = 1; end;

  %Get number of points
  n = size(pnts,1);
  
  %Get size of feature image
  [row,col] = size(normal_p);
  
  %Reshape feature_image to a single column vector
  target = reshape(normal_p,row*col,1);
  
  %Generate vectors for the possible places to search in the x/y directions
  %neighbouring a point
  delta_x = -max_delta_x:resol_x:max_delta_x;
  num_states = size(delta_x,2);

  states_x = round(pnts(:,1))*ones(1,num_states) + ones(n,1)*delta_x;
  states_y = round(pnts(:,2))*ones(1,num_states);

  % take care of boundary cases
  states_x = min(max(states_x,1),col);
  states_y = min(max(states_y,1),row);

  %states_i is the indices of all possible points we will investigate in
  %this iteration
  states_i = (states_x-1)*row + states_y;
  
  %Now set states_x and states_y based on normal_x and normal_y
  pos_x = normal_x(states_i);
  pos_y = normal_y(states_i);
  
  Smat = zeros(n,num_states^2);
  Imat = zeros(n,num_states^2);

  %-----------------------------------------------------------------------
  % forward pass
  
  %Copy the values of the feature image into Smat for all possible states
  %of the first point
  for v2 = 1:num_states,
    Smat(1,(v2-1)*num_states+1:v2*num_states) = -target(states_i(1,:))';
  end
    
  %For the 2nd to last but one points compute the energy sum associated with
  %possible states of each point
  for k = 2:n-1,

    for v2 = 1:num_states, 
        for v1 = 1:num_states,
            
            %What is does v0_domain do?
            v0_domain = 1:num_states;
            
            %Find the minimum of a sum?             
            continuity = ...
                + alpha * ( ...
                    (states_x(k,v1)-states_x(k-1,v0_domain)).^2 ...
                +   (states_y(k,v1)-states_y(k-1,v0_domain)).^2) ...
                ...
                + beta * (...
                    (pos_x(k+1,v2)-2*pos_x(k,v1) + pos_x(k-1,v0_domain)).^2 ...
                +   (pos_y(k+1,v2)-2*pos_y(k,v1) + pos_y(k-1,v0_domain)).^2);
            
            energy_sum = ...
                Smat(k-1,(v1-1)*num_states+v0_domain) + continuity;
            
            %Compute the minimum energy, and the state that achieved this
            [y,i] = min(energy_sum);
            
            Smat(k,(v2-1)*num_states+v1) = y-target(states_i(k,v1));
            Imat(k,(v2-1)*num_states+v1) = i;
            
        end
    end
    %display(['Continuity = ', num2str(min(continuity)), ', Energy = ', num2str(target(states_i(k,v1)))]);
  end

  %For the last point compute the energy sum associated with all possible
  %states
  for v1 = 1:num_states,
    v0_domain = 1:num_states;
    [y,i] = min( Smat(n-1,(v1-1)*num_states+v0_domain) ...
        + alpha * ( ...
            (states_x(n,v1)-states_x(n-1,v0_domain)).^2 + ...
            (states_y(n,v1)-states_y(n-1,v0_domain)).^2));
    Smat(n,v1) = y-target(states_i(n,v1));
    Imat(n,v1) = i;
  end;

  %Compute the minimal energy of the final point and the state that
  %achieved it
  [e, final_i] = min(Smat(n,1:num_states));


  %------------------------------------------------------------------------
  % backward pass

  snake_pnts = zeros(n,2);

  %Fix the final point to the one that achieved the minimum on the foward
  %pass
  snake_pnts(n,:) = [states_x(n,final_i), states_y(n,final_i)];
  
  v1 = final_i; 
  v2 = 1;
  
  %Working backwards through the points...
  for k = n-1:-1:1
      
    v = Imat(k+1,(v2-1)*num_states+v1);
    v2 = v1; 
    v1 = v;
    snake_pnts(k,:) = [states_x(k,v1),states_y(k,v1)];
  end;
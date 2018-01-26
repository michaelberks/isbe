function [opt_ellipse_xy, opt_ellipse_params] = optimise_ellipse(feature_im, init_ellipse, opt_options)
%OPTIMISE_ELLIPSE *Insert a one line summary here*
%   [opt_ellipse_xy, opt_ellipse_params] = optimise_ellipse(feat_im, init_ellipse, opt_params)
%
% Inputs:
%      feature_im - fitted ellipse maximise score interpolated at ellipse
%      points
%
%      init_ellipse - either 1x5 vector of ellipse parameters [rx, ry, x0,
%      y0, theta] or [nx2] vector of xy coordinates defining an outline to
%      which an ellipse will be fitted
%
%      opt_params - optional parameters for the optimisation
%
%
% Outputs:
%      opt_ellipse_xy - xy co-ordinates of fitted ellipse
%
%      opt_ellipse_params - parameters of fitted ellipse
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 03-Aug-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if size(init_ellipse,2) == 2
    %Assume init ellipse is an outline of points
    exy = init_ellipse;
    init_ellipse = zeros(1,5);
    [init_ellipse(1), init_ellipse(2), init_ellipse(3), init_ellipse(4), x_axis] = ...
        fit_ellipse(exy);
    
    init_ellipse(5) = atan2(x_axis(2), x_axis(1));
end

if ~exist('opt_options', 'var') || isempty(opt_options)
    % %Set options for optimisation
    opt_options = optimoptions('fmincon',...
        'Display', 'final',...
        'MaxIter', 20);
end

%Estimate the perimeter of the initial ellipse
a = init_ellipse(1);
b = init_ellipse(2);
n_pts = ceil(pi*( 3*(a+b) - sqrt((3*a+b)+(a+3*b)) ));

%Estimate some upper and lower bounds
lower_bounds = [0.9*a 0.9*b init_ellipse(3:4)-1 init_ellipse(5)-pi/6];
upper_bounds = [1.1*a 1.1*b init_ellipse(3:4)+1 init_ellipse(5)+pi/6];

%objective function, optimising model params in signal space 
obj_fun = @(x)fit_ellipse_to_im(x,...
    feature_im, n_pts); %fixed variables in the objective function
        
opt_ellipse_params = fmincon(obj_fun, init_ellipse, [], [], [], [],...
    lower_bounds, upper_bounds, [], opt_options);

r_x = opt_ellipse_params(1);
r_y = opt_ellipse_params(2);
x0 = opt_ellipse_params(3);
y0 = opt_ellipse_params(4);
x_axis = [cos(opt_ellipse_params(5)) sin(opt_ellipse_params(5))];
[ex ey] = ellipse(r_x, r_y, x0 , y0, x_axis, 0, n_pts);
opt_ellipse_xy = [ex(1:end-1)' ey(1:end-1)'];

function [score] = fit_ellipse_to_im(ellipse_params, feature_im, n_pts)

r_x = ellipse_params(1);
r_y = ellipse_params(2);
x0 = ellipse_params(3);
y0 = ellipse_params(4);
x_axis = [cos(ellipse_params(5)) sin(ellipse_params(5))];

[ex ey] = ellipse(r_x, r_y, x0 , y0, x_axis, 0, n_pts);
ellipse_feat = interp2(feature_im, ex(1:end-1), ey(1:end-1), 'bilinear');
score = -sum(ellipse_feat);
%plot(ex, ey);


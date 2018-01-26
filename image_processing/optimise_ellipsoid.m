function [opt_center, opt_radii, opt_axes, opt_e3_params, opt_e3_vol] =...
    optimise_ellipsoid(feature_im, im_x, im_y, im_z, init_ellipsoid, varargin)
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
args = u_packargs(varargin, 0, ...
    'fixed_idx', 13:15,...
    'lower_bounds', [],...
    'upper_bounds', [],...
    'max_centre_offsets', [2 2 1],...
    'max_radii_scaling', [0.1 0.1 0.1],...
    'weights', [],...
    'num_itr', 200,...
    'opt_options', [],...
    'debug', false);
clear varargin;

if size(init_ellipsoid,2) == 2
    %Assume init ellipse is an outline of points
    exy = init_ellipsoid;
    [center, radii, e_axes] = ellipsoid_fit( exy );
    init_ellipsoid = [e_axes(:)' center(:)' radii(:)'];
elseif length(init_ellipsoid) == 10
    [center, radii, e_axes] = ellipsoid_params(init_ellipsoid);
    init_ellipsoid = [e_axes(:)' center(:)' radii(:)'];
end

if isempty(args.fixed_idx)
    opt_idx = 1:12;
else
    opt_idx = setdiff(1:15, args.fixed_idx);
end
if isempty(args.opt_options)
    % %Set options for optimisation
    args.opt_options = optimoptions('fmincon',...
        'Display', 'final',...
        'MaxIter', args.num_itr);
end

%Estimate some upper and lower bounds
%Axes should not have upper or lower bounds, they are contrained by
%their format
if isempty(args.lower_bounds)
    
    args.lower_bounds = [-Inf(1,9) init_ellipsoid(10:12)-args.max_centre_offsets];
    if length(args.fixed_idx) == 15
        args.lower_bounds = ...
            [args.lower_bounds init_ellipsoid(13:15).*(1-args.max_radii_scaling)];
    end
end
if isempty(args.upper_bounds)
    
    args.upper_bounds = [Inf(1,9) init_ellipsoid(10:12)+args.max_centre_offsets];
    if length(args.fixed_idx) == 15
        args.upper_bounds = ...
            [args.upper_bounds init_ellipsoid(13:15).*(1+args.max_radii_scaling)];
    end
end

if isempty(args.weights)
    e_axes = reshape(init_ellipsoid(1:9),3,3);
    e_center = init_ellipsoid(10:12);
    e_radii = init_ellipsoid(13:15);
    v = ellipsoid_general(e_center, e_radii, e_axes);
    [exyz] = compute_ellipsoid_volume(v, im_x, im_y, im_z);
    args.weights = 1 / mean(exyz(exyz > 0));
end

%objective function, optimising model params in signal space 
obj_fun = @(x)fit_ellipsoid_to_im(x,...
    feature_im, im_x, im_y, im_z,...
    init_ellipsoid(args.fixed_idx), args.fixed_idx, opt_idx, args.weights); %fixed variables in the objective function

init_axes = reshape(init_ellipsoid(1:9),3,3);
con_fun = @(x)constrain_axes(x, init_axes);
        
A = [];
b = [];
Aeq = [];
beq = [];

ep = fmincon(obj_fun, init_ellipsoid(opt_idx), A, b, Aeq, beq,...
    args.lower_bounds, args.upper_bounds, con_fun, args.opt_options);

op = init_ellipsoid;
op(opt_idx) = ep;
opt_axes = reshape(op(1:9),3,3);
opt_center = op(10:12);
opt_radii = op(13:15);

[opt_e3_params] = ellipsoid_general(opt_center, opt_radii, opt_axes);

opt_e3_vol = compute_ellipsoid_volume(opt_e3_params, im_x, im_y, im_z);

function [score] = fit_ellipsoid_to_im(ellipse_params, feature_im, im_x, im_y, im_z, fixed_params, fixed_idx, opt_idx, weights)

ep = zeros(15,1);
ep(opt_idx) = ellipse_params;
ep(fixed_idx) = fixed_params;
e_axes = reshape(ep(1:9),3,3);
e_center = ep(10:12);
e_radii = ep(13:15);
ev = ellipsoid_general(e_center, e_radii, e_axes);
[wxyz] = weight_ellipsoid_volume(ev, im_x, im_y, im_z, weights);

score = -sum(feature_im(:).*wxyz(:));
% figure; 
% for i_slice = 1:11; 
%     subplot(3,4,i_slice); 
%     imgray(feature_im(:,:,i_slice).*wxyz(:,:,i_slice)); 
% end

function [c, ceq] = constrain_axes(ellipse_params, init_axes)

e_axes = reshape(ellipse_params(1:9),3,3);
id_mat = e_axes*e_axes' - eye(3);
rot_mat = init_axes*e_axes' - eye(3);

c = [sum(id_mat(:).^2) - 1e-6; sum(rot_mat(:).^2) - pi/12];
ceq = det(e_axes) - 1;


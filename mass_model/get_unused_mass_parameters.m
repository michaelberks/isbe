function [unused_model] = get_unused_mass_parameters(mass_model_name, unused_list, mass_dir)
%GET_UNUSED_MASS_PARAMETERS *Insert a one line summary here*
%   [unused_model] = get_unused_mass_parameters(mass_model,unused_list)
%
% Inputs:
%      mass_model- *Insert description of input variable here*
%
%      unused_list- *Insert description of input variable here*
%
%
% Outputs:
%      unused_model- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also: GENERATE_MASS_AM
%
% Created: 10-Dec-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

if nargin < 3
    mass_dir = 'C:\isbe\dev\masses\';
end
load(mass_model_name);

N = length(unused_list);
target_shape = reshape(mass_model.mean_shape,[],2);

mean_shape  = mass_model.mean_shape;
mean_scale  = mass_model.mean_scale;
mean_tex    = mass_model.mean_tex;

P_shape     = mass_model.P_shape;
P_scale     = mass_model.P_scale;
P_tex       = mass_model.P_tex;

P_com       = mass_model.P_com;

W_shape     = mass_model.W_shape;
W_tex       = mass_model.W_tex;
W_scale     = mass_model.W_scale;

k_shape     = size(P_shape, 2);
k_tex       = size(P_tex, 2);

size_shape_vec = length(mean_shape) / 2;

s_x = mass_model.scaled_mean_shape(1:size_shape_vec);% + mean_centre(1);
s_y = mass_model.scaled_mean_shape(size_shape_vec+1:end);% + mean_centre(2);

%Define points to be interpolated by TPS - as row vectors
i_x = mass_model.mean_shape_pl(:,1)';
i_y = mass_model.mean_shape_pl(:,2)';

X_shape = zeros(N, 2*size_shape_vec);
X_tex = zeros(N, length(i_x));
X_scale = zeros(N,1);
X_com = zeros(N, k_shape+k_tex+1);

rotations = zeros(2,2,N);
translations = zeros(N,2);
origins = zeros(N,1);

B_shape = zeros(k_shape, N);
B_tex = zeros(k_tex, N);
B_scale = zeros(1, N);
B_com = zeros(size(P_com, 2), N);

[shapes_unaligned] =...
    get_shapes_from_masses(unused_list, size_shape_vec, 'mass_path', mass_dir);

for ii = 1:N
    source_shape = reshape(shapes_unaligned(ii,:), [], 2);
    min_dist = inf;
    for oo = 1:500
        shift_shape = circshift(source_shape, oo - 1);
        [dd Z t] = mb_procrustes(target_shape, shift_shape);
        if dd < min_dist
            min_dist = dd;
            shape = Z;
            transform = t;
            origin = oo;
        end
    end
    X_shape(ii, :) = shape(:)';
    X_scale(ii) = transform.b;
    rotations(:,:,ii) = transform.T;
    translations(ii,:) = transform.c(1,:) / transform.b;
    origins(ii) = origin;
    
    
    mass = u_load([mass_dir, unused_list(ii).name]);
    z_x = source_shape(:,1)';
    z_y = source_shape(:,2)'; 
    
    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
        'transform', 'spline');
    [pts] = geom_transformpoints([i_x; i_y], T);

    clear T;
    
    %Get texture vector
    X_tex(ii,:) = interp2(mass.subtract_ROI, pts(1,:), pts(2,:));
    X_tex(ii,isnan(X_tex(ii,:))) = 0;
    
    b_shape = P_shape' * (X_shape(ii,:) - mean_shape)';
    b_tex = P_tex' * (X_tex(ii,:) - mean_tex)';
    b_scale = P_scale' * (X_scale(ii) - mean_scale)';
    x_c = [W_shape*b_shape; W_tex*b_tex; W_scale*b_scale]';
    b_c = P_com' * x_c';
    
    X_com(ii,:)     = x_c;
    B_shape(:,ii)   = b_shape;
    B_tex(:,ii)     = b_tex;
    B_scale(:,ii)   = b_scale;
    B_com(:,ii)     = b_c;
end
X_tex(isnan(X_tex)) = 0;

unused_model.X_shape    = X_shape;
unused_model.X_tex      = X_tex;
unused_model.X_scale    = X_scale;
unused_model.X_com      = X_com;

unused_model.B_shape    = B_shape;
unused_model.B_tex      = B_tex;
unused_model.B_scale    = B_scale;
unused_model.B_com      = B_com;

unused_model.mass_model_name = mass_model_name;
unused_model.unused_list = unused_list;

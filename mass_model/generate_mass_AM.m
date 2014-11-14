% GENERATE_MASS_AM implementation of Caulkin's method for modelling mass
%           appearance 
%    [mass_model model_id] = generate_mass_AM(mass_files, file_out,...
%         size_shape_vec, size_tex_vec, varargin)
%
%    inputs:
%       
%
%    outputs:
%       
%
%    notes:
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [mass_model model_id]...
    = generate_mass_AM(mass_files, file_out, varargin)

args = u_packargs(varargin, 0,...
    'mass_path', 'C:\isbe\dev\masses\',...
    'size_shape_vec', 500,...
    'size_tex_vec', 50000,...
    'spline', 'isbe',...
    'shiftOrigin', 1);
mass_path = args.mass_path;
size_shape_vec = args.size_shape_vec;
size_tex_vec = args.size_tex_vec;

% Create ID structure to save model info
model_id.mass_files = mass_files;
model_id.args = args;

%
% Create matrix of shape vectors X_shape from data
%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('1: Method started');
mass_model.progress = '1: Method started';
N = length(mass_files);

[shapes_unaligned] =...
    get_shapes_from_masses(mass_files, size_shape_vec, 'mass_path', mass_path);

% Align the shapes - X_shape is the matrix of aligned shapes, X_scale the
% corresponding scale parameters
[X_shape, X_scale, mean_target mass_model.rotations, ...
    mass_model.translations, mass_model.origins] =...
    align_shapes(shapes_unaligned,...
    'area', size_tex_vec, 'shiftOrigin', args.shiftOrigin);
%
% Compute shape model from shape vectors
%%%%%%%%%%%%%%%%%%%%%%%%
[mean_shape, P_shape, B_shape, L_shape] = pca(X_shape, 0.98);
[mean_scale, P_scale, B_scale, L_scale] = pca(X_scale, 0.98);

mass_model.shapes_unaligned = shapes_unaligned;
mass_model.X_shape = X_shape;
mass_model.X_scale = X_scale;

mass_model.mean_shape = mean_shape;
mass_model.P_shape = P_shape;
mass_model.B_shape = B_shape;
mass_model.L_shape = L_shape;
mass_model.mean_scale = mean_scale;
mass_model.P_scale = P_scale;
mass_model.B_scale = B_scale;
mass_model.L_scale = L_scale;
mass_model.mean_target = mean_target;

display('2: Shape model complete');
mass_model.progress = '2: Shape model complete';
save(file_out, 'mass_model', 'model_id');

display('3: Background subtraction complete');
mass_model.progress = '3: Background subtraction complete';
save(file_out, 'mass_model', 'model_id');

%
% Compute texture model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create mean shape pixel list
%%%%%%%%%%

if isempty(size_tex_vec)
    scale_factor = 1;
else
    scale_factor = sqrt(size_tex_vec / ...
        polyarea(mean_shape(1:size_shape_vec),...
        mean_shape(size_shape_vec + 1: end)));
end

scaled_mean_shape = mean_shape*scale_factor;

% Set mean centre to leave sufficient room at edges
mean_centre(1) = 101 - floor(min(scaled_mean_shape(1:size_shape_vec))); % x
mean_centre(2) = 101 - floor(min(scaled_mean_shape(size_shape_vec+1:end))); %y

% Move scaled_mean_shape to mean centre
scaled_mean_shape(1:size_shape_vec) = scaled_mean_shape(1:size_shape_vec)...
    + mean_centre(1);
scaled_mean_shape(size_shape_vec+1:end) = scaled_mean_shape(size_shape_vec+1:end)...
    + mean_centre(2);
mass_model.scaled_mean_shape = scaled_mean_shape;

% Create a rectangular array big enough for the mean shape to fit in
mean_row = ceil(max(scaled_mean_shape(size_shape_vec+1:end))) + 100; % y
mean_col = ceil(max(scaled_mean_shape(1:size_shape_vec))) + 100; % x

mass_model.mean_centre = mean_centre;
mass_model.mean_row = mean_row;
mass_model.mean_col = mean_col;

mean_bw = poly2mask(scaled_mean_shape(1:size_shape_vec),...
            scaled_mean_shape(size_shape_vec+1:end),...
            mean_row, mean_col);

%for ii = 1:49; mean_bw = imdilate(mean_bw, strel('disk', 1)); end

[mean_shape_pl(:,2) mean_shape_pl(:,1)] = find(mean_bw);
mass_model.mean_shape_pl = mean_shape_pl;

display('4: Mean pixel list created');
mass_model.progress = '4: Mean pixel list created';
save(file_out, 'mass_model', 'model_id');

% Warp each shape to mean shape and extract texture vector
%%%%%%%%%%

%Define source points for TPS - as row vectors
s_x = scaled_mean_shape(1:size_shape_vec);% + mean_centre(1);
s_y = scaled_mean_shape(size_shape_vec+1:end);% + mean_centre(2);

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';

if strcmp(args.spline, 'mab')
    tps_L_inv = tps_weights(s_x, s_y, 'biharmTPS');
end

X_tex = zeros(N, length(i_x));
for ii=1:N
    
    display(['warping mean to shape ', int2str(ii)]);
    
    temp = load([mass_path, mass_files(ii).name]);
    mass = temp.mass; clear temp;
    
    %Define displacement to target points
    z_x = shapes_unaligned(ii, 1:size_shape_vec);% + mass.mass_centroid(1);
    z_y = shapes_unaligned(ii, size_shape_vec+1:end);% + mass.mass_centroid(2);  
    
    switch args.spline
        case 'mab'
            %Compute displacement of interpolated points
            pts = tps_warp(tps_L_inv, s_x, s_y, z_x, z_y, i_x, i_y, 'biharmTPS');
            
        case 'isbe'
            T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                'transform', 'spline');
            [pts] = geom_transformpoints([i_x; i_y], T);
            
            clear T;
    end       
    
    %Get texture vector
    X_tex(ii,:) = interp2(mass.subtract_ROI, pts(1,:), pts(2,:));

%     helpful visual check on spline warps    
%     if ~rem(ii, 20)
%         temp = zeros(mean_row, mean_col);
%         temp(sub2ind([mean_row, mean_col], i_y, i_x)) = uint8(X_tex(ii,:));
%         figure('WindowStyle', 'docked'); hold on;
%         image(temp); axis image; colormap(gray(256));
%         figure; hold on;
%         image(mass.subtract_ROI); axis image; colormap(gray(256));
%         plot(pts(1,1:10:end), pts(2,1:10:end), 'r.');
%         plot(mass.mass_outline(:,1), mass.mass_outline(:,2));
%         plot(z_x, z_y, 'g.');
%         clear temp;
%     end
    clear mass pts;
    
end
clear s_x s_y z_x z_y i_x i_y tps_L_inv;

X_tex(isnan(X_tex)) = 0;
mass_model.X_tex = X_tex;

display('5: Shapes warped to mean shape');
mass_model.progress = '5: Shapes warped to mean shape';
save(file_out, 'mass_model', 'model_id');
%
%Compute model of texture vectors
%%%%%%%%%%%%%%%%%%%%%
[mean_tex, P_tex, B_tex, L_tex] = pca(X_tex, 0.98);%, 0);

mass_model.mean_tex = mean_tex;
mass_model.P_tex = P_tex;
mass_model.B_tex = B_tex;
mass_model.L_tex = L_tex;

display('6: Texture model complete');
mass_model.progress = '6: Texture model complete';
save(file_out, 'mass_model', 'model_id');

%
%Calculate weights for combined model
%%%%%%%%%%%%%%%%%%%%%%
k_shape = length(L_shape);
k_tex   = length(L_tex);

W_shape = k_shape / sum(sqrt(L_shape));
W_tex   = k_tex / sum(sqrt(L_tex));
W_scale = 1 / sqrt(L_scale); %length L_scale = 1

mass_model.W_shape = W_shape;
mass_model.W_tex = W_tex;
mass_model.W_scale = W_scale;

combined_data = [W_shape*B_shape; W_tex*B_tex; W_scale*B_scale]';
mass_model.combined_data = combined_data;

[mean_com, P_com, B_com, L_com] = pca(combined_data, 0.98);%, 0);

mass_model.mean_com = mean_com;
mass_model.P_com = P_com;
mass_model.B_com = B_com;
mass_model.L_com = L_com;

display('7: Combined model complete, function successful!');
mass_model.progress = '7: Combined model complete, function successful!';
save(file_out, 'mass_model', 'model_id');


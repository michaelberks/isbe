function change_model_shape(mass_model, file_out, mass_files, size_shape_vec, varargin)

default.mass_path = 'C:\isbe\dev\masses\';
args = u_packargs(varargin, 0, default);
mass_path = args.mass_path;

%
% Create matrix of shape vectors X_shape from data
%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('1: Method started');
mass_model.progress = '1: Method started';

[shapes_unaligned, mass_areas] =...
    get_shapes_from_masses(mass_files, size_shape_vec, 'mass_path', mass_path);

% Align the shapes - X_shape is the matrix of aligned shapes, X_scale the
% corresponding scale parameters
[X_shape, X_scale, mean_target] =...
    align_shapes(shapes_unaligned, 'area', mean(mass_areas));
%
% Compute shape model from shape vectors
%%%%%%%%%%%%%%%%%%%%%%%%
[mean_shape, P_shape, B_shape, L_shape] = pca(X_shape, 0.98);
[mean_scale, P_scale, B_scale, L_scale] = pca(X_scale, 0.98);

mass_model.shapes_unaligned = shapes_unaligned;
mass_model.X_shape = X_shape;
mass_model.X_scale = X_shape;

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
save(file_out, 'mass_model');

display('3: Background subtraction complete');
mass_model.progress = '3: Background subtraction complete';
save(file_out, 'mass_model');

%
% Compute texture model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create mean shape pixel list
%%%%%%%%%%


scale_factor = sqrt(length(mass_model.mean_shape_pl) / ...
    polyarea(mean_shape(1:size_shape_vec),...
    mean_shape(size_shape_vec + 1: end)));

mass_model.scale_factor = scale_factor;

% Set mean centre to leave sufficient room at edges
mass_model.mean_centre(1) = mass_model.mean_off_c;
mass_model.mean_centre(2) = mass_model.mean_off_r;
mass_model = rmfield(mass_model, 'mean_off_c');
mass_model = rmfield(mass_model, 'mean_off_r');

%
%Calculate weights for combined model
%%%%%%%%%%%%%%%%%%%%%%
k_shape = length(L_shape);

W_shape = k_shape / sum(sqrt(L_shape));
W_scale = 1 / sqrt(L_scale); %length L_scale = 1

mass_model.W_shape = W_shape;
mass_model.W_scale = W_scale;

combined_data = [W_shape*B_shape; mass_model.W_tex*mass_model.B_tex; W_scale*B_scale]';
mass_model.combined_data = combined_data;

[mean_com, P_com, B_com, L_com] = pca(combined_data, 0.98);%, 0);

mass_model.mean_com = mean_com;
mass_model.P_com = P_com;
mass_model.B_com = B_com;
mass_model.L_com = L_com;

mass_model = rmfield(mass_model, 'mean_c');
mass_model = rmfield(mass_model, 'P_c');
mass_model = rmfield(mass_model, 'B_c');
mass_model = rmfield(mass_model, 'L_c');

display('7: Combined model complete, function successful!');
mass_model.progress = '7: Combined model complete, function successful!';
save(file_out, 'mass_model');
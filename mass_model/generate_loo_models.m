% GENERATE_LOO_MODELS
%    [com_error ind_error] =... 
%       generate_loo_models(varargin)
%
%    inputs:
%
%       optional:
%
%    outputs:
%
%    notes: 
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks
function generate_loo_models(varargin)

    args = u_packargs(varargin, 0,...
        'method', 'quick',...
        'mass_files', [],...
        'mass_model', [],...
        'mass_path', 'C:\isbe\dev\masses\',...
        'model_path', [],...
        'size_shape_vec', 500,...
        'size_tex_vec', 50000,...
        'spline', 'isbe');
    
    switch args.method
        
        case 'quick'
            if isempty(args.mass_model)
                error('Please supply a model for quick LOO model generation');
            end
            
            N = size(args.mass_model.X_shape, 1);
            for ii = 1:N
                mass_model = args.mass_model;
                mass_model.shapes_unaligned(ii,:) = [];
                mass_model.X_shape(ii,:) = [];
                mass_model.X_scale(ii) = [];
                
                [mean_shape, P_shape, B_shape, L_shape] = pca(mass_model.X_shape, 0.98);
                [mean_scale, P_scale, B_scale, L_scale] = pca(mass_model.X_scale, 0.98);

                mass_model.mean_shape = mean_shape;
                mass_model.P_shape = P_shape;
                mass_model.B_shape = B_shape;
                mass_model.L_shape = L_shape;
                mass_model.mean_scale = mean_scale;
                mass_model.P_scale = P_scale;
                mass_model.B_scale = B_scale;
                mass_model.L_scale = L_scale;
                
                mass_model.X_tex(ii,:) = [];
                [mean_tex, P_tex, B_tex, L_tex] = pca(mass_model.X_tex, 0.98);%, 0);

                mass_model.mean_tex = mean_tex;
                mass_model.P_tex = P_tex;
                mass_model.B_tex = B_tex;
                mass_model.L_tex = L_tex;

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
                model_id.info = 'loo model created by quick method';
                model_name = [args.model_path, 'model', zerostr(ii,3)]; 
                save(model_name, 'mass_model', 'model_id');
            end
        case 'full'
            if isempty(args.mass_files)
                error('Please supply list of masses for full LOO model generation');
            end
            N = length(args.mass_files);
            for ii = 1:N
                loo_files = args.mass_files;
                loo_files(ii) = [];
                model_name = [args.model_path, 'model', zerostr(ii,3)]; 
                generate_mass_AM(loo_files, model_name, args)
            end
    end
     





function [] = compute_apex_candidate_displacements(varargin)
%SELECT_VESSELS_FROM_CANDIDATES *Insert a one line summary here*
%   [] = select_vessels_from_candidates(varargin)
%
% SELECT_VESSELS_FROM_CANDIDATES uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 20-Nov-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, 0, ...
    {'image_names'},...
    'data_dir', [nailfoldroot 'data/rsa_study/master_set/'],...
    'vessel_centre_dir', 'vessel_centres',...
    'displacement_dir', 'displacements',...
    'candidates_dir', 'rescores',...
    'min_candidates', 3,...
    'initial_thresh', 0.3,...
    'grid_spacing', 8, ...
    'use_snake', 0,...
    'alpha', 1,...
    'beta', 1,...
    'snake_delta_y', 30,...
    'snake_res_y', 1,...
    'snake_spacing', 50,...
    'do_plot',  0);            
clear varargin;

%Form full directory paths and create folder for HoGs
vessel_centre_dir = [args.data_dir '/' args.vessel_centre_dir '/'];
displacement_dir = [args.data_dir '/' args.displacement_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];
create_folder(displacement_dir);

num_images = length(args.image_names);

%Loop though each image
for i_im = 1:num_images
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);
    
    im_name = args.image_names{i_im} ;  
    
    load([candidates_dir im_name '_candidates.mat'], 'candidate_xy', 'candidate_rescores');
    load([vessel_centre_dir im_name '_vc.mat'], 'nrows', 'ncols');

    valid_candidates = candidate_rescores > args.initial_thresh; %#ok
    if sum(valid_candidates) < args.min_candidates
        candidate_displacements = inf(size(candidate_rescores)); %#ok
    else
        %Create grid of points spaced over the image
        xx = 1:args.grid_spacing:ncols;
        yy = (1:args.grid_spacing:nrows)';
        grid_x = repmat(xx, length(yy), 1);
        grid_y = repmat(yy, 1, length(xx));
   
        %Compute weighted kernel estimates of the spatial distribution of
        %candidates over this grid
        [location_distribution] = build_2d_kernel_distribution(...
           candidate_xy(valid_candidates,:),...
           [grid_x(:) grid_y(:)],...
           candidate_rescores(valid_candidates,:), 0); %#ok
   
        location_distribution.D_f = reshape(location_distribution.D_f, size(grid_x));
   
        %Compute the max y-coordinate for each x-coordinate of the
        %candidates, and the displacement to this value
        [~, y_max] = max(location_distribution.D_f);
        x_i = xx;
        y_i = yy(y_max);
        
        if args.use_snake

            x_sampled = linspace(1, x_i(end), x_i(end) / (args.snake_spacing));
            y_sampled = interp1(x_i, y_i, x_sampled, 'linear');

            [snake_pnts] = mb_snake(...
                [x_sampled(:) y_sampled(:)] / args.grid_spacing,...
                args.alpha, args.beta, 0, 1, args.snake_delta_y, args.snake_res_y,...
                location_distribution.D_f * 1e6);
            
            x_i = snake_pnts(:,1)*args.grid_spacing;
            y_i = snake_pnts(:,2)*args.grid_spacing;     
           
        end
        
        candidate_polyfit = interp1(x_i, y_i, candidate_xy(:,1), 'linear');
        candidate_displacements = candidate_xy(:,2) - candidate_polyfit;  %#ok 
        
        if args.do_plot
            figure; imgray(location_distribution.D_f);
            plot(x_i/args.grid_spacing, y_i/args.grid_spacing, 'y');
            plot(candidate_xy(:,1)/args.grid_spacing, candidate_xy(:,2)/args.grid_spacing, 'rx');
            plot(...
                [candidate_xy(:,1) candidate_xy(:,1)]'/args.grid_spacing,...
                [candidate_xy(:,2) candidate_polyfit]'/args.grid_spacing, '-');
        end
    end    
    
    save([displacement_dir im_name '_dis.mat'], 'candidate_displacements');

end
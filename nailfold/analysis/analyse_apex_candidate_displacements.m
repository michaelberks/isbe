function [] = analyse_apex_candidate_displacements(varargin)
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
    'candidates_dir', 'rescores',...
    'label_dir', 'labels',...
    'min_candidates', 6,...
    'initial_thresh', 0.3,...
    'grid_spacing', 8, ...
    'poly_n', 5,...
    'do_plot',  0);
clear varargin;

%Form full directory paths and create folder for HoGs
vessel_centre_dir = [args.data_dir '/' args.vessel_centre_dir '/'];
label_dir = [args.data_dir '/' args.label_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];

num_images = length(args.image_names);

%Loop though each image
for i_im = 1:num_images
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);
    
    im_name = args.image_names{i_im} ;  
    
    load([candidates_dir im_name '_candidates.mat'], 'candidate_xy', 'candidate_rescores');
    load([vessel_centre_dir im_name '_vc.mat'], 'nrows', 'ncols');
    load([label_dir im_name '_label.mat'], 'candidates_class', 'candidate_labels');
    
    [candidate_xy sort_idx] = sortrows(candidate_xy);
    candidate_rescores = candidate_rescores(sort_idx);
    candidate_labels = candidate_labels(sort_idx);
    %distal_candidates = candidate_labels == 2;
    
    valid_candidates = candidate_rescores > args.initial_thresh;
    if sum(valid_candidates) >= args.min_candidates %&& sum(distal_candidates) >= args.min_candidates;
        
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
           candidate_rescores(valid_candidates,:), 0);
   
        location_distribution.D_f = reshape(location_distribution.D_f, size(grid_x));
        
        %Compute the max y-coordinate for each x-coordinate of the
        %candidates, and the displacement to this value
        [~, y_max] = max(location_distribution.D_f);
        x_m = xx;
        y_m = yy(y_max);
                                
        max_mask = ...
            location_distribution.D_f(2:end-1,:) > location_distribution.D_f(1:end-2,:) &...
            location_distribution.D_f(2:end-1,:) > location_distribution.D_f(3:end,:);
        
        max_mask = [zeros(1,size(max_mask,2)); double(max_mask); zeros(1,size(max_mask,2))];
        imfilter(max_mask, fspecial('gaussian'));
        
%         x_i = candidate_xy(distal_candidates,1);
%         y_i = candidate_xy(distal_candidates,2);
%         x_i = [xx(1); x_i(:); xx(end)];
%         y_i = [y_i(1); y_i(:); y_i(end)];
        
        %[pp, s, mu] = polyfit(x_i, y_i, args.poly_n); 
        %y_d = polyval(pp, xx, s, mu);
        
%         try
%             [y_d] = ... c
%                 nag_smooth_fit_spline_parest('C', 'U', x_i, y_i, 0, 0, 1e-3, 'u', 1e7);
%         
%         catch err
%             display(err);
%             pp = pchip(x_i, y_i);
%             y_d = ppval(pp, xx);
%             x_i = xx;
%         end
%         

%         
%         try
%             [y_sampled_s] = ... c
%                 nag_smooth_fit_spline_parest('C', 'U', x_sampled, y_sampled, 0, 0, 1e-3, 'u', 1e6);
%         
%         catch err
%             display(err);
%             y_sampled_s = [];
%         end
        
        
        if args.use_snake
            alpha = 1;
            beta = 1;
            max_delta_x = 0;
            resol_x = 1;
            max_delta_y = 30;
            resol_y = 1;
            x_sampled = linspace(1, x_m(end), x_m(end) / (8*args.grid_spacing));
            y_sampled = interp1(x_m, y_m, x_sampled, 'linear');

            [snake_pnts,e] = mb_snake(...
                [x_sampled(:) y_sampled(:)] / args.grid_spacing,...
                alpha, beta, max_delta_x, resol_x, max_delta_y, resol_y,...
                location_distribution.D_f * 1e6);
            
%             [snake_pnts,e] = mb_snake(...
%                 [x_sampled(:) y_sampled(:)] / args.grid_spacing,...
%                 alpha, beta, max_delta_x, resol_x, max_delta_y, resol_y,...
%                 max_mask);
            
            x_snake = snake_pnts(:,1)*args.grid_spacing;
            y_snake = snake_pnts(:,2)*args.grid_spacing;                    
          
           
        end
        
        candidate_polyfit = interp1(x_m, y_m, candidate_xy(:,1), 'linear');
        candidate_displacements = candidate_xy(:,2) - candidate_polyfit;  %#ok 
        
        if args.do_plot
            figure; 
            %subplot(2,1,1); 
            imgray(location_distribution.D_f);
            plot(x_m/args.grid_spacing, y_m/args.grid_spacing, 'y--');
            %plot(x_i/args.grid_spacing, y_d/args.grid_spacing, 'g');
            cc = 'rgb';
            for i_a = 1:3
                plot(...
                   candidate_xy(candidate_labels==i_a,1) / args.grid_spacing,... 
                   candidate_xy(candidate_labels==i_a,2) / args.grid_spacing,...
                   [cc(i_a) 'x']);
            end
            
            %if ~isempty(y_sampled_s)
            %    plot(x_sampled / args.grid_spacing, y_sampled_s / args.grid_spacing, 'm');
            %end
             
            %subplot(2,1,2); 
            %imgray(max_mask);
            if exist('x_snake', 'var')
                plot(x_snake/args.grid_spacing, y_snake/args.grid_spacing, 'g');
                %plot(x_sampled/args.grid_spacing, y_sampled/args.grid_spacing, 'gx');
            end
        end
    end    
end
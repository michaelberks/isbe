function [] = select_vessels_from_candidates_set(varargin)
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
    'start_i', 1,...
    'end_i', [],...
    'data_dir', [nailfoldroot 'data/rsa_study/test_half/'],...
    'fov_mask_dir', 'fov_masks\',...
    'image_dir', 'images\',...
    'im_rot_dir', 'image_rotations\',...
    'input_dir', 'apex_maps\frog\local_maxima\',...
    'output_dir', 'apex_maps\frog\local_maxima\selected_candidates\',...
    'width_dir', 'predictions\width\rf_regression\297037\',...
    'vessel_centre_dir', [],...
    'min_candidates', 3,...
    'max_candidates', 15,...
    'poly_n', 5,...
    'initial_thresh', 0.3,...
    'class_map', [],...
    'weak_vessel_thresh', 0.3,...
    'strong_vessel_thresh', 0.8,...
    'bad_image_thresh', 0.6,...
    'upper_ydist', -140,...
    'lower_ydist', 90,...
    'angle_discard_thresh', pi/2.5,...
    'angle_keep_thresh', 40,...
    'do_distal_sub', 0, ...
    'do_fill_gaps', 1,...
    'do_final_cull', 1,...
    'do_width', 1,...
    'do_plot', 0,...
    'do_save', 1);
clear varargin;

fov_mask_dir = [args.data_dir args.fov_mask_dir];
image_dir = [args.data_dir args.image_dir];
im_rot_dir = [args.data_dir args.im_rot_dir];
input_dir = [args.data_dir args.input_dir];
output_dir = [args.data_dir args.output_dir];
width_dir = [args.data_dir args.width_dir];

create_folder(output_dir);

if ~isempty(args.vessel_centre_dir)
    vessel_centre_dir = [args.data_dir args.vessel_centre_dir];
else
    vessel_centre_dir = [];
end

im_list = dir([image_dir '*.mat']);
if isempty(args.end_i)
    args.end_i = length(im_list);
end

display(['Processing images ' num2str(args.start_i) ' to ' num2str(args.end_i)]);
for i_im = args.start_i:args.end_i
    display(['Processing image : ' num2str(i_im)]);
    
    im_name = im_list(i_im).name(1:6);
    load([input_dir im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');   
    rot_struc = load([im_rot_dir im_name '.mat'],...
        'centres_r', 'rot_mat', 'ncolsr', 'nrowsr', 'ncols', 'nrows', 'bad');
    
    if args.do_plot
        load([image_dir im_name '.mat']);
        load([fov_mask_dir im_name '_f_mask.mat']);
    
        if args.do_width
            vessel_width = u_load([width_dir im_name '_pred.mat']);
            vessel_width = imfilter(vessel_width, fspecial('gaussian', [10 2]));
            candidate_widths = interp2(vessel_width, candidate_xy(:,1), candidate_xy(:,2)); %#ok
        else
            candidate_widths = 15*ones(size(candidate_xy(:,1))); %#ok
        end
    
        figure; imgray(nailfold);
        caxis([min(nailfold(fov_mask))+5 max(nailfold(fov_mask))-5]);
    end
    
    if sum(candidate_scores > args.initial_thresh) < args.min_candidates
        if args.do_plot
            title({[im_name ': very few candidates found, image may be of ungradeable quality'];...
                'No attempt made to find distal row. Vessel detections at 95% specificty shown'});
        end
        
        kept = candidate_scores > args.bad_image_thresh;
        non_distal = false(size(kept));
        fell_at_the_last = [];
        intermediate_selections = []; %#ok
        candidate_class_probs = []; %#ok
        candidate_class = []; %#ok
        candidate_displacements = []; %#ok
        status = -1; %#ok
        
    elseif rot_struc.bad
        if args.do_plot
            title({[im_name ': mosaic shape suggests frames were not properly registered.'];...
                'No attempt made to find distal row. Vessel detections at 95% specificty shown'});
        end
        
        kept = candidate_scores > args.bad_image_thresh;
        non_distal = false(size(kept));
        fell_at_the_last = [];
        intermediate_selections = []; %#ok
        candidate_class_probs = []; %#ok
        candidate_class = []; %#ok
        candidate_displacements = []; %#ok
        status = -2; %#ok
    else
        
        if ~isempty(vessel_centre_dir)
            load([vessel_centre_dir im_name '_vc.mat']);
        else
            vessel_centre = [];
        end
           
        %Make sure candidates are in descending order of score
        if ~issorted(candidate_scores(end:-1:1))
            [candidate_scores cs_i] = sort(candidate_scores, 'descend');
            candidate_xy = candidate_xy(cs_i,:);
        end

        [kept non_distal intermediate_selections candidate_displacements candidate_class_probs candidate_class] = ...
            select_vessels_from_candidates_old(candidate_xy, candidate_scores, rot_struc,...
                'max_candidates', args.max_candidates,...
                'poly_n', args.poly_n,...
                'initial_thresh', args.initial_thresh,...
                'class_map', args.class_map,...
                'weak_vessel_thresh', args.weak_vessel_thresh,...
                'strong_vessel_thresh', args.strong_vessel_thresh,...
                'upper_ydist', args.upper_ydist,...
                'lower_ydist', args.lower_ydist,...
                'angle_discard_thresh', args.angle_discard_thresh,...
                'angle_keep_thresh', args.angle_keep_thresh,...
                'do_distal_sub', args.do_distal_sub, ...
                'do_fill_gaps', args.do_fill_gaps,...
                'do_final_cull', args.do_final_cull, ...
                'vessel_centre', vessel_centre); %#ok
        status = 0; %#ok
        fell_at_the_last = intermediate_selections(:,3) & ~kept;
        
        if args.do_plot
            title([im_name ': ' num2str(sum(kept)) ' distal apices found']);
        end
    end    
    
    if args.do_plot
        if any(kept)
            cx = candidate_xy(kept,1);
            cy = candidate_xy(kept,2);
            cw = candidate_widths(kept);
            for i_c = 1:length(cx)
                plot_circle(gca, cx(i_c), cy(i_c), cw(i_c)/2, 'r-');
            end

        end
        if any(fell_at_the_last)
            plot(candidate_xy(fell_at_the_last,1), candidate_xy(fell_at_the_last,2), 'bo', 'markersize', 8);
        end
        if any(non_distal)
            plot(candidate_xy(non_distal,1), candidate_xy(non_distal,2), 'gx', 'markersize', 4);
        end
        if any(~kept & ~non_distal)
            plot(candidate_xy(~kept & ~non_distal,1), candidate_xy(~kept & ~non_distal,2), 'cx', 'markersize', 2);
        end
    end
    
    if args.do_save
        save([output_dir im_name '_candidates'],...
            'candidate_xy', 'candidate_scores', 'kept', 'non_distal', 'intermediate_selections',...
            'candidate_displacements', 'candidate_class_probs', 'candidate_class', 'status');   
    end

end


%-------------------------------------------------------------------------
function plot_circle(a_h, cx, cy, r, linestyle)

pi_pts = linspace(0, 2*pi, 20);
circ_x = cx + r*cos(pi_pts);
circ_y = cy + r*sin(pi_pts);
plot(a_h, circ_x, circ_y, linestyle);
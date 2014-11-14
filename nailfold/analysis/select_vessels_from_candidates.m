function [] = ...
    select_vessels_from_candidates(varargin)
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
    {'image_names',...
    'class_map'},...
    'data_dir', [nailfoldroot 'data/rsa_study/master_set/'],...
    'selected_dir', 'selected_apexes',...
    'displacement_dir', 'displacements',...
    'rescore_dir', 'rescores',...
    'prob_dir',             'rf_classification/296655',...
    'width_dir',            'rf_regression/297037',...
    'strong_vessel_thresh', 0.8,...
    'angle_discard_thresh', pi/2.5,...
    'do_final_cull', 1,...
    'do_post_merge', 1,...
    'merge_sigma',   1,...
    'merge_dist_thresh', 60,...
    'merge_connect_thresh', 0.5,...
    'merge_n_pts', 20,...
    'do_plot', 0);
clear varargin;

rescore_dir = [args.data_dir '/' args.rescore_dir '/'];
displacement_dir = [args.data_dir '/' args.displacement_dir '/'];
selected_dir = [args.data_dir '/' args.selected_dir '/'];
if args.do_post_merge
    prob_dir = [args.data_dir 'predictions/detection/' args.prob_dir '/'];
    width_dir = [args.data_dir 'predictions/width/' args.width_dir '/'];
end
create_folder(selected_dir);

num_images = length(args.image_names);

if args.do_post_merge
    g = gaussian_filters_1d(args.merge_sigma);
    g = g / sum(g);
end

%Loop though each image
for i_im = 1:num_images
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);
    
    im_name = args.image_names{i_im} ;  
    
    load([displacement_dir im_name '_dis.mat'], 'candidate_displacements');
    load([rescore_dir im_name '_candidates.mat'], 'candidate_xy', 'candidate_rescores');


    if all(candidate_displacements < inf)
        %Use class probs model to determine whether a candidate should be kept as a
        %distal or non-distal vessel
        status = 0; %#ok
        candidate_class = interp2(args.class_map.x, args.class_map.y, args.class_map.post_class,...
            candidate_rescores, candidate_displacements, 'linear', 0); %#ok
        candidate_class_probs = interp2(args.class_map.x, args.class_map.y, args.class_map.post_probs,...
            candidate_rescores, candidate_displacements, 'linear', 0); %#ok 

        selected_distal = candidate_class == 1;
        selected_non_distal = candidate_class == 2;
        
        if args.do_post_merge
            vessel_prob = u_load([prob_dir im_name '_pred.mat']);
            vessel_prob = conv2(g', g, vessel_prob, 'same');
            vessel_width = u_load([width_dir im_name '_pred.mat']);
            vessel_width = conv2(g', g, vessel_width, 'same');
            
            candidate_widths = interp2(vessel_width,...
                candidate_xy(:,1), candidate_xy(:,2), 'nearest'); %#ok
            
            if args.do_plot
                figure; imgray(vessel_prob);
                for i_can = 1:length(candidate_rescores)
                    if selected_distal(i_can)
                        text(candidate_xy(i_can,1), candidate_xy(i_can,2), num2str(candidate_rescores(i_can)), 'color', 'r');
                        text(candidate_xy(i_can,1), candidate_xy(i_can,2)+10, num2str(candidate_displacements(i_can)), 'color', 'r');
                    elseif selected_non_distal(i_can)
                        text(candidate_xy(i_can,1), candidate_xy(i_can,2), num2str(candidate_rescores(i_can)), 'color', 'g');
                        text(candidate_xy(i_can,1), candidate_xy(i_can,2)+10, num2str(candidate_displacements(i_can)), 'color', 'g');
                    else
                        text(candidate_xy(i_can,1), candidate_xy(i_can,2), num2str(candidate_rescores(i_can)), 'color', 'b');
                        text(candidate_xy(i_can,1), candidate_xy(i_can,2)+10, num2str(candidate_displacements(i_can)), 'color', 'b');
                    end
                end
            end
            
            [merged_with] = post_merge_candidates(candidate_xy, candidate_rescores, candidate_widths, ...
                 vessel_prob, args.merge_dist_thresh, args.merge_connect_thresh, args.merge_n_pts, 0);
             
            for i_can = 1:length(merged_with)
                
                if merged_with(i_can)
                    
                    if selected_distal(i_can)
                        
                        if selected_distal(merged_with(i_can))
                            %Two distals merged together, discard the one with
                            %lower score
                            selected_distal(i_can) = 0;
                            
                        elseif selected_non_distal(merged_with(i_can))
                            %Distal, merged with a non-distal. The
                            %non-distal has higher score, so must have been
                            %reject on its location - do nothing
                            
                        else
                            %Distal, merged with something not selected as
                            %a non-distal, but has a higher vessel score
                            %than it. Swap the distal status
                            selected_distal(i_can) = 0;
                            selected_distal(merged_with(i_can)) = 1;
                        end
                    elseif selected_non_distal(i_can)
                        
                        if selected_distal(merged_with(i_can))
                            %Non-distal, merged with a distal. Discard the
                            %non-distal altogether
                            selected_non_distal(i_can) = 0;
                            
                        elseif selected_non_distal(merged_with(i_can))
                            %2 non-distals merged together, keep the one
                            %with the higher score
                            selected_non_distal(i_can) = 0;
                            
                        else
                            %Non-istal, merged with something not selected as
                            %a non-distal, but has a higher vessel score
                            %than it. Swap the distal status
                            selected_non_distal(i_can) = 0;
                            selected_non_distal(merged_with(i_can)) = 1;
                        end
                    end
                end
            end
            
        end
        
        fell_at_the_last = false(size(candidate_rescores));
        %Finally, we discard any vessels that lie on top of one another
        if args.do_final_cull
            if any(selected_distal)
                %Which candidates belong to the distal row?
                [distal_idx] = select_distal_candidates(...
                    candidate_xy(selected_distal,:), candidate_rescores(selected_distal,:),...
                    10, args.angle_discard_thresh, 0, 0);

                %Work out which ones we've discarded (with indices relative to the
                %original list)
                distal_idx_tf = true(sum(selected_distal),1);
                distal_idx_tf(distal_idx) = 0;
                fell_at_the_last(selected_distal) = distal_idx_tf;

                %
                selected_distal(fell_at_the_last) = 0; %#ok
                selected_non_distal(fell_at_the_last) = 1; %#ok
            end
        end
    else
        %Minimum number of candidates not met to estimate displacements,
        %apply strong threshold
        selected_distal = candidate_rescores > args.strong_vessel_thresh; %#ok
        selected_non_distal = false(size(candidate_rescores)); %#ok
        candidate_class = []; %#ok
        candidate_class_probs = []; %#ok
        fell_at_the_last = []; %#ok
        status = -1; %#ok
    end
    
    save([selected_dir im_name '_sel.mat'],... 
        'selected_distal', 'selected_non_distal', 'candidate_class_probs', 'candidate_class', 'status');
end

    
    
    
    

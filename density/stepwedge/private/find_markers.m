function [markers, region_lims, region] = find_markers( mammo , marker_template, large_image, left_breast, reduction_factor, plot_results, debug)
%FIND_MARKERS : Highest order function used to find markers given the original image
% 
%
% Inputs:
%			mammo				image that you want to find markers in [2D array]
%			plot_results 		1 if you want to display found markers on image, 0 if not [#]
%			extra_plots			1 if you want to display extra intermediate images used for analysis of result, 0 if not [#]
%
% Outputs:
%			image				the original image that has now been cut down in size by removing the white frame around its edge [2D array]
%			markers 			x-y position of the markers found by the software [2D array] 
%								x-position/y-position/strip number/region ID number output
%			x_s 				gives the x-split of the edge strips taken for analysis within a single image (varies with image size) [1D array]
%			y_s 				gives the y-split of the edge strips taken for analysis within a single image (varies with image size) [1D array]
%			delete_range		number of pixel colomns or rows removed from the image edges during "remove_white_frame.m"
%			region 				which set of small regions contains the most found markers [#]
%
% Example:
%
% Notes:	I think it would be best to make the white frame around the image a mask rather than removing it. As the code stands at the moment
%			the marker co-ordinates are in the frame of the cut image and not the original. 
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com

exclusion_zone = 101*reduction_factor;
exclusion_bar = round(200*reduction_factor);
padding = 100;
edge_remove = round(padding*reduction_factor);
exclusion_thresh = 0.4;
corr_thresh = 0.4;
do_separation_check = 0;

markers_per_region = 1;
white_thresh = 250;


[rows cols] = size(mammo);
edge_c = [(1:cols)'; ones(rows-2,1); (1:cols)'; ones(rows-2,1)];
edge_r = [ones(cols,1); (2:rows-1)'; ones(cols,1); (2:rows-1)'];
anon_mask = bwselect(~mammo, edge_c, edge_r, 4);
edge_mask = bwselect(mammo > white_thresh*256, edge_c, edge_r, 4);
valid_mask = ~(edge_mask | anon_mask);

if isempty(large_image)
    large_image = length(mammo)>6000;
end

if  large_image%repeat for small images, not large ones (due to regions)
    R=1;
    markers_to_find = 8;
else
    R=2;
    markers_to_find = 6;
end

if plot_results
    fig_h = zeros(markers_to_find,R);
end
    
for i_r = 1:R %two region definitions
    
    %Prepare image and find the allow_map--------------------------------------
    [region_lims] = region_definitions([rows cols], i_r, large_image, left_breast, padding); %Define regions
    
    %Set up containers for markers,  Each marker has the data:
    % x-position/y-position/strip number/region number
    markers_raw = zeros(7, markers_per_region, markers_to_find);
    
    for i_m = 1:markers_to_find
        
        if i_m > markers_to_find/2
            strip_num = 2;
            pair_num = i_m - markers_to_find/2;
        else
            strip_num = 1;
            pair_num = i_m;
        end

        %Get row,column subscripts of this region
        region_c = region_lims(i_m,1):region_lims(i_m,2);
        region_r = region_lims(i_m,3):region_lims(i_m,4);
        
        %Extract region from image and mask
        valid_mask_i = imresize(valid_mask(region_r, region_c), reduction_factor);
        mammo_i = imresize(mammo(region_r, region_c), reduction_factor);
        
        %Apply normxcorr
        corr_i = mb_normxcorr2(marker_template, mammo_i, [], valid_mask_i);
        
        h_ex = imopen(corr_i > exclusion_thresh, strel('rectangle', [1 exclusion_bar]));
        v_ex = imopen(corr_i > exclusion_thresh, strel('rectangle', [exclusion_bar 1]));
        maxima_mask_i = valid_mask_i & ~(h_ex | v_ex);% & g2d_mask;
        maxima_mask_i([1:edge_remove end-edge_remove+1:end],:) = 0;
        maxima_mask_i(:,[1:edge_remove end-edge_remove+1:end]) = 0;
        
        %Take local maxima or corr map
        [maxima_pos maxima_vals] = local_image_maxima(corr_i, exclusion_zone, maxima_mask_i, corr_thresh,0); 
        max_found = min(markers_per_region, size(maxima_pos,1));
        
        if plot_results
            fig_h(i_m,i_r) = figure('visible', 'off');
            subplot(2,3,1); imgray(mammo_i);
            if ~isempty(maxima_pos)
                plot(maxima_pos(1:max_found,1), maxima_pos(1:max_found,2), 'gx');
            end
            
            subplot(2,3,2); imgray(valid_mask_i);
            if ~isempty(maxima_pos)
                plot(maxima_pos(1:max_found,1), maxima_pos(1:max_found,2), 'gx');
            end
            
            subplot(2,3,3); imgray(corr_i);
            if ~isempty(maxima_pos)
                plot(maxima_pos(1:max_found,1), maxima_pos(1:max_found,2), 'gx');
            end
            
            subplot(2,3,4); imgray(corr_i > exclusion_thresh);
            if ~isempty(maxima_pos)
                plot(maxima_pos(1:max_found,1), maxima_pos(1:max_found,2), 'gx');
            end
            
            subplot(2,3,5); imgray(maxima_mask_i);
            if ~isempty(maxima_pos)
                plot(maxima_pos(1:max_found,1), maxima_pos(1:max_found,2), 'gx');
            end
            
            subplot(2,3,6); imgray(corr_i > corr_thresh);
            if ~isempty(maxima_pos)
                plot(maxima_pos(1:max_found,1), maxima_pos(1:max_found,2), 'gx');
            end
            
        end                
        
        %Correct coordinates to image frame
        if max_found
            markers_raw(1,1:max_found,i_m) = maxima_pos(1:max_found,1)'/reduction_factor + region_lims(i_m,1) - 1;
            markers_raw(2,1:max_found,i_m) = maxima_pos(1:max_found,2)'/reduction_factor + region_lims(i_m,3) - 1;
            markers_raw(3,1:max_found,i_m) = strip_num;
            markers_raw(4,1:max_found,i_m) = pair_num;   
            markers_raw(7,1:max_found,i_m) = maxima_vals(1:max_found,1)';  
        end
    end
    
    markers_raw = markers_raw(:,:)';
	
    %--------------------------------------------------------------------------
	
	%Do separation check
    if do_separation_check && any(markers_raw(:,1))
        temp_markers = single_separation_check(markers_raw);%separation check (discrimination)        
    else
        temp_markers = markers_raw;
    end
	%--------------------------------------------------------------------------
	
    if i_r == 1
        region = 1;
        markers = temp_markers;
        num_markers = sum(temp_markers(:,7));%sum((temp_markers(:,5) == temp_markers(:,6)) & temp_markers(:,1) );
    else
        %Count how many markers we have in region 2
        num_markers2 = sum(temp_markers(:,7));%sum((temp_markers(:,5) == temp_markers(:,6)) & temp_markers(:,1) );
        
        %If more markers in 2, use R2 markers, otherwise stay as we are
        if num_markers2 > num_markers
            
            markers = temp_markers;
            region = 2;
            
            if plot_results
                for i_m = 1:markers_to_find
                   close(fig_h(i_m,1));
                end
                fig_h(:,1) = [];
            end
        elseif plot_results            
            for i_m = 1:markers_to_find
               close(fig_h(i_m,2));
            end
            fig_h(:,2) = [];
        end
    end
end
region_lims = region_definitions([rows cols], region, large_image, left_breast, padding);

%--------------------------------------------------------------------------

%Plot results if require by function input
if plot_results

    figure; imgray(mammo); a1 = gca;
    not_found = ~temp_markers(:,1);
    not_tested_idx = ~not_found & ~markers(:,5) & ~markers(:,6);
    failed_idx = markers(:,5) == 0 & ~(not_tested_idx | not_found);
    partial_idx = markers(:,5) < markers(:,6) & markers(:,5);
    passed_idx = ~(not_tested_idx | failed_idx | partial_idx);
    
    plot(markers(not_tested_idx,1), markers(not_tested_idx,2), 'y^', 'MarkerSize',10);
    plot(markers(failed_idx,1), markers(failed_idx,2), 'rx', 'MarkerSize',10);
    plot(markers(partial_idx,1), markers(partial_idx,2), 'b+', 'MarkerSize',10);
    plot(markers(passed_idx,1), markers(passed_idx,2), 'co', 'MarkerSize',10);
    
    for i_m = 1:markers_to_find
        plot(a1, region_lims(i_m,[1 2 2 1 1]), region_lims(i_m,[3 3 4 4 3]));
        if debug || failed_idx(i_m)
            set(fig_h(i_m), 'visible', 'on');
        else
            close(fig_h(i_m));
        end
    end
end
%--------------------------------------------------------------------------



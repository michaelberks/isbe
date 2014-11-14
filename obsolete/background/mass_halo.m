function [mean_halo var_halo seg_area] = mass_halo(mass, mass_model, shape_idx, num_segs, n1, n2, varargin)

    args = u_packargs(varargin, 0,...
    'startingPts', 1,...
    'if_plot', 0);
    
    im_loci = mass.mass_ROI;
    [rows cols] = size(im_loci);
    clear mass;
    
    mass_outline = mass_model.P_shape(:,1:15)*mass_model.B_shape(1:15,shape_idx)...
            + mass_model.mean_shape';
    mass_outline = reshape(mass_outline, [], 2);
    
    rotation = mass_model.rotations(:,:,shape_idx);
    translation = repmat(mass_model.translations(shape_idx,:),...
        size(mass_outline, 1), 1);
    origin = mass_model.origins(shape_idx);
    scale = mass_model.X_scale(shape_idx);
%     mass_outline = (circshift(mass_outline/scale, 1-origin)...
%             - translation)*inv(rotation);
        
    inner_ring = (circshift(mass_outline/(scale-n1), 1-origin)...
        - translation)*inv(rotation);
    outer_ring = (circshift(mass_outline/(scale -n2), 1-origin)...
        - translation)*inv(rotation);
    
    if args.if_plot
        figure('WindowStyle', 'docked'); hold on;
        imagesc(im_loci); axis('image'); colormap(gray(256));
        plot(mass_model.shapes_unaligned(shape_idx, 1:end/2),...
            mass_model.shapes_unaligned(shape_idx,end/2+1:end), 'b:');
        plot(inner_ring(:,1), inner_ring(:,2), 'gx');
        plot(outer_ring(:,1), outer_ring(:,2), 'rx');
        
    end
    
%     c_to_i = (outer_ring(o_i(1),:) - mass_centroid) /...
%         sqrt(sum((outer_ring(o_i(1),:) - mass_centroid).^2));
% 
%     p = repmat(outer_ring(o_i(1),:) - n2*c_to_i,...
%         size(inner_ring, 1), 1);
% 
%     [md m_idx] = min(sum((inner_ring - p).^2, 2));

    mean_halo = zeros(num_segs, length(args.startingPts));
    var_halo = zeros(num_segs, length(args.startingPts));
    seg_area = zeros(num_segs, length(args.startingPts));
    
    for ii = 1:length(args.startingPts)
        start_pt = args.startingPts(ii);
        o_ring = circshift(outer_ring, start_pt-1);
        o_ring(end+1,:) = o_ring(1,:);
        i_ring = circshift(inner_ring, start_pt-1);
        i_ring(end+1,:) = i_ring(1,:);
        
        i_i = round(linspace(1, size(i_ring, 1), num_segs+1));
        o_i = i_i;
        
%         o_ring = circshift(outer_ring, start_pt-1);
%         o_ring(end+1,:) = o_ring(1,:);
%         
%         o_i = round(linspace(1, size(o_ring, 1), num_segs+1));
%         
%         [md m_idx] = min(sum((inner_ring - ...
%             repmat(o_ring(o_i(1),:), size(inner_ring,1), 1)).^2, 2));
% 
%         i_ring = circshift(inner_ring, -m_idx);
%         i_ring(end+1,:) = i_ring(1,:);
%         
%         i_i(1) = size(i_ring, 1);
%         i_i(num_segs+1) = 1;
% 
%     %     if args.if_plot; 
%     %         plot([outer_ring(o_i(1),1), p(1,1)],...
%     %             [outer_ring(o_i(1),2), p(1,2)], 'cx'); 
%     %     end
% 
%         for jj = 2:num_segs
%     %         c_to_i = (outer_ring(o_i(jj),:) - mass_centroid) /...
%     %             sqrt(sum((outer_ring(o_i(jj),:) - mass_centroid).^2));
%     % 
%     %         p = repmat(outer_ring(o_i(jj),:) - n2*c_to_i,...
%     %             size(inner_ring, 1), 1);
%     %         if args.if_plot; 
%     %             plot([outer_ring(o_i(jj),1), p(1,1)],...
%     %                 [outer_ring(o_i(jj),2), p(1,2)], 'cx'); 
%     %         end
%     %         
%     %         [md m_idx] = min(sum((inner_ring - p).^2, 2));
% 
%             [md m_idx] = min(sum((i_ring -...
%                 repmat(o_ring(o_i(jj),:), size(i_ring,1),1) ).^2, 2));
%             i_i(jj) = m_idx;
%         end
        %i_i = fliplr(i_i);

        for jj = 1:num_segs 

            seg_poly = [i_ring(i_i(jj), :);...
                        o_ring(o_i(jj):o_i(jj+1),:);...
                        i_ring(i_i(jj+1):-1:i_i(jj),:)];

            %display(['segment area = ',...
            %    num2str(polyarea(seg_poly(:,1), seg_poly(:,2)))]);
            mask = poly2mask(seg_poly(:,1),...
                seg_poly(:,2), rows, cols);
            mean_halo(jj, ii) = mean(im_loci(mask));
            var_halo(jj, ii) = var(im_loci(mask));
            seg_area(jj, ii) = polyarea(seg_poly(:,1), seg_poly(:,2));
            
            if args.if_plot
                plot(seg_poly(:,1), seg_poly(:,2), 'r');
%                 plot([o_ring(o_i(jj),1), mass_centroid(1)],...
%                     [o_ring(o_i(jj),2), mass_centroid(2)], 'g');
%                 plot([i_ring(i_i(jj),1), mass_centroid(1)],...
%                     [i_ring(i_i(jj),2), mass_centroid(2)], 'y');
            end
        end
    end
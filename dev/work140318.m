data_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\';
centre_dir = [data_dir 'vessel_centres\full_centres\'];
apex_dir = [data_dir 'apex_clusters_merged\'];

load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat', 'miccai_selection');
im_names = image_id_data.im_names(miccai_selection.validation);

grid_spacing = 8;
%%
for i_im = 1:20

    %Load in the gt data
    load([apex_dir im_names{i_im} '_apex_clusters.mat']);

    %Load in the image size
    load([centre_dir im_names{i_im} '_vc.mat'], 'nrows', 'ncols');
    
    %For each marker
    num_vessels = length(vessels.cluster_members);
    
    figure;
    marker_dist = cell(2,1);
    for i_ma = 1:2
        ma = vessels.markers(i_ma);
        marked_by = false(num_vessels,1);
        ma_xy = zeros(num_vessels,2);
        ma_shapes = cell(num_vessels,1); 
        for i_ve = 1:num_vessels
            ma_i = find(ma == vessels.cluster_members{i_ve});
            if ~isempty(ma_i) && ~strcmpi(vessels.cluster_shapes{i_ve}{ma_i}, 'NonDistal')
                marked_by(i_ve) = 1;
                ma_xy(i_ve,:) = vessels.cluster_xy{i_ve}(ma_i,:) / 2;
                ma_shapes{i_ve} = vessels.cluster_shapes{i_ve}{ma_i};
            end
        end
        ma_xy(~marked_by,:) = [];
        ma_shapes(~marked_by,:) = [];
        %display(ma_shapes);
        
        xx = 1:grid_spacing:ncols;
        yy = (1:grid_spacing:nrows)';
        grid_x = repmat(xx, length(yy), 1);
        grid_y = repmat(yy, 1, length(xx));
   
        %Compute weighted kernel estimates of the spatial distribution of
        %candidates over this grid
        [location_distribution] = build_2d_kernel_distribution(...
           ma_xy,...
           [grid_x(:) grid_y(:)]);
   
        marker_dist{i_ma} = reshape(location_distribution.D_f, size(grid_x));
        
        
           
%         if i_ma == 1
%             subplot(2,1,i_ma); 
%             imgray(marker_dist{i_ma});
%             plot(ma_xy(:,1)/grid_spacing, ma_xy(:,2)/grid_spacing, 'rx');
%         else
%             plot(ma_xy(:,1)/grid_spacing, ma_xy(:,2)/grid_spacing, 'gx');
%             subplot(2,1,i_ma);
%             imgray(marker_dist{1} - marker_dist{2});
%             
%             max_val = max(max(marker_dist{1}(:)), max(marker_dist{2}(:)));
%             caxis([-max_val max_val]);
%         end

        subplot(2,1,i_ma); 
        imgray(marker_dist{i_ma});
        plot(ma_xy(:,1)/grid_spacing, ma_xy(:,2)/grid_spacing, 'rx');
        
        if i_ma == 2
            [D_b] = compute_bhattacharyya_distance(...
                marker_dist{1}, marker_dist{2});
        
            subplot(2,1,1);
            title(['D_b = ' num2str(D_b, 3)]);
        end

    end
end
%%
apex_areas = [0 0];
image_widths = 0;
for i_im = 1:length(im_names)

    %Load in the gt data
    load([apex_dir im_names{i_im} '_apex_clusters.mat']);

    %Load in the image size
    load([centre_dir im_names{i_im} '_vc.mat'], 'nrows', 'ncols');
    image_widths = image_widths + ncols;
    
    %For each marker
    num_vessels = length(vessels.cluster_members);
    
    marker_dist = cell(2,1);
    for i_ma = 1:2
        ma = vessels.markers(i_ma);
        marked_by = false(num_vessels,1);
        ma_xy = zeros(num_vessels,2);
        ma_widths = zeros(num_vessels,1);
        ma_shapes = cell(num_vessels,1); 
        for i_ve = 1:num_vessels
            ma_i = find(ma == vessels.cluster_members{i_ve});
            if ~isempty(ma_i) %&& ~strcmpi(vessels.cluster_shapes{i_ve}{ma_i}, 'NonDistal')
                marked_by(i_ve) = 1;
                ma_xy(i_ve,:) = vessels.cluster_xy{i_ve}(ma_i,:) / 2;
                ma_shapes{i_ve} = vessels.cluster_shapes{i_ve}{ma_i};
                ma_widths(i_ve,:) = vessels.cluster_widths{i_ve}(ma_i) / 4;
            end
        end
        ma_xy(~marked_by,:) = [];
        ma_shapes(~marked_by,:) = [];
        ma_widths(~marked_by,:) = [];

        apex_areas(i_ma) = apex_areas(i_ma) + sum(pi*ma_widths.^2);        

    end
end

% 2d inter-level coefficient stuff
% Build oriented stepdedges at angles 0:pi/20:19pi/20;
angles = 0:pi/20:19*pi/20;
pos_edge_256 = [ones(128,256); 2*ones(128,256)];
neg_edge_256 = [2*ones(128,256); ones(128,256)];

pos_edges_128 = zeros(128,128,length(angles));
neg_edges_128 = zeros(128,128,length(angles));

pos_lines_256 = ones(256); pos_lines_256(128,:) = 2;
neg_lines_256 = 2*ones(256); neg_lines_256(128,:) = 1;

pos_lines_128 = zeros(128,128,length(angles));
neg_lines_128 = zeros(128,128,length(angles));

for ii = 1:length(angles);
    
    temp1 = imrotate(pos_edge_256, 180*angles(ii)/pi, 'crop');
    temp2 = imrotate(neg_edge_256, 180*angles(ii)/pi, 'crop');
    temp3 = imrotate(pos_lines_256, 180*angles(ii)/pi, 'crop');
    temp4 = imrotate(neg_lines_256, 180*angles(ii)/pi, 'crop');
    
    pos_edges_128(:,:,ii) = temp1(65:192, 65:192);
    neg_edges_128(:,:,ii) = temp2(65:192, 65:192);
    pos_lines_128(:,:,ii) = temp3(65:192, 65:192);
    neg_lines_128(:,:,ii) = temp4(65:192, 65:192);

    clear temp*;
    
end

%%
o = 1;
for edge_angle = [3 6 9 13 16 19] %edges aligned at sub-band orientations
    
    dt_pos_edge = dtwavexfm2(pos_edges_128(:,:,edge_angle), 5, 'near_sym_b','qshift_b');
    %[ilp_pos_edge icp_pos_edge] = mb_dual_tree_transform(dt_pos_edge);
    dt_pos_edge = dt_to_full_image(dt_pos_edge);
    ilp_pos_edge = ilp_same_size(dt_pos_edge);
    
    dt_neg_edge = dtwavexfm2(neg_edges_128(:,:,edge_angle), 5, 'near_sym_b','qshift_b');
    dt_neg_edge = dt_to_full_image(dt_neg_edge);
    ilp_neg_edge = ilp_same_size(dt_neg_edge);
    %[ilp_neg_edge icp_neg_edge] = mb_dual_tree_transform(dt_neg_edge);
    
    figure('Name', 'Edge');
    for level = 2:4
        dim = 2^level;
        st = (dim+1)/2;
        
        %[max_icp_pos] = max(icp_pos_edge{level}, [], 3);
        %[max_icp_neg] = max(icp_neg_edge{level}, [], 3);
        
        [max_ilp_pos band_idx_pos] = max(ilp_pos_edge(:,:,:,level), [], 3);
        [max_ilp_neg band_idx_neg] = max(ilp_neg_edge(:,:,:,level), [], 3);     
        
        subplot(2,3,level-1);
        %imagesc(pos_edges_128(:,:,edge_angle)); axis image; colormap(gray(256)); hold on;
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_pos), -imag(max_ilp_pos), 'r');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_pos_opp), -imag(max_ilp_pos_opp), 'b');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_icp_pos), -imag(max_icp_pos), 'g');
        %image(complex2rgb(ilp_pos_edge(:,:,o,level))); axis image;
        image(complex2rgb(max_ilp_pos)); axis image;
        subplot(2,3,level+2);
        %imagesc(neg_edges_128(:,:,edge_angle)); axis image; colormap(gray(256)); hold on;
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_neg), -imag(max_ilp_neg), 'r');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_neg_opp), -imag(max_ilp_neg_opp), 'b');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_icp_neg), -imag(max_icp_neg), 'g');
        %image(complex2rgb(ilp_neg_edge(:,:,o,level))); axis image;
        image(complex2rgb(max_ilp_neg)); axis image;    
    end
    
    dt_pos_line = dtwavexfm2(pos_lines_128(:,:,edge_angle), 5, 'near_sym_b','qshift_b');
    %[ilp_pos_line icp_pos_line] = mb_dual_tree_transform(dt_pos_line);  
    dt_pos_line = dt_to_full_image(dt_pos_line);
    ilp_pos_line = ilp_same_size(dt_pos_line);
    
    dt_neg_line = dtwavexfm2(neg_lines_128(:,:,edge_angle), 5, 'near_sym_b','qshift_b');
    [ilp_neg_line icp_neg_line] = mb_dual_tree_transform(dt_neg_line);
    dt_neg_line = dt_to_full_image(dt_neg_line);
    ilp_neg_line = ilp_same_size(dt_neg_line);
    
    figure('Name', 'Line');
    for level = 2:4
        dim = 2^level;
        st = (dim+1)/2;
        %[max_icp_p] = max(icp_pos_line{level}, [], 3);
        %[max_icp_n] = max(icp_neg_line{level}, [], 3);
        
        [max_ilp_pos band_idx_p] = max(ilp_pos_line(:,:,:,level), [], 3);
        [max_ilp_neg band_idx_n] = max(ilp_neg_line(:,:,:,level), [], 3);
        
        subplot(2,3,level-1);
        %figure;
        %imagesc(pos_lines_128(:,:,edge_angle)); axis image; colormap(gray(256)); hold on;
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_p), -imag(max_ilp_p), 0, 'r');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_p_opp), -imag(max_ilp_p_opp), 0, 'b:');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_icp_p), -imag(max_icp_p), 'g:');
        %image(complex2rgb(ilp_pos_line(:,:,o,level))); axis image;
        image(complex2rgb(max_ilp_pos)); axis image;
        subplot(2,3,2+level);
        %figure;
        %imagesc(neg_lines_128(:,:,edge_angle)); axis image; colormap(gray(256)); hold on;
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_n), -imag(max_ilp_n), 0, 'r');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_n_opp), -imag(max_ilp_n_opp), 0, 'b:');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_icp_n), -imag(max_icp_n), 'g:');
        %image(complex2rgb(ilp_neg_line(:,:,o,level))); axis image;
        image(complex2rgb(max_ilp_neg)); axis image;
    end
    o = o+1;
end
%%
o = 1;
for edge_angle = [3 6 9 13 16 19] %edges aligned at sub-band orientations
    

    [phase_map_pos_edge] = dt_phase_congruency(pos_edges_128(:,:,edge_angle), 2, 5);
    %[phase_map_neg_edge] = dt_phase_congruency(neg_edges_128(:,:,edge_angle), 2, 5);
    
    figure('Name', 'Edge');   
    for band = 1:6
        subplot(2,3,band);
        image(complex2rgb(phase_map_pos_edge(:,:,band))); axis image;
    end
    
    [phase_map_pos_line] = dt_phase_congruency(pos_lines_128(:,:,edge_angle), 2, 5);
    %[phase_map_neg_line] = dt_phase_congruency(neg_lines_128(:,:,edge_angle), 2, 5);
    
    figure('Name', 'Line');
    for band = 1:6
        subplot(2,3,band);
        image(complex2rgb(phase_map_pos_line(:,:,band))); axis image;
    end
    
    o = o+1;
end
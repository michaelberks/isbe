% 2d inter-level coefficient stuff
% Build oriented stepdedges at angles 0:pi/20:19pi/20;
angles = 0:pi/20:19*pi/20;
pos_edge_256 = [ones(128,256); 2*ones(128,256)];
neg_edge_256 = [2*ones(128,256); ones(128,256)];

pos_edges_128 = zeros(128,128,length(angles));
neg_edges_128 = zeros(128,128,length(angles));

for ii = 1:length(angles);
    
    temp = imrotate(pos_edge_256, 180*angles(ii)/pi, 'crop');
    pos_edges_128(:,:,ii) = temp(65:192, 65:192);
    clear temp;
    temp = imrotate(neg_edge_256, 180*angles(ii)/pi, 'crop');
    neg_edges_128(:,:,ii) = temp(65:192, 65:192);
    clear temp;
    %figure; image(step_edges_128(:,:,ii)); colormap(gray(256)); axis image;
end
%%

for edge_angle = 1:20;%[3 6 9 13 16 19] %edges aligned at sub-band orientations
    
    [icp_pos_edge] = mb_full_image_icp(pos_edges_128(:,:,edge_angle), 2, 4);
    [icp_neg_edge] = mb_full_image_icp(pos_edges_128(:,:,edge_angle), 2, 4);
    
%     dt_pos_edge = dtwavexfm2(pos_edges_128(:,:,edge_angle), 5, 'near_sym_b','qshift_b');
%     [ilp_pos_edge icp_pos_edge] = mb_dual_tree_transform(dt_pos_edge);
%     
%     dt_neg_edge = dtwavexfm2(neg_edges_128(:,:,edge_angle), 5, 'near_sym_b','qshift_b');
%     [ilp_neg_edge icp_neg_edge] = mb_dual_tree_transform(dt_neg_edge);
    
    figure('Name', 'Edge');
    for level = 1:3
        
        [max_icp_pos] = max(icp_pos_edge(:,:,:,level), [], 3);
        [max_icp_neg] = max(icp_neg_edge(:,:,:,level), [], 3);
%         [max_icp_pos] = max(icp_pos_edge{level}, [], 3);
%         [max_icp_neg] = max(icp_neg_edge{level}, [], 3);
        
        subplot(2,3,level);
        image(complex2rgb(abs(max_icp_pos) .* exp(2*i*angle(max_icp_pos)))); axis image;
        subplot(2,3,level+3);
        image(complex2rgb(abs(max_icp_neg) .* exp(2*i*angle(max_icp_neg)))); axis image;    
    end
end
%%
normal_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
target_region = (double(imread(['C:\isbe\dev\background\images\normal1024\', normal_list(24).name])));
[r c] = size(target_region);
[target_orientation_map] = mb_dt_orientation_map('Image', target_region);
    
[xx yy] = meshgrid(1:c, 1:r');
circle = (xx - c/2).^2 + (yy - r/2).^2 < 512^2;

[tout rout hist_0] = weighted_complex_rose(target_orientation_map(circle), 360);
figure; polar(tout, rout);
for aa = 30:30:360

    pad_vec = ceil(0.5*(sqrt(r^2 + c^2) - [r c]));
    image_pad = padarray(target_region, pad_vec, 'symmetric');
    
    %Rotate
    angle_region = imrotate(image_pad, aa, 'bilinear', 'crop');
    
    %Take central region
    angle_region = angle_region(pad_vec(1)+(1:r), pad_vec(2)+(1:c));
    
    %Compute orientation map and histogram
    [angle_orientation_map] = mb_dt_orientation_map('Image', angle_region);
    [tout rout hist_a] = weighted_complex_rose(angle_orientation_map(circle), 360);
    figure; polar(tout, rout);
    
    orientation_corr = zeros(360,1);
    for o = 1:360
        orientation_corr(o) = ...
            corr2(hist_a, circshift(hist_0,o));
    end
    
    %Take rotation angle as the shift that produces maximum correlation
    [dummy rotation_angle] = max(orientation_corr);
    clear dummy;
    display(rotation_angle);
end
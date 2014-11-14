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
for edge_angle = 1:20%[3 6 9 13 16 19] %edges aligned at sub-band orientations
    
    dt_pos_edge = dtwavexfm2(pos_edges_128(:,:,edge_angle), 5, 'near_sym_b','qshift_b');
    [ilp_pos_edge icp_pos_edge] = mb_dual_tree_transform(dt_pos_edge);
    %icp_pos_edge = inter_coefficient_product(dt_pos_edge);
    
    dt_neg_edge = dtwavexfm2(neg_edges_128(:,:,edge_angle), 5, 'near_sym_b','qshift_b');
    [ilp_neg_edge icp_neg_edge] = mb_dual_tree_transform(dt_neg_edge);
    %icp_neg_edge = inter_coefficient_product(dt_neg_edge);
    
    figure;
    for level = 1:3
        dim = 2^level;
        st = (dim+1)/2;
        
        [max_icp_pos] = max(icp_pos_edge{level}, [], 3);
        [max_icp_neg] = max(icp_neg_edge{level}, [], 3);
        
        [max_ilp_pos band_idx_pos] = max(ilp_pos_edge{level}, [], 3);
        [max_ilp_neg band_idx_neg] = max(ilp_neg_edge{level}, [], 3);
        
        band_idx_pos_opp = mod(band_idx_pos+2, 6) + 1;
        band_idx_neg_opp = mod(band_idx_neg+2, 6) + 1;
        
        max_ilp_pos_opp = zeros(size(ilp_pos_edge{level}(:,:,1)));
        max_ilp_neg_opp = zeros(size(ilp_neg_edge{level}(:,:,1)));
        
        for band = 1:6
            temp = ilp_pos_edge{level}(:,:,band);
            max_ilp_pos_opp(band_idx_pos_opp==band) = temp(band_idx_pos_opp==band);
            
            temp = ilp_neg_edge{level}(:,:,band);
            max_ilp_neg_opp(band_idx_neg_opp==band) = temp(band_idx_neg_opp==band);
        end
        
        subplot(2,3,level);
        imagesc(pos_edges_128(:,:,edge_angle)); axis image; colormap(gray(256)); hold on;
        quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_pos), -imag(max_ilp_pos), 'r');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_pos_opp), -imag(max_ilp_pos_opp), 'b');
        quiver(st:dim:129-st, st:dim:129-st, real(max_icp_pos), -imag(max_icp_pos), 'g');
%         image(complex2rgb(max_icp_pos)); axis image;
        subplot(2,3,level+3);
        imagesc(neg_edges_128(:,:,edge_angle)); axis image; colormap(gray(256)); hold on;
        quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_neg), -imag(max_ilp_neg), 'r');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_neg_opp), -imag(max_ilp_neg_opp), 'b');
        quiver(st:dim:129-st, st:dim:129-st, real(max_icp_neg), -imag(max_icp_neg), 'g');
%         image(complex2rgb(max_icp_neg)); axis image;
            
    end
    %display(round(20*angle(max_icp_pos)/pi));
    %display(abs(max_ilp_pos));
    %display(abs(max_ilp_pos_opp));
end
%%
pos_lines_256 = ones(256); pos_lines_256(128,:) = 2;
neg_lines_256 = 2*ones(256); neg_lines_256(128,:) = 1;

%Make some slightly different lines, with a gradual profile
% pos_lines_256 = ones(256); 
% pos_lines_256([126 130],97:160) = 2;
% pos_lines_256([127 129],97:160) = 3;
% pos_lines_256(128,97:160) = 3.5;
% neg_lines_256 = 3.5*ones(256); 
% neg_lines_256([126 130],97:160) = 3;
% neg_lines_256([127 129]:160) = 2;
% neg_lines_256(128,97:160) = 1;

pos_lines_128 = zeros(128,128,length(angles));
neg_lines_128 = zeros(128,128,length(angles));

for ii = 1:length(angles);
    
    temp1 = imrotate(pos_lines_256, 180*angles(ii)/pi, 'crop');
    temp2 = imrotate(neg_lines_256, 180*angles(ii)/pi, 'crop');
    pos_lines_128(:,:,ii) = temp1(65:192, 65:192);
    neg_lines_128(:,:,ii) = temp2(65:192, 65:192);
    %figure;
    %subplot(1,2,1); imagesc(pos_lines_128(:,:,ii)); colormap(gray(256)); axis image;
    %subplot(1,2,2); imagesc(neg_lines_128(:,:,ii)); colormap(gray(256)); axis image;
    clear temp*;
    
end
%%
for line_angle = 1:20 %[3 6 9 13 16 19] %lines aligned at sub-band orientations
    
    dt_pos_line = dtwavexfm2(pos_lines_128(:,:,line_angle), 5, 'near_sym_b','qshift_b');
    [ilp_pos_line icp_pos_line] = mb_dual_tree_transform(dt_pos_line);
    %icp_pos_line = inter_coefficient_product4(dt_pos_line);
    
    dt_neg_line = dtwavexfm2(neg_lines_128(:,:,line_angle), 5, 'near_sym_b','qshift_b');
    [ilp_neg_line icp_neg_line] = mb_dual_tree_transform(dt_neg_line);
    %icp_neg_line = inter_coefficient_product4(dt_neg_line);
    
    figure;
    for level = 2:4
        dim = 2^level;
        st = (dim+1)/2;
        [max_icp_p] = max(icp_pos_line{level}, [], 3);
        [max_icp_n] = max(icp_neg_line{level}, [], 3);
        
        [max_ilp_p band_idx_p] = max(ilp_pos_line{level}, [], 3);
        [max_ilp_n band_idx_n] = max(ilp_neg_line{level}, [], 3);
        
        band_idx_p_opp = mod(band_idx_p+2, 6) + 1;
        band_idx_n_opp = mod(band_idx_n+2, 6) + 1;
        
        max_ilp_p_opp = zeros(size(ilp_pos_line{level}(:,:,1)));
        max_ilp_n_opp = zeros(size(ilp_neg_line{level}(:,:,1)));
        
        for band = 1:6
            temp_p = ilp_pos_line{level}(:,:,band);
            temp_n = ilp_neg_line{level}(:,:,band);
            max_ilp_p_opp(band_idx_p_opp==band) = temp_p(band_idx_p_opp==band);
            max_ilp_n_opp(band_idx_n_opp==band) = temp_n(band_idx_n_opp==band);
        end
        
        %scale the ILP coefficients for quiver (so we can use absolute scaling
        %function in quiver)
        sf_p = dim / max(abs(max_ilp_p(:)));
        sf_n = dim / max(abs(max_ilp_n(:)));
        
        max_ilp_p = max_ilp_p*sf_p;
        max_ilp_p_opp = max_ilp_p_opp*sf_p;
        max_ilp_n = max_ilp_n*sf_n;
        max_ilp_n_opp = max_ilp_n_opp*sf_n;
        
        subplot(2,3,level-1);
        %figure;
        imagesc(pos_lines_128(:,:,line_angle)); axis image; colormap(gray(256)); hold on;
        quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_p), -imag(max_ilp_p), 0, 'r');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_p_opp), -imag(max_ilp_p_opp), 0, 'b:');
        quiver(st:dim:129-st, st:dim:129-st, real(max_icp_p), -imag(max_icp_p), 'g:');
        subplot(2,3,2+level);
        %figure;
        imagesc(neg_lines_128(:,:,line_angle)); axis image; colormap(gray(256)); hold on;
        quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_n), -imag(max_ilp_n), 0, 'r');
        %quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_n_opp), -imag(max_ilp_n_opp), 0, 'b:');
        quiver(st:dim:129-st, st:dim:129-st, real(max_icp_n), -imag(max_icp_n), 'g:');    
    end
    
%     display(abs(max_ilp_p));
%     display(abs(max_ilp_p_opp));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Filled circled
[xx yy] = meshgrid(1:128, 1:128);
pos_circle_128 = ((xx - 64).^2 + (yy - 64).^2 < 25) + 1;
figure; imagesc(pos_circle_128); colormap(gray); axis image;

dt_pos_circle = dtwavexfm2(pos_circle_128, 5, 'near_sym_b','qshift_b');
[ilp_pos_circle icp_pos_circle] = mb_dual_tree_transform(dt_pos_circle);
%%

for level = 1:4
    dim = 2^level;
    st = (dim+1)/2;
    [max_icp_p] = max(icp_pos_circle{level}, [], 3);
    [max_ilp_p band_idx_p] = max(ilp_pos_circle{level}, [], 3);
    
    band_idx_p_opp = mod(band_idx_p+2, 6) + 1;
    max_ilp_p_opp = zeros(size(ilp_pos_circle{level}(:,:,1)));

    for band = 1:6
        temp_p = ilp_pos_circle{level}(:,:,band);
        max_ilp_p_opp(band_idx_p_opp==band) = temp_p(band_idx_p_opp==band);
    end
    
    sf_p = dim / max(abs(max_ilp_p(:)));
    max_ilp_p = max_ilp_p*sf_p;
    max_ilp_p_opp = max_ilp_p_opp*sf_p;
        
    figure;
    imagesc(pos_circle_128); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_p), -imag(max_ilp_p), 0, 'r');
    quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_p_opp), -imag(max_ilp_p_opp), 0, 'b:');
    quiver(st:dim:129-st, st:dim:129-st, real(max_icp_p), -imag(max_icp_p), 'g:');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Open circle
[xx yy] = meshgrid(1:128, 1:128);
open_circle_128 = (((xx - 64).^2 + (yy - 64).^2 > 30^2) & ((xx - 64).^2 + (yy - 64).^2 <= 32^2)) + 1;
figure; imagesc(open_circle_128); colormap(gray); axis image;

dt_open_circle = dtwavexfm2(open_circle_128, 5, 'near_sym_b','qshift_b');
[ilp_open_circle icp_open_circle] = mb_dual_tree_transform(dt_open_circle);

%%

for level = 1:4
    dim = 2^level;
    st = (dim+1)/2;
    [max_icp_p] = max(icp_open_circle{level}, [], 3);
    [max_ilp_p band_idx_p] = max(ilp_open_circle{level}, [], 3);
    
    band_idx_p_opp = mod(band_idx_p+2, 6) + 1;
    max_ilp_p_opp = zeros(size(ilp_open_circle{level}(:,:,1)));

    for band = 1:6
        temp_p = ilp_open_circle{level}(:,:,band);
        max_ilp_p_opp(band_idx_p_opp==band) = temp_p(band_idx_p_opp==band);
    end
    
    sf_p = dim / max(abs(max_ilp_p(:)));
    max_ilp_p = max_ilp_p*sf_p;
    max_ilp_p_opp = max_ilp_p_opp*sf_p;
        
    figure;
    imagesc(open_circle_128); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_p), -imag(max_ilp_p), 0, 'r');
    quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_p_opp), -imag(max_ilp_p_opp), 0, 'b:');
    quiver(st:dim:129-st, st:dim:129-st, real(max_icp_p), -imag(max_icp_p), 'g:');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try the new constructs on Lenna
lenna = u_load('C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna.mat');
lenna90 = imrotate(lenna, 90);
lenna180 = imrotate(lenna, 180);
lenna270 = imrotate(lenna, 270);

% compute dual-tree for each image
dt_lenna = dtwavexfm2(lenna, 5, 'near_sym_b','qshift_b');
dt_lenna90 = dtwavexfm2(lenna90, 5, 'near_sym_b','qshift_b');
dt_lenna180 = dtwavexfm2(lenna180, 5, 'near_sym_b','qshift_b');
dt_lenna270 = dtwavexfm2(lenna270, 5, 'near_sym_b','qshift_b');

%
% compute ILP and ICP coefficients for each image
[ilp_0 icp_0] = mb_dual_tree_transform(dt_lenna);
[ilp_90 icp_90] = mb_dual_tree_transform(dt_lenna90);
[ilp_180 icp_180] = mb_dual_tree_transform(dt_lenna180);
[ilp_270 icp_270] = mb_dual_tree_transform(dt_lenna270);

for level = 2:4
    
    dim = 2^level;
    st = (dim+1)/2;
        
    [max_icp_0 band_idx_0] = max(icp_0{level}, [], 3);
    [max_icp_90 band_idx_90] = max(icp_90{level}, [], 3);
    [max_icp_180 band_idx_180] = max(icp_180{level}, [], 3);
    [max_icp_270 band_idx_270] = max(icp_270{level}, [], 3);
    
    max_ilp_0 = zeros(size(ilp_0{level}(:,:,1)));
    max_ilp_90 = zeros(size(ilp_0{level}(:,:,1)));
    max_ilp_180 = zeros(size(ilp_0{level}(:,:,1)));
    max_ilp_270 = zeros(size(ilp_0{level}(:,:,1)));    

    for band = 1:6
        temp_0 = ilp_0{level}(:,:,band);
        temp_90 = ilp_90{level}(:,:,band);
        temp_180 = ilp_180{level}(:,:,band);
        temp_270 = ilp_270{level}(:,:,band);
        
        max_ilp_0(band_idx_0==band) = temp_0(band_idx_0==band);
        max_ilp_90(band_idx_90==band) = temp_90(band_idx_90==band);
        max_ilp_180(band_idx_180==band) = temp_180(band_idx_180==band);
        max_ilp_270(band_idx_270==band) = temp_270(band_idx_270==band);
    end
    
    % Show the images
    figure;
    subplot(2,2,1); imagesc(lenna); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:257-st, st:dim:257-st, real(max_ilp_0), -imag(max_ilp_0), 'r');
    quiver(st:dim:257-st, st:dim:257-st, real(max_icp_0), -imag(max_icp_0));
    
    subplot(2,2,2); imagesc(lenna90); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:257-st, st:dim:257-st, real(max_ilp_90), -imag(max_ilp_90), 'r');
    quiver(st:dim:257-st, st:dim:257-st, real(max_icp_90), -imag(max_icp_90));
    
    subplot(2,2,3); imagesc(lenna180); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:257-st, st:dim:257-st, real(max_ilp_180), -imag(max_ilp_180), 'r');
    quiver(st:dim:257-st, st:dim:257-st, real(max_icp_180), -imag(max_icp_180));
    
    subplot(2,2,4); imagesc(lenna270); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:257-st, st:dim:257-st, real(max_ilp_270), -imag(max_ilp_270), 'r');
    quiver(st:dim:257-st, st:dim:257-st, real(max_icp_270), -imag(max_icp_270));
    
end

%%
%Check for rot invariance later...











%%
dt_pi10 = dtwavexfm2(step_edges_128(:,:,3), 4, 'near_sym_b','qshift_b');

w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; % Nominally j * pi/2, but reduced a bit due to asymmetry of subband freq responses.
p = [1 3]/4;

temp = cpxinterp2(dt_pi10{4}(:,:,1), p-0.5, w(1,:),'spline');
dt_pi10_4_1i = temp(1:size(dt_pi10{3},1),1:size(dt_pi10{3},2));

figure;
subplot(1,2,1); subimage(complex2rgb(dt_pi10{4}(:,:,1))); axis image;
title('Level 4 co-efficients: Magnitude and phase');
subplot(1,2,2); subimage(complex2rgb(dt_pi10_4_1i)); axis image;
title('Interpolated level 4 co-efficients: Magnitude and phase');
%%
figure;
for ang = 1:20
    ls = 'b';
    if any([3 6 9 13 16 19] == ang)
        ls = 'r';
    end
    compass(cos(angles(ang)), sin(angles(ang)), ls); hold on;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = -1 + i*1;
b = -1 + i*-1;

c1 = interp1([0 2], [a b], 1);
c2 = interp1([0 2], abs([a b]), 1)*exp(i*interp1([0 2], angle([a b]), 1));
c3 = 2*ifft([fft([a b]) 0 0]);
c3 = c3(2);

figure; 
compass(a, 'r'); hold on;
compass(b, 'b');
compass(c1, 'g');
compass(c2, 'm');
compass(c3, 'c');

axis equal


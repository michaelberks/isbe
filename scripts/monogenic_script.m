%... monogenic_script
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
    
    [local_amp_p, local_phase_p, local_ori_p] = monogenic(pos_edges_128(:,:,edge_angle), 2, 4, 2, 0.65);
    [local_amp_n, local_phase_n, local_ori_n] = monogenic(neg_edges_128(:,:,edge_angle), 2, 4, 2, 0.65);
    
    figure;
    %subplot(2,3,1); imagesc(pos_edges_128(:,:,edge_angle)); axis image; colormap(gray(256)); hold on;
    
    for level = 1:3
        subplot(2,3,level);
        image(complex2rgb(local_amp_p(:,:,level).*exp(i*local_ori_p(:,:,level)))); axis image;
        subplot(2,3,level+3);
        image(complex2rgb(local_amp_n(:,:,level).*exp(i*local_ori_n(:,:,level)))); axis image;     
    end
    %display(round(20*angle(max_icp_pos)/pi));
    %display(abs(max_ilp_pos));
    %display(abs(max_ilp_pos_opp));
end
%%
pos_lines_256 = ones(256); pos_lines_256(128,:) = 2;
neg_lines_256 = 2*ones(256); neg_lines_256(128,:) = 1;

pos_lines_128 = zeros(128,128,length(angles));
neg_lines_128 = zeros(128,128,length(angles));

for ii = 1:length(angles);
    
    temp1 = imrotate(pos_lines_256, 180*angles(ii)/pi, 'crop');
    temp2 = imrotate(neg_lines_256, 180*angles(ii)/pi, 'crop');
    pos_lines_128(:,:,ii) = temp1(65:192, 65:192);
    neg_lines_128(:,:,ii) = temp2(65:192, 65:192);
    clear temp*;
    
end
%%
for line_angle = 1:20 %[3 6 9 13 16 19] %lines aligned at sub-band orientations
    
    [local_amp_p, local_phase_p, local_ori_p] = monogenic(pos_lines_128(:,:,line_angle), 2, 4, 2, 0.65);
    [local_amp_n, local_phase_n, local_ori_n] = monogenic(neg_lines_128(:,:,line_angle), 2, 4, 2, 0.65);
    
    figure;
    for level = 1:3
        subplot(2,3,level);
        image(complex2rgb(local_amp_p(:,:,level).*exp(i*local_ori_p(:,:,level)))); axis image;
        subplot(2,3,level+3);
        image(complex2rgb(local_amp_n(:,:,level).*exp(i*local_ori_n(:,:,level)))); axis image;     
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Create pos line, neg line, pos edge, neg_edge
[local_amp_pl, local_phase_pl, local_ori_pl] = monogenic([zeros(128,63) ones(128,2) zeros(128,63)], 4, 4, 2, 0.65);
[local_amp_nl, local_phase_nl, local_ori_nl] = monogenic([ones(128,63) zeros(128,2) ones(128,63)], 4, 4, 2, 0.65);
[local_amp_pe, local_phase_pe, local_ori_pe] = monogenic([zeros(128,64) ones(128,64)], 4, 4, 2, 0.65);
[local_amp_ne, local_phase_ne, local_ori_ne] = monogenic([ones(128,64) zeros(128,64)], 4, 4, 2, 0.65);
%%
phase_c_pl = 0;
phase_c_nl = 0;
phase_c_pe = 0;
phase_c_ne = 0;

f1 = figure; hold all;
f2 = figure; hold all;
f3 = figure; hold all;
f4 = figure; hold all;
for level = 2:5
%     figure;
%     subplot(2,2,1); image(complex2rgb(local_amp_pl(:,:,level).*exp(i*local_ori2_pl(:,:,level)))); axis image;
%     subplot(2,2,2); image(complex2rgb(local_amp_nl(:,:,level).*exp(i*local_ori2_nl(:,:,level)))); axis image;
%     subplot(2,2,3); image(complex2rgb(local_amp_pe(:,:,level).*exp(i*local_ori2_pe(:,:,level)))); axis image;
%     subplot(2,2,4); image(complex2rgb(local_amp_ne(:,:,level).*exp(i*local_ori2_ne(:,:,level)))); axis image; 
    
    figure(f1); plot(1:128, local_ori_pl(64,:,level)); legend show
    figure(f2); plot(1:128, local_ori_nl(64,:,level)); legend show
    figure(f3); plot(1:128, local_ori_pe(64,:,level)); legend show
    figure(f4); plot(1:128, local_ori_ne(64,:,level)); legend show
    
    phase_c_pl = phase_c_pl + local_amp_pl(:,:,level).*exp(i*local_ori_pl(:,:,level));
    phase_c_nl = phase_c_nl + local_amp_nl(:,:,level).*exp(i*local_ori_nl(:,:,level));
    phase_c_pe = phase_c_pe + local_amp_pe(:,:,level).*exp(i*local_ori_pe(:,:,level));
    phase_c_ne = phase_c_ne + local_amp_ne(:,:,level).*exp(i*local_ori_ne(:,:,level));
    
    figure;
    subplot(2,2,1); imagesc(abs(local_amp_pl(:,:,level))); axis image;
    subplot(2,2,2); imagesc(abs(local_amp_nl(:,:,level))); axis image;
    subplot(2,2,3); imagesc(abs(local_amp_pe(:,:,level))); axis image;
    subplot(2,2,4); imagesc(abs(local_amp_ne(:,:,level))); axis image;
end
figure;
subplot(2,2,1); image(complex2rgb(phase_c_pl)); axis image;
subplot(2,2,2); image(complex2rgb(phase_c_nl)); axis image;
subplot(2,2,3); image(complex2rgb(phase_c_pe)); axis image;
subplot(2,2,4); image(complex2rgb(phase_c_ne)); axis image;

figure;
subplot(2,2,1); imagesc(abs(phase_c_pl)); axis image; colormap(gray(256));
subplot(2,2,2); imagesc(abs(phase_c_nl)); axis image; colormap(gray(256));
subplot(2,2,3); imagesc(abs(phase_c_pe)); axis image; colormap(gray(256));
subplot(2,2,4); imagesc(abs(phase_c_ne)); axis image; colormap(gray(256));

%%
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Create radial lines image
%--------------------------------------------------------------------------
angles = 0:18:162;
pos_line512 = zeros(512); pos_line512([17:224 289:496], 256:257) = 1;
radial_lines = zeros(512);
for ii = 1:length(angles)
    radial_lines = radial_lines + imrotate(pos_line512, angles(ii), 'crop');
end
radial_lines = 100*(radial_lines+1);
figure; image(radial_lines); axis image; colormap(gray(256));

dt_radial_lines = dtwavexfm2(radial_lines, 5, 'near_sym_b','qshift_b');
[ilp_radial_lines icp_radial_lines] = mb_dual_tree_transform(dt_radial_lines);
[max_icp_radial_lines] = max(icp_radial_lines{2}, [], 3);
figure; image(complex2rgb(max_icp_radial_lines)); axis image;
%%
[local_amp_rad, local_phase_rad, local_ori_rad] = monogenic(radial_lines, 4, 4, 2, 0.65);
for level = 1:4
    figure; 
    subplot(1,2,1); image(complex2rgb(local_amp_rad(:,:,level).*exp(i*local_ori_rad(:,:,level)))); axis image;
    subplot(1,2,2); image(complex2rgb(local_amp_rad(:,:,level).*exp(2*i*local_phase_rad(:,:,level)))); axis image;
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
test_image = load('M:\chen\data\testimage_contrast1to8_exprnd_sin\image001.mat');
[local_amp_test, local_phase_test, local_ori_test] = monogenic(test_image.image, 4, 4, 2, 0.65);
phase_cong_test = 0;
for level = 2:5
    figure; 
    subplot(1,2,1); image(complex2rgb(local_amp_test(:,:,level).*exp(i*local_ori_test(:,:,level)))); axis image;
    subplot(1,2,2); image(complex2rgb(local_amp_test(:,:,level).*exp(2*i*local_phase_test(:,:,level)))); axis image;
    phase_cong_test = phase_cong_test + local_amp_test(:,:,level).*exp(i*local_ori_test(:,:,level));
end
figure; image(complex2rgb(phase_cong_test)); axis image;

pos_lines = imag(phase_cong_test);
pos_lines(pos_lines < 0) = 0;
figure; imagesc(pos_lines); axis image; colormap(gray(256));

[pos_lines2, phase_cong] = monogenic_phase_cong(test_image.image, 4, 4, 2, 0.65, 'pos_line');
figure; imagesc(pos_lines2); axis image; colormap(gray(256));
%%
test_image = load('M:\chen\data\testimage_contrast1to8_exprnd_sin\image001.mat');
for min_wav = 1:4
    for onf = 0.65
        for lev = 2:6
            figure; imagesc(monogenic_phase_cong(test_image.image, lev, min_wav, 2, onf)); axis image; colormap(gray(256));
            title([num2str(lev) ' levels, Onf = ' num2str(onf), ', Min wave length = ', num2str(min_wav)]);
        end
    end
end



%%
    seed_image = imread('http://graphics.cs.cmu.edu/people/efros/research/NPS/images/fill-bread.gif');

 	sample_image = imread('http://graphics.cs.cmu.edu/people/efros/research/NPS/images/fill-bread_w11.gif');

	filled_image = seed_image ~= 0;

	syn_image7 = cjr_efros_tex_syn('SampleImage', sample_image,...
        'SeedImage', seed_image, 'FilledImage', filled_image, 'WindowSize', 7);
    
%%
    syn_image11 = cjr_efros_tex_syn('SampleImage', sample_image,...
        'SeedImage', seed_image, 'FilledImage', filled_image, 'WindowSize', 11);
    syn_image15 = cjr_efros_tex_syn('SampleImage', sample_image,...
        'SeedImage', seed_image, 'FilledImage', filled_image, 'WindowSize', 15);
    syn_image19 = cjr_efros_tex_syn('SampleImage', sample_image,...
        'SeedImage', seed_image, 'FilledImage', filled_image, 'WindowSize', 19);
    
    figure;
    subplot(2,2,1)
    title('Window = 7');
    imagesc(syn_image7); colormap gray; axis image;
    subplot(2,2,2)
    title('Window = 11');
    imagesc(syn_image11); colormap gray; axis image;
    subplot(2,2,3)
    title('Window = 15');
    imagesc(syn_image15); colormap gray; axis image;
    subplot(2,2,4)
    title('Window = 19');
    imagesc(syn_image19); colormap gray; axis image;
%%
    syn_image11 = cjr_efros_tex_syn('SampleImage', sample_image,...
        'SeedImage', seed_image, 'FilledImage', filled_image,...
        'WindowSize', 11, 'Debug', 0);
    figure, title('Window = 11'), imagesc(syn_image11); colormap gray; axis image;


%% Make outlines
    s_data(100).mass_outline = [];
    s_data(100).mass_ROI = [];
    
    f1 = figure;
    f2 = figure;
    
    for ii = 1:100
        a = 100 + 20*randn(1); b = 100 + 20*randn(1);
        
        [xx yy] = ellipse(a, b, 2, 2, [1 0], 1);
        s_data(ii).mass_outline = [xx', yy'];
        
        mass_ROI = 128 + 32*randn(1)+ 4*randn(600, 600);
        mass_bw = (a-b)*poly2mask(xx+300, yy+300, 600, 600);
        s_data(ii).mass_ROI = uint8(mass_ROI + mass_bw);
        s_data(ii).mass_centroid = [mean(xx) + 300, mean(yy)+300];
        
    end
%%

%     shape_vec = m_data(13).mass_outline...
%         + repmat(m_data(13).mass_centroid,...
%         size(m_data(13).mass_outline,1), 1);
%     
%     %size_shape_vec = length(shape_vec);
%     
%     mass_ROI = m_data(13).mass_ROI;
%     
%     [r c] = size(mass_ROI);
%     clear data_in;
%     
%     %compute BW mask and list of belonging to the shape
%     shape_bw = roipoly(mass_ROI, shape_vec(:,1), shape_vec(:,2));
%     shape_pl = regionprops(bwlabel(shape_bw, 4), 'PixelList');
%     tex_vec = shape_pl.PixelList;
%     
%     %dilate the shape mask by n1
%     mask1 = shape_bw;
%     for jj = 1:20
%         mask1 = imdilate(mask1, strel('disk', 1));
%     end
%     
%     %dilate the shape mask by n2
%     mask2 = shape_bw;
%     for jj = 1:80
%         mask2 = imdilate(mask2, strel('disk', 1));
%     end
%     
%     %subtract mask1 from mask2, result is ring in shape of mass, n1 pixels
%     % from the shape border, n2-n1 pixels wide
%     mask3 = mask2 - mask1; clear mask2 mask1;
%     
%     %create a mask of equally spaced dots (as defined by input: spacing) 
%     mask4 = zeros([r c]);
%     mask4(1:20:end, 1:20:end) = 1;
%     
%     %mask of landmark pts is intersection of spaced dots and ring of shape
%     mask5 = mask3 & mask4; clear mask3 mask4
%     [landmark_y landmark_x] = ind2sub(size(mask5), find(mask5));
%     
%     image(255*ones(size(shape_bw))); colormap gray; axis image; hold on;
%     plot(shape_vec(:,1), shape_vec(:,2));
%     plot(landmark_x, landmark_y, 'rx');

%     %aviobj = avifile('K:\isbe\project\misc\ss07\com_modes.avi','fps', 10);
%     for ii = 1:20
%         open(strcat('K:\isbe\project\misc\ss07\com_modes\', num2str(ii), '.fig'));
%         saveas(gcf, strcat('K:\isbe\project\misc\ss07\com_modes\', num2str(ii), '.jpg'));
%         %frame = getframe(gcf);
%         close gcf;
%         %aviobj = addframe(aviobj,frame);
%     end
%     for ii = 1:19
%         kk = 20 - ii;
%         open(strcat('K:\isbe\project\misc\ss07\com_modes\', num2str(ii), '.fig'));
%         saveas(gcf, strcat('K:\isbe\project\misc\ss07\com_modes\', num2str(ii), '.jpg'));
%         %frame = getframe(gcf);
%         close gcf;
%         %aviobj = addframe(aviobj,frame);
%     end
%     %aviobj = close(aviobj); clear aviobj;

% s_x = mass_model.mean_shape(1:200) + mass_model.mean_off_c;
% s_y = mass_model.mean_shape(201:end) + mass_model.mean_off_r;
% 
% %Define points to be interpolated by TPS - as row vectors
% i_x = mass_model.mean_shape_pl(:,1)';
% i_y = mass_model.mean_shape_pl(:,2)';
% 
% tps_L_inv = tps_weights(s_x, s_y);
% 
% %Define displacement to target points
% z_x = mass_outline(:,1)';
% z_y = mass_outline(:,2)';
% 
% %Compute displacement of interpolated points
% f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
% f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
% 
% %Get texture vector
% X_tex = interp2(double(tex_data(3).subtract_ROI), f_x, f_y);
% 
% temp = zeros(mass_model.mean_r, mass_model.mean_c);
% temp(sub2ind([mass_model.mean_r, mass_model.mean_c], i_y, i_x)) = uint8(X_tex);
% figure;
% imagesc(temp); colormap gray; axis image; clear temp;
% 
% clear s_x s_y z_x z_y i_x i_y f_x f_y tps_L_inv tex_data;
%
% % figure;
% % for ii = 1:5;
% %     [r c] = size(tex_data(ii).subtract_ROI);
% %     t_im1 = bg_im(end-r+1:end, end-c+1:end) + uint8(tex_data(ii).subtract_ROI);
% %     t_im2 = imfilter(uint8(tex_data(ii).subtract_ROI), gf, 'replicate');
% %     t_im2 = bg_im(end-r+1:end, end-c+1:end) + t_im2;
% %     subplot(3, 5, ii)
% %     image(uint8(m_data(ii).mass_ROI)); axis image; colormap(gray(256));
% %     subplot(3, 5, ii+5)
% %     image(t_im1); axis image; colormap(gray(256));
% %     subplot(3, 5, ii+10)
% %     image(t_im2); axis image; colormap(gray(256));
% %     clear t_im;
% % end
% % figure;
% % for ii = 6:10;
% %     [r c] = size(tex_data(ii).subtract_ROI);
% %     t_im1 = bg_im(end-r+1:end, end-c+1:end) + uint8(tex_data(ii).subtract_ROI);
% %     t_im2 = imfilter(uint8(tex_data(ii).subtract_ROI), gf, 'replicate');
% %     t_im2 = bg_im(end-r+1:end, end-c+1:end) + t_im2;
% %     subplot(3, 5, ii-5)
% %     image(uint8(m_data(ii).mass_ROI)); axis image; colormap(gray(256));
% %     subplot(3, 5, ii)
% %     image(t_im1); axis image; colormap(gray(256));
% %     subplot(3, 5, ii+5)
% %     image(t_im2); axis image; colormap(gray(256));
% %     clear t_im;
% % end
% % figure;
% % for ii = 11:15;
% %     [r c] = size(tex_data(ii).subtract_ROI);
% %     t_im1 = bg_im(end-r+1:end, end-c+1:end) + uint8(tex_data(ii).subtract_ROI);
% %     t_im2 = imfilter(uint8(tex_data(ii).subtract_ROI), gf, 'replicate');
% %     t_im2 = bg_im(end-r+1:end, end-c+1:end) + t_im2;
% %     subplot(3, 5, ii-10)
% %     image(uint8(m_data(ii).mass_ROI)); axis image; colormap(gray(256));
% %     subplot(3, 5, ii-5)
% %     image(t_im1); axis image; colormap(gray(256));
% %     subplot(3, 5, ii)
% %     image(t_im2); axis image; colormap(gray(256));
% %     clear t_im;
% % end
% % figure;
% % for ii = 16:20;
% %     [r c] = size(tex_data(ii).subtract_ROI);
% %     t_im1 = bg_im(end-r+1:end, end-c+1:end) + uint8(tex_data(ii).subtract_ROI);
% %     t_im2 = imfilter(uint8(tex_data(ii).subtract_ROI), gf, 'replicate');
% %     t_im2 = bg_im(end-r+1:end, end-c+1:end) + t_im2;
% %     subplot(3, 5, ii-15)
% %     image(uint8(m_data(ii).mass_ROI)); axis image; colormap(gray(256));
% %     subplot(3, 5, ii-10)
% %     image(t_im1); axis image; colormap(gray(256));
% %     subplot(3, 5, ii-5)
% %     image(t_im2); axis image; colormap(gray(256));
% %     clear t_im;
% % end
%
% % fid = fopen('C:\isbe\mammograms\transfer\annotations\anno_list.txt');
% % file_names = textscan(fid, '%s');
% % clear fid mass_files;
% % 
% % figure,
% % for ii = 1:size(file_names{1}, 1)
% %     mass = load(file_names{1}{ii});
% %     spicules = mass.mass_spicules;
% %     mass_ROI = mass.mass_ROI;
% %     mass_outline = mass.mass_outline;    
% %     clear mass;
% %     
% %     subplot(4,5,ii);
% %     imagesc(mass_ROI); colormap gray; axis image; hold on;
% %     plot([mass_outline(:,1); mass_outline(1,1)],...
% %          [mass_outline(:,2); mass_outline(1,2)], 'LineWidth', 1.25); 
% %     for jj = 1:length(spicules)
% %         s = spicules(jj).outline;
% %         plot(s(:,1), s(:,2), 'r', 'LineWidth', 1.25);
% %     end
% % end

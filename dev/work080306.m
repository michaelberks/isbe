% Take some pyramid sub-band with string structures in and see if structure
% really is conditioned from course to fine levels

real_line = double(rgb2gray(imread('C:\isbe\dev\background\images\normal_2\normal016.bmp')));
[pyr_real p_sizes] = buildSFpyr(double(real_line), 5, 4);
pyramid_real = mb_change_pyramid_form(pyr_real, p_sizes);
%%
for level = 2:5
    for ori = 1:5
        figure; imagesc(pyramid_real{level,ori}); axis image; colormap(jet(256));
    end
end
%%

% Ok, sub-band level = 5, ori = 4 has strong linear structure. So do
% texture synthesis on {4,4} to compare un/conditioned methods

%select subregion to synthesis (i.e. around major structure)
figure; imagesc(pyramid_real{6,4}); axis image; colormap(jet(256));
filled_image = roipoly;

%%
clear syn_args
%args for both methods
syn_args.PixelOrder = mb_compute_pixel_order(~filled_image);
syn_args.Uint8Mode = 0;
syn_args.ForceCluster = 0;
syn_args.TargetImage = pyramid_real{6, 3};

% %args for unconditioned synthesis
% syn_args.PathToTextureGMM = ...
%             'C:/isbe/dev/background/results/normal_2/normal_2_model_4_4.mat';
% 
% [tex_un cluster_image_un] = ...
%     mb_gmm_tex_synthesis(syn_args);

%args for conditioned synthesis
syn_args.SampleImage2 = pyramid_real{7, 1};
syn_args.PathToTextureGMM = ...
            'C:/isbe/dev/background/results/mass_2_2levels/mass_2_20_model_6_3.mat';

[tex_con cluster_image_con] = ...
    mb_gmm_tex_synthesis2(syn_args);

figure; imagesc(pyramid_real{6,3}); axis image; colormap(jet(256));
figure; imagesc(tex_con); axis image; colormap(jet(256));
%%
figure; imagesc(tex_un); axis image; colormap(jet(256));
figure; imagesc(tex_con); axis image; colormap(jet(256));
%%
figure; imagesc(cluster_image_un); axis image; colormap(jet(256));
figure; imagesc(cluster_image_con); axis image; colormap(jet(256));
%%
figure; imagesc(real_line); axis image; colormap(gray(256));
filled_image = roipoly;
%%
pyr_args.TargetImage = real_line;
pyr_args.FilledImage = ~filled_image; %#ok
pyr_args.ModelName = 'mass_2_model';
pyr_args.ModelDir = [mberksroot, 'background/results/mass_2_2levels/'];
pyr_args.CutOffLevel = 5;
pyr_args.ConditionLevels = 1;

[synthesised_image pyramid cluster_image] = mb_gmm_pyr_synthesis(pyr_args);
%%
for ii = [3 4 6 8 12 15 16 20 24 34 36 40 51 52 53]
    i1 = imread(['C:\isbe\dev\background\images\mass_2\mass', zerostr(ii,3), '.bmp']);
    figure('Name', ['mass', zerostr(ii,3)]); imagesc(i1); colormap(gray(256)); axis image;
    clear i1;
end
%%
for ii = [3 4 6 8 12 15 16 20]% 24 34 36 40 51 52 53]
    i1 = imread(['C:\isbe\dev\background\images\mass_2\mass', zerostr(ii,3), '.bmp']);
    i1 = imfilter(i1, fspecial('gaussian', 25, 5), 'symmetric');
    [response orientations] = line_operator_conv(i1, 12, 11, 15);
    nms = non_maximal_supp(response, orientations);
%     [rows cols] = find(nms > 0.4);
%     
%     figure('Name', ['mass', zerostr(ii,3)]); imagesc(i1); colormap(gray(256)); axis image;
%     hold on;
%     plot(cols, rows, 'r.', 'MarkerSize', 2);
    
    figure; imagesc(nms); axis image; colormap(gray(256));
    clear i1 response orientations response;
end
%%
for ii = [3 4 6 8 12 15 16 20]% 24 34 36 40 51 52 53]
    i1 = imread(['C:\isbe\dev\background\images\mass_2\mass', zerostr(ii,3), '.bmp']);
    figure('Name', ['mass', zerostr(ii,3)]); imagesc(i1); colormap(gray(256)); axis image;
    
    [response] = gabor_filter(i1, 4, 8, 12, 4)
    figure; imagesc(response); axis image; colormap(gray(256));
    clear i1 response;
end
%%
for ii = 1:9
    load([mberksroot, 'background/syn/mass_2_2levels/conditioned_synthesis_circle', zerostr(ii, 3)]);
%     [p pp] = mb_change_pyramid_form(pyramid);
%     synthesised_image = reconSFpyr(p, pp);
%     save([mberksroot, 'background/syn/mass_2_2levels/conditioned_synthesis_circle', zerostr(ii, 3)],...
%         'synthesised_image', 'pyramid', 'cluster_image');
    i1 = double(imread(['C:\isbe\dev\background\images\normal_2\normal', zerostr(ii, 3), '.bmp']));
    figure; 
    subplot(1,2,1); image(uint8(synthesised_image)); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(i1-synthesised_image); axis image; colormap(gray(256));
    clear synthesised_image pyramid cluster_image

end

%%
for ii = 1:9
    load([mberksroot, 'background/syn/mass_2_2levels/conditioned_synthesis_circle', zerostr(ii, 3)]);
    figure; 
    imagesc(uint8(synthesised_image)); axis image; colormap(gray(256));
    clear synthesised_image pyramid cluster_image

end
%%
for ii = 1:9
    load([mberksroot, 'background/syn/mass_2_2levels/conditioned_synthesis_circle', zerostr(ii, 3)]);
    i1 = double(imread(['C:\isbe\dev\background\images\normal_2\normal', zerostr(ii, 3), '.bmp']));
    figure;
    subplot(1,2,1); image(uint8(synthesised_image)); axis image; colormap(gray(256));
    subplot(1,2,2); image(uint8(i1)); axis image; colormap(gray(256));
    clear synthesised_image pyramid cluster_image

end
%%
for ii = 1:15
    pyr = u_load([mberksroot,...
        'background/pyramid/normal_2/normal', zerostr(ii, 3), '.bmp_pyramid.mat']);
    figure; 
    subplot(1,2,1); imagesc(pyr{7,1}); colormap(gray(256)); axis image;
    subplot(1,2,2); imagesc(pyr{6,1}); colormap(gray(256)); axis image;
    clear pyr;
end

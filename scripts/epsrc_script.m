% Script to produce images and interesting results for the EPSRC meeting
% This contains a good summary of work done on pyramid modelling/synthesis

% Pyramid decompostions of different classes of regions:

%%
%1. Normal region
pyramid = u_load([mberksroot, 'background/pyramid/normal512/o04_104LML_1024_4493_1688_pyramid.mat']);
%%
for lev = 2:6
    for ori = 1:4
        %figure; imagesc(pyramid{lev,ori}); axis image; caxis([-4 4]);
        write_im_from_colormap(pyramid{lev,ori},...
            ['C:\isbe\epsrc_meeting\figures\o04_104LML', num2str(lev), '_', num2str(ori), '.bmp'],...
            jet(256), [-4 4]);
    end
end
%%
[pyr_vec pyr_sizes] = mb_change_pyramid_form(pyramid);
recon_image123456 = mb_reconSFpyr(pyr_vec, pyr_sizes, (1:6)');
recon_image23456 = mb_reconSFpyr(pyr_vec, pyr_sizes, (2:6)');
recon_image3456 = mb_reconSFpyr(pyr_vec, pyr_sizes, (3:6)');
recon_image456 = mb_reconSFpyr(pyr_vec, pyr_sizes, (4:6)');
recon_image56 = mb_reconSFpyr(pyr_vec, pyr_sizes, (5:6)');
recon_image6 = mb_reconSFpyr(pyr_vec, pyr_sizes, (6:6)');
recon_image012345 = mb_reconSFpyr(pyr_vec, pyr_sizes, (0:5)');
recon_image01234 = mb_reconSFpyr(pyr_vec, pyr_sizes, (0:4)');
recon_image0123 = mb_reconSFpyr(pyr_vec, pyr_sizes, (0:3)');
recon_image012 = mb_reconSFpyr(pyr_vec, pyr_sizes, (0:2)');
recon_image01 = mb_reconSFpyr(pyr_vec, pyr_sizes, (0:1)');
recon_image1 = mb_reconSFpyr(pyr_vec, pyr_sizes, 1);
recon_image2 = mb_reconSFpyr(pyr_vec, pyr_sizes, 2);
recon_image3 = mb_reconSFpyr(pyr_vec, pyr_sizes, 3);
recon_image4 = mb_reconSFpyr(pyr_vec, pyr_sizes, 4);
recon_image5 = mb_reconSFpyr(pyr_vec, pyr_sizes, 5);
write_im_from_colormap(recon123456, 'C:\isbe\epsrc_meeting\figures\n_recon123456.bmp');
write_im_from_colormap(recon_image123456, 'C:\isbe\epsrc_meeting\figures\n_recon123456.bmp');
write_im_from_colormap(recon_image23456, 'C:\isbe\epsrc_meeting\figures\n_recon23456.bmp');
write_im_from_colormap(recon_image3456, 'C:\isbe\epsrc_meeting\figures\n_recon3456.bmp');
write_im_from_colormap(recon_image456, 'C:\isbe\epsrc_meeting\figures\n_recon456.bmp');
write_im_from_colormap(recon_image56, 'C:\isbe\epsrc_meeting\figures\n_recon56.bmp');
write_im_from_colormap(recon_image6, 'C:\isbe\epsrc_meeting\figures\n_recon6.bmp');
write_im_from_colormap(recon_image012345, 'C:\isbe\epsrc_meeting\figures\n_recon012345.bmp');
write_im_from_colormap(recon_image01234, 'C:\isbe\epsrc_meeting\figures\n_recon01234.bmp');
write_im_from_colormap(recon_image0123, 'C:\isbe\epsrc_meeting\figures\n_recon0123.bmp');
write_im_from_colormap(recon_image012, 'C:\isbe\epsrc_meeting\figures\n_recon012.bmp');
write_im_from_colormap(recon_image01, 'C:\isbe\epsrc_meeting\figures\n_recon01.bmp');
write_im_from_colormap(recon_image1, 'C:\isbe\epsrc_meeting\figures\n_recon1.bmp');
write_im_from_colormap(recon_image2, 'C:\isbe\epsrc_meeting\figures\n_recon2.bmp');
write_im_from_colormap(recon_image3, 'C:\isbe\epsrc_meeting\figures\n_recon3.bmp');
write_im_from_colormap(recon_image4, 'C:\isbe\epsrc_meeting\figures\n_recon4.bmp');
write_im_from_colormap(recon_image5, 'C:\isbe\epsrc_meeting\figures\n_recon5.bmp');

%%
%2. Subtracted mass region

%%
%3. Mass region
load C:\isbe\dev\files\bg_files.mat
mass_006 = u_load(['C:\isbe\dev\masses\', bg_files(6).name]);
mass_006 = imresize(double(mass_006.mass_ROI), 0.5);

[pyr_vec006, pyr_bands006] = mb_buildSFpyr(mass_006, 5, 3);

pyr006 =  mb_change_pyramid_form(pyr_vec006, pyr_bands006);
recon_image1 = mb_reconSFpyr(pyr_vec006, pyr_bands006, 1);
recon_image2 = mb_reconSFpyr(pyr_vec006, pyr_bands006, 2);
recon_image3 = mb_reconSFpyr(pyr_vec006, pyr_bands006, 3);
recon_image4 = mb_reconSFpyr(pyr_vec006, pyr_bands006, 4);
recon_image5 = mb_reconSFpyr(pyr_vec006, pyr_bands006, 5);
figure; imagesc(recon_image1); axis image; colormap(gray(256));
figure; imagesc(recon_image2); axis image; colormap(gray(256));
figure; imagesc(recon_image3); axis image; colormap(gray(256));
figure; imagesc(recon_image4); axis image; colormap(gray(256));
figure; imagesc(recon_image5); axis image; colormap(gray(256));
write_im_from_colormap(recon_mass006_1, 'C:\isbe\dev\background\misc\recon\recon1.bmp');
write_im_from_colormap(recon_mass006_2, 'C:\isbe\dev\background\misc\recon\recon2.bmp');
write_im_from_colormap(recon_mass006_3, 'C:\isbe\dev\background\misc\recon\recon3.bmp');
write_im_from_colormap(recon_mass006_4, 'C:\isbe\dev\background\misc\recon\recon4.bmp');
write_im_from_colormap(recon_mass006_5, 'C:\isbe\dev\background\misc\recon\recon5.bmp');

for lev = 2:6
    for ori = 1:4
        figure; imagesc(pyr006{lev, ori}); axis image; caxis([-4 4]);
    end
end

%%
%4. Micro-calcs
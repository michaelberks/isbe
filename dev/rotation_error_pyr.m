function rot_errors = rot_pyramid_errors(region_main)
%%

if ~isa(region_main, 'double')
    region_main = double(region_main);
end

figure; imagesc(region_main); axis image; colormap(gray(256));

[p p_sizes] = mb_buildSFpyr(region_main, 5, 3);
pyr_main = mb_change_pyramid_form(p, p_sizes); clear p;

%%
if 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute errors for 90 deg rotation
region_90 = imrotate(region_main, 90, 'nearest', 'crop');
[p p_sizes] = mb_buildSFpyr(region_90, 5, 3);
pyr_90 = mb_change_pyramid_form(p, p_sizes);

pyr_recon = pyr_main;
%%
rot_errors.r90.subbands = zeros(5, 4);

for level = 2:6
    lim(1) = 2^(10-level) + 1;
    lim(2) = 3*2^(10-level);
    
    x = 0;
    if level > 2
        x = 1;
    end

    for ori = 1:4

        ori_2_4 = rem(ori+1, 4) + 1;
        test_region_main = pyr_main{level,ori_2_4}(lim(1):lim(2),(lim(1):lim(2))+x);
        test_region90 = imrotate(pyr_90{level,ori}(lim(1):lim(2),(lim(1):lim(2))), -90, 'crop');
        
        if ori > 2
            test_region90 = -test_region90;
        end
        
        rot_errors.r90.subbands(level-1, ori) = ...
            mean(abs(test_region_main(:) - test_region90(:))) / mean(abs(test_region_main(:)));
        
        pyr_recon{level,ori_2_4}((lim(1):lim(2)),(lim(1):lim(2))+x) = ...
            test_region90;

%         figure;
%         subplot(2,2,1); imagesc(test_region_main(1:end,1:end)); axis image; a1 = gca; caxis([-2 2]); colorbar; 
%         subplot(2,2,2); imagesc(test_region90(1:end,1:end)); axis image; a2 = gca; caxis([-2 2]); colorbar;
%         subplot(2,2,3); imagesc(abs(test_region90(1:end,1:end) - test_region_main(1:end,1:end))); axis image; a3 = gca; caxis([-2 2]); colorbar;
%         linkaxes([a1 a2 a3]);
    end
end

p = mb_change_pyramid_form(pyr_recon); clear pyr_recon;
recon_90 = mb_reconSFpyr(p, p_sizes); clear p;
rot_errors.r90.recon = ...
    region_main(257:768, 257:768) - recon_90(257:768, 257:768);
figure; imagesc(recon_90 - region_main); axis image;
clear *90*
%%
end

%%
if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute errors for 90 deg rotation
region_180 = imrotate(region_main, 180, 'nearest', 'crop');
[p p_sizes] = mb_buildSFpyr(region_180, 5, 3);
pyr_180 = mb_change_pyramid_form(p, p_sizes); clear p;

rot_errors.r180.subbands = zeros(5, 4);

pyr_recon = pyr_main;
%%
for level = 2:6
    lim(1) = 2^(10-level) + 1;
    lim(2) = 3*2^(10-level);
    
    x = 0; y = 0;
    if level > 2
        x = 1;
        y = 1;
    end
    
    for ori = 1:4

        test_region_main = pyr_main{level,ori}((lim(1):lim(2))+y,(lim(1):lim(2))+x);
        test_region180 = -imrotate(pyr_180{level,ori}(lim(1):lim(2),(lim(1):lim(2))), -180, 'crop');
        
        
        rot_errors.r180.subbands(level-1, ori) = ...
            mean(abs(test_region_main(:) - test_region180(:))) / mean(abs(test_region_main(:)));
        
        pyr_recon{level,ori}((lim(1):lim(2))+y,(lim(1):lim(2))+x) = ...
            test_region180;

%         figure;
%         subplot(2,2,1); imagesc(test_region_main(1:end,1:end)); axis image; a1 = gca; caxis([-2 2]); colorbar; 
%         subplot(2,2,2); imagesc(test_region180(1:end,1:end)); axis image; a2 = gca; caxis([-2 2]); colorbar;
%         subplot(2,2,3); imagesc(abs(test_region180(1:end,1:end) - test_region_main(1:end,1:end))); axis image; a3 = gca; caxis([-2 2]); colorbar;
%         linkaxes([a1 a2 a3]);
    end
end

p = mb_change_pyramid_form(pyr_recon); clear pyr_recon;
recon_180 = mb_reconSFpyr(p, p_sizes); clear p;
rot_errors.r180.recon = ...
    region_main(257:768, 257:768) - recon_180(257:768, 257:768);
figure; imagesc(recon_180 - region_main); axis image; colormap(jet(256));
clear *180*
%%

end
%%
if 1
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute errors for 90 deg rotation
region_270 = imrotate(region_main, 270, 'nearest', 'crop');
[p p_sizes] = mb_buildSFpyr(region_270, 5, 3);
pyr_270 = mb_change_pyramid_form(p, p_sizes); clear p;

pyr_recon = pyr_main;
%%
rot_errors.r270.subbands = zeros(5, 4);

for level = 2:6
    lim(1) = 2^(10-level) + 1;
    lim(2) = 3*2^(10-level);
    
    x = 0;
    if level > 2
        x = 1;
    end

    for ori = 1:4

        ori_2_4 = rem(ori+1, 4) + 1;
        test_region_main = pyr_main{level,ori_2_4}((lim(1):lim(2)) + x,(lim(1):lim(2)));
        test_region270 = imrotate(pyr_270{level,ori}(lim(1):lim(2),(lim(1):lim(2))), -270, 'crop');
        
        if ori < 3
            test_region270 = -test_region270;
        end
        
        rot_errors.r270.subbands(level-1, ori) = ...
            mean(abs(test_region_main(:) - test_region270(:))) / mean(abs(test_region_main(:)));
        
        pyr_recon{level,ori_2_4}((lim(1):lim(2)) + x,(lim(1):lim(2))) = ...
            test_region270;

%         figure;
%         subplot(2,2,1); imagesc(test_region_main(1:end,1:end)); axis image; a1 = gca; caxis([-2 2]); colorbar; 
%         subplot(2,2,2); imagesc(test_region270(1:end,1:end)); axis image; a2 = gca; caxis([-2 2]); colorbar;
%         subplot(2,2,3); imagesc(abs(test_region270(1:end,1:end) - test_region_main(1:end,1:end))); axis image; a3 = gca; caxis([-2 2]); colorbar;
%         linkaxes([a1 a2 a3]);
    end
end

p = mb_change_pyramid_form(pyr_recon); clear pyr_recon;
recon_270 = mb_reconSFpyr(p, p_sizes); clear p;
rot_errors.r270.recon = ...
   region_main(257:768, 257:768) - recon_270(257:768, 257:768);
figure; imagesc(recon_270 - region_main); axis image; colormap(jet(256));
clear *270*
%%
end
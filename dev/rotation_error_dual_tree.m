function rot_errors = rotation_error_dual_tree(region_main)
%%
%
% MB: This is designed to work with patches of a set size 512x512 - through
% laziness at the moment. Should be modified to test regions of any size!

if ~isa(region_main, 'double')
    region_main = double(region_main);
end

figure; imagesc(region_main); axis image; colormap(gray(256));

[Yh_main] = dtwavexfm2(region_main, 5, 'near_sym_b','qshift_b');

%%
if 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute errors for 90 deg rotation
region_90 = imrotate(region_main, 90, 'nearest', 'crop');
[Yh_90] = dtwavexfm2(region_90, 5, 'near_sym_b','qshift_b');

Yh_recon = Yh_main;
%%
rot_errors.r90.subbands = zeros(5, 6);

% Limits for each level are
% 1: 256: 65 - 192
% 2: 128: 33 - 96
% 3: 64:  17 - 48
% 4: 32:  9  - 24
% 5: 16:  5  - 12
for level = 1:5
    lim(1) = 2^(7-level) + 1;
    lim(2) = 3*2^(7-level);
    

    for ori = 1:6

        ori_2_4 = rem(ori+2, 6) + 1;
        test_region_main = Yh_main{level}(lim(1):lim(2),(lim(1):lim(2)),ori_2_4);
        test_region90 = imrotate(Yh_90{level}(lim(1):lim(2),(lim(1):lim(2)),ori), -90, 'crop');
        
        
        if ori < 4
            test_region90 = imag(test_region90) + i*real(test_region90);
        else
            test_region90 = -imag(test_region90) + i*real(test_region90);
        end
        
        rot_errors.r90.subbands(level, ori) = ...
            mean(abs(test_region_main(:) - test_region90(:))) / mean(abs(test_region_main(:)));
        
        Yh_recon{level}((lim(1):lim(2)),(lim(1):lim(2)), ori_2_4) = ...
            test_region90;
%         if level > 1
%             figure;
%             subplot(2,2,1); imagesc(abs(test_region_main(1:end,1:end))); axis image; a1 = gca; colorbar; 
%             subplot(2,2,2); imagesc(abs(test_region90(1:end,1:end))); axis image; a2 = gca; colorbar;
%             subplot(2,2,3); imagesc(angle(test_region_main(1:end,1:end))); axis image; a3 = gca; colorbar;
%             subplot(2,2,4); imagesc(angle(test_region90(1:end,1:end))); axis image; a4 = gca; colorbar;
%             linkaxes([a1 a2 a3 a4]);
%         end
    end
end

recon_90 = dtwaveifm2(Yh_recon,'near_sym_b','qshift_b');
rot_errors.r90.recon = ...
    region_main(129:384, 129:384) - recon_90(129:384, 129:384);
figure; imagesc(recon_90 - region_main); axis image;
clear *90*
%%
end

%%
if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute errors for 180 deg rotation
region_180 = imrotate(region_main, 180, 'nearest', 'crop');
[Yh_180] = dtwavexfm2(region_180, 5, 'near_sym_b','qshift_b');

Yh_recon = Yh_main;
%%
rot_errors.r180.subbands = zeros(5, 6);

for level = 1:5
    lim(1) = 2^(7-level) + 1;
    lim(2) = 3*2^(7-level);
    

    for ori = 1:6

        test_region_main = Yh_main{level}(lim(1):lim(2),(lim(1):lim(2)),ori);
        test_region180 = imrotate(Yh_180{level}(lim(1):lim(2),(lim(1):lim(2)),ori), -180, 'crop');
        
        
        if ori < 4
            test_region180 = -real(test_region180) + i*imag(test_region180);
        else
            test_region180 = real(test_region180) - i*imag(test_region180);
        end
        
        rot_errors.r180.subbands(level, ori) = ...
            mean(abs(test_region_main(:) - test_region180(:))) / mean(abs(test_region_main(:)));
        
        Yh_recon{level}((lim(1):lim(2)),(lim(1):lim(2)), ori) = ...
            test_region180;
%         if level > 1
%             figure;
%             subplot(2,2,1); imagesc(abs(test_region_main(1:end,1:end))); axis image; a1 = gca; colorbar; 
%             subplot(2,2,2); imagesc(abs(test_region180(1:end,1:end))); axis image; a2 = gca; colorbar;
%             subplot(2,2,3); imagesc(angle(test_region_main(1:end,1:end))); axis image; a3 = gca; colorbar;
%             subplot(2,2,4); imagesc(angle(test_region180(1:end,1:end))); axis image; a4 = gca; colorbar;
%             linkaxes([a1 a2 a3 a4]);
%         end
    end
end

recon_180 = dtwaveifm2(Yh_recon,'near_sym_b','qshift_b');
rot_errors.r180.recon = ...
    region_main(129:384, 129:384) - recon_180(129:384, 129:384);
figure; imagesc(recon_180 - region_main); axis image; colormap(jet(256));
clear *180*
%%

end
%%
if 1
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute errors for 270 deg rotation
region_270 = imrotate(region_main, 270, 'nearest', 'crop');
[Yh_270] = dtwavexfm2(region_270, 5, 'near_sym_b','qshift_b');

Yh_recon = Yh_main;
%%
rot_errors.r270.subbands = zeros(5, 6);

for level = 1:5
    lim(1) = 2^(7-level) + 1;
    lim(2) = 3*2^(7-level);

    for ori = 1:6

        ori_2_4 = rem(ori+2, 6) + 1;
        test_region_main = Yh_main{level}(lim(1):lim(2),(lim(1):lim(2)),ori_2_4);
        test_region270 = imrotate(Yh_270{level}(lim(1):lim(2),(lim(1):lim(2)),ori), -270, 'crop');
        
        
        if ori < 4
            test_region270 = imag(test_region270) - i*real(test_region270);
        else
            test_region270 = imag(test_region270) + i*real(test_region270);
        end
        
        rot_errors.r270.subbands(level, ori) = ...
            mean(abs(test_region_main(:) - test_region270(:))) / mean(abs(test_region_main(:)));
        
        Yh_recon{level}((lim(1):lim(2)),(lim(1):lim(2)), ori_2_4) = ...
            test_region270;
%         if level > 1
%             figure;
%             subplot(2,2,1); imagesc(abs(test_region_main(1:end,1:end))); axis image; a1 = gca; colorbar; 
%             subplot(2,2,2); imagesc(abs(test_region270(1:end,1:end))); axis image; a2 = gca; colorbar;
%             subplot(2,2,3); imagesc(angle(test_region_main(1:end,1:end))); axis image; a3 = gca; colorbar;
%             subplot(2,2,4); imagesc(angle(test_region270(1:end,1:end))); axis image; a4 = gca; colorbar;
%             linkaxes([a1 a2 a3 a4]);
%         end
    end
end

recon_270 = dtwaveifm2(Yh_recon,'near_sym_b','qshift_b');
rot_errors.r270.recon = ...
   region_main(129:384, 129:384) - recon_270(129:384, 129:384);
figure; imagesc(recon_270 - region_main); axis image; colormap(jet(256));
clear *270*
%%
end
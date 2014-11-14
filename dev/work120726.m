ret = rgb2gray( u_load('C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\images\01_test.mat') );

f_mask = u_load('C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\foveal_masks\01_test_f_mask.mat');
v_mask = u_load('C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\vessel_masks\01_test_v_mask.mat');

%%
[r c] = size(ret);
responses = zeros(r,c,4);
orientations = zeros(r,c,4);
for ss = 1:4
    [responses(:,:,ss) orientations(:,:,ss)] = gaussian_2nd_derivative_line(ret, 2^(ss-1));
end

[~, scales_idx] = max(abs(responses), [], 3);
scales = 2.^(scales_idx-1);
    
%%
f1 = figure; 
subplot(1,2,1); imgray(ret); aor = gca;
subplot(1,2,2); imgray(v_mask); aom = gca;
linkaxes([aor aom]);
%%
scale_cols = 'rgby';
for ii = 1:3
    figure(f1);
    [cx cy ret_g] = improfile('bicubic');
    [ret_s] = improfile(scales, cx, cy, 'nearest');
    [ret_m] = improfile(v_mask, cx, cy, 'nearest');
    plot(aor, cx, cy, 'r');

    %
    figure; 
    subplot(3,1,1); hold all; title('Intensity profile'); ag = gca;
    subplot(3,1,2); hold all; title('G" response profile'); ar = gca;
    subplot(3,1,3); hold all; title('G" scale profile'); as = gca;
    
    for ss = 1:4
        [ret_r] = improfile(responses(:,:,ss), cx, cy, 'bicubic');
        [ret_o] = improfile(orientations(:,:,ss), cx, cy, 'nearest');
        plot(ar, ret_r);
        dx1 = cx + 5*cos(ret_o);
        dx2 = cx - 5*cos(ret_o);
        dy1 = cy - 5*sin(ret_o);
        dy2 = cy + 5*sin(ret_o);
        plot(aom, [dx1 dx2]', [dy1 dy2]', scale_cols(ss));
    end
    legend(ar, {'\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8'});
    plot(ar, ret_m*8, 'k--');
    plot(ag, ret_g); plot(ag, ret_m*8 + mean(ret_g), 'k--');  
    plot(as, ret_s); plot(as, ret_m*8, 'k--');
end

 
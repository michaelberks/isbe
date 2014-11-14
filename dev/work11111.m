test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\';
prob_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\';
results_dir = 'C:\isbe\asymmetry_project\experiments\line_detection\';
param_dir = {'44412', '44413', '44414', '44415', '44477'};

for ii = 1:5
    [orientation_errors] = compute_image_orientation_errors(...
         test_dir, [prob_dir param_dir{ii} '\'], 'centre_line');
end
%%
for ii = 1:5
    load([prob_dir param_dir{ii} '\centre_lineerrors\ori_errors.mat']);
    display([mean(abs(orientation_errors)) median(abs(orientation_errors))]);
end
%%
for ii = 1:5
    [orientation_errors] = compute_image_orientation_errors(...
         test_dir, [prob_dir param_dir{ii} '\'], 'all_line');
     display(param_dir{ii});
     display([mean(abs(orientation_errors)) median(abs(orientation_errors))]);
end
%%
[orientation_errors] = compute_image_orientation_errors(...
         test_dir, ['Z:\asym\data\synthetic_lines\real512\results\clover_analytic\orientations\'], 'all_line');
display('clover');
display([mean(abs(orientation_errors)) median(abs(orientation_errors))]);
%%
test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\';
for s = 4%[1 2 4 8]
    for w = 3:5%2:5
        clover_dir = ['C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clovert_w' num2str(w) '_s' num2str(s)];
        mkdir([clover_dir '\orientations\']);
        mkdir([clover_dir '\lines\']);
        for ii = 1:100
            load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\image' zerostr(ii,3) '.mat']);
            [line_map, ori_map] = gaussian_clover_line_exp(test_image, s, s*w);
%             save_uint8([clover_dir '\orientations\image' zerostr(ii,3) '_class.mat'], ori_map);
%             save_uint8([clover_dir '\lines\image' zerostr(ii,3) '_class.mat'], line_map);
            save([clover_dir '\orientations\image' zerostr(ii,3) '_class.mat'], 'ori_map');
        end
        [orientation_errors] = compute_image_orientation_errors(...
             test_dir, [clover_dir '\orientations\'], 'centre_line');
    display(clover_dir);
    display([mean(abs(orientation_errors)) median(abs(orientation_errors))]);
    end
end
%%
figure; hold all;
for s = 4%[1 2 4 8]
    for w = 2:5
        load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w' num2str(w) '_s' num2str(s) ...
            '\orientations\all_lineerrors\ori_errors.mat']);
        
        orientation_errors = sort(abs(orientation_errors));
        ori_cdf = zeros(101,1);
        for jj = 1:100
            x = ceil(jj*size(orientation_errors, 1)/100);
            ori_cdf(jj+1) = orientation_errors(x,1);
        end
        plot(ori_cdf, (0:100)/100, 'linewidth', 2);
    end
end
%%
for w = 2:5
    clover_dir = ['C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w' num2str(w) '\scales\'];
    [scales] =...
        get_image_classifications(test_dir, clover_dir, 'all_line');
    scales = round(scales);
    display(clover_dir);
    display([sum(scales == 1) sum(scales == 2) sum(scales == 4) sum(scales == 8)]);
end
%%
s = 4;
w = 3;
clover_dir = ['C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_pad_s' num2str(s)];
mkdir([clover_dir '\orientations\']);
mkdir([clover_dir '\lines\']);
for ii = 1:100
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\image' zerostr(ii,3) '.mat']);
    [line_map, ori_map] = gaussian_clover_line(test_image, s, w);
    save_uint8([clover_dir '\orientations\image' zerostr(ii,3) '_class.mat'], ori_map);
    save_uint8([clover_dir '\lines\image' zerostr(ii,3) '_class.mat'], line_map);
end
[orientation_errors] = compute_image_orientation_errors(...
     test_dir, [clover_dir '\orientations\'], 'all_line');
display(clover_dir);
display([mean(abs(orientation_errors)) median(abs(orientation_errors))]);
%%
sigma = 4;
figure;
for w = [5 4 3 2]
    width = round(w*sigma);
    sigmasq = sigma^2;

    x	= (-width:width);
    g	= exp(-0.5* (x.*x)/sigmasq) / sqrt(2*pi);
    dg	= -x/sigmasq .* g;
    ddg	= (-1/sigmasq * g) - (x/sigmasq .* dg);
    
    subplot(1,3,1); hold all; plot(x, g, 'linewidth', 2);
    subplot(1,3,2); hold all; plot(x, dg, 'linewidth', 2);
    subplot(1,3,3); hold all; plot(x, ddg, 'linewidth', 2);
end
legend({'26.2', '24.7', '18.8', '19.2'});
%%
% [ori_errors3] = compute_image_orientation_errors(...
%      'C:\isbe\asymmetry_project\data\synthetic_lines\real512\',...
%      'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\orientations\', 'centre_line');
% [ori_errors4] = compute_image_orientation_errors(...
%      'C:\isbe\asymmetry_project\data\synthetic_lines\real512\',...
%      'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w4_s4\orientations\', 'centre_line');
% [ori_errors5] = compute_image_orientation_errors(...
%      'C:\isbe\asymmetry_project\data\synthetic_lines\real512\',...
%      'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w5_s4\orientations\', 'centre_line');
 
ori_errors3 = ...
     u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\orientations\centre_lineerrors\ori_errors.mat');
ori_errors4 = ...
     u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w4_s4\orientations\centre_lineerrors\ori_errors.mat');
ori_errors5 = ...
     u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w5_s4\orientations\centre_lineerrors\ori_errors.mat');
ori_errors3t = ...
     u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clovert_w3_s4\orientations\centre_lineerrors\ori_errors.mat');
ori_errors4t = ...
     u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clovert_w4_s4\orientations\centre_lineerrors\ori_errors.mat');
ori_errors5t = ...
     u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clovert_w5_s4\orientations\centre_lineerrors\ori_errors.mat');
 
% [ori_errors_rf] = compute_image_orientation_errors(...
%      'C:\isbe\asymmetry_project\data\synthetic_lines\real512\',...
%      'Z:\asym\data\synthetic_lines\real512\results\286712\', 'centre_line');
% [ori_errors_g2] = compute_image_orientation_errors(...
%      'C:\isbe\asymmetry_project\data\synthetic_lines\real512\',...
%      'Z:\asym\data\synthetic_lines\real512\results\10762\', 'centre_line');
 
% figure; plot(mod(ori3,pi), ori_errors3, 'rx');
% figure; plot(mod(ori4,pi), ori_errors4, 'rx');
% figure; plot(mod(true_oris, 180), ori_errors3, 'rx');
% figure; plot(mod(true_oris, 180), ori_errors_rf, 'rx');
% figure; plot(mod(true_oris, 180), ori_errors_g2, 'rx');

err_dist3 = hist3([mod(true_oris, 180) ori_errors3], [180 180]);
err_dist5 = hist3([mod(true_oris, 180) ori_errors5], [180 180]);
err_dist3t = hist3([mod(true_oris, 180) ori_errors3t], [180 180]);
err_dist5t = hist3([mod(true_oris, 180) ori_errors5t], [180 180]);

% figure;
% subplot(2,1,1); imagesc(err_dist3');
% subplot(2,1,2); imagesc(err_dist5');

figure;
subplot(2,1,1); imagesc(log(err_dist3' ./ repmat(sum(err_dist3'),180,1)));
subplot(2,1,2); imagesc(log(err_dist5' ./ repmat(sum(err_dist5'),180,1)));

figure;
subplot(2,1,1); imagesc(log(err_dist3t' ./ repmat(sum(err_dist3t'),180,1)));
subplot(2,1,2); imagesc(log(err_dist5t' ./ repmat(sum(err_dist5t'),180,1)));

err_dist3l = log(err_dist3' ./ repmat(sum(err_dist3'),180,1)); err_dist3l(err_dist3l == -inf) = min(err_dist3l(err_dist3l ~= -inf));
err_dist5l = log(err_dist5' ./ repmat(sum(err_dist5'),180,1)); err_dist5l(err_dist5l == -inf) = min(err_dist5l(err_dist5l ~= -inf));
err_dist3tl = log(err_dist3t' ./ repmat(sum(err_dist3t'),180,1)); err_dist3tl(err_dist3tl == -inf) = min(err_dist3tl(err_dist3tl ~= -inf));
err_dist5tl = log(err_dist5t' ./ repmat(sum(err_dist5t'),180,1)); err_dist5tl(err_dist5tl == -inf) = min(err_dist5tl(err_dist5tl ~= -inf));

write_im_from_colormap(err_dist3l, 'C:\isbe\asymmetry_project\data\misc\gauss_err_dist_log3.bmp', jet(256));
write_im_from_colormap(err_dist5l, 'C:\isbe\asymmetry_project\data\misc\gauss_err_dist_log5.bmp', jet(256));
write_im_from_colormap(err_dist3tl, 'C:\isbe\asymmetry_project\data\misc\gauss_err_dist_log3t.bmp', jet(256));
write_im_from_colormap(err_dist5tl, 'C:\isbe\asymmetry_project\data\misc\gauss_err_dist_log5t.bmp', jet(256));

%%
[ori3] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\orientations\', 'centre_line');
[ori4] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w4_s4\orientations\', 'centre_line');
[ori_rf] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'Z:\asym\data\synthetic_lines\real512\results\286712\', 'centre_line');
ori3c = exp(2*i*ori3);
ori4c = exp(2*i*ori4);

mean_ori = angle((ori3c + ori4c)/2) / 2;
ori_diff = mb_mod(ori3 - ori4, pi);

figure; plot(mean_ori, ori_diff, 'rx');
figure; hist(mod(mean_ori,pi), 180);
figure; hist(mod(ori3,pi), 180);
figure; hist(mod(ori4,pi), 180);
figure; hist(mod(true_oris,180), 180);
%%
[theta_a3] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\orientations\theta_a\', 'all_line');
[theta_a4] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w4_s4\orientations\theta_a\', 'all_line');
theta_a3c = exp(2*i*theta_a3);
theta_a4c = exp(2*i*theta_a4);

mean_theta = angle((theta_a3c + theta_a4c)/2) / 2;
theta_diff = mb_mod(theta_a3 - theta_a4, pi);

figure; plot(mean_theta, theta_diff, 'r.');
%%
[theta_b3] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\orientations\theta_b\', 'all_line');
[theta_b4] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w4_s4\orientations\theta_b\', 'all_line');
theta_b3c = exp(2*i*theta_b3);
theta_b4c = exp(2*i*theta_b4);

mean_theta = angle((theta_b3c + theta_b4c)/2) / 2;
theta_diff = mb_mod(theta_b3 - theta_b4, pi);

figure; plot(mean_theta, theta_diff, 'rx');
%%
test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\';
for s = 4%[1 2 4 8]
    for w = 3%3:4
        clover_dir = ['C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w' num2str(w) '_s' num2str(s)];
        mkdir([clover_dir '\orientations\theta_a\']);
        mkdir([clover_dir '\orientations\theta_b\']);
        mkdir([clover_dir '\lines\response_a\']);
        mkdir([clover_dir '\lines\response_b\']);
        for ii = 1:100
            load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\image' zerostr(ii,3) '.mat']);
            [dummy dummy dummy theta_a theta_b response_a response_b] = gaussian_clover_line_exp(test_image, s, s*w);
            save([clover_dir '\orientations\theta_a\' zerostr(ii,3) '_class.mat'], 'theta_a');
            save([clover_dir '\orientations\theta_b\' zerostr(ii,3) '_class.mat'], 'theta_b');
            save([clover_dir '\lines\response_a\' zerostr(ii,3) '_class.mat'], 'response_a');
            save([clover_dir '\lines\response_b\' zerostr(ii,3) '_class.mat'], 'response_b');
        end

    end
end

%%
true_oris = [];
for ii = 1:100
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\labels\label' zerostr(ii,3) '.mat']);
    label = label_centre & label < 2;
    true_oris = [true_oris; label_orientation(label)];
end
true_oris = mod(true_oris, 180);
save C:\isbe\asymmetry_project\data\misc\true_oris.mat true_oris
            
    
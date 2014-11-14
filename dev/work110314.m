results_dir_list = dir('Z:\data\synthetic_lines\real512\results\');

image_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512';

for ii = 3:length(results_dir_list)
    results_dir = ['Z:\data\synthetic_lines\real512\results\' results_dir_list(ii).name '\'];
    error_dir = [results_dir 'errors\'];
    try
        [ori_errors] = compute_image_orientation_errors(...
            image_dir,...
            results_dir);
    catch
        display(['Couldn''t compute errors for ' results_dir_list(ii).name]);
        continue;
    end
    display(['Errors computed for ' results_dir_list(ii).name]);
    mkdir(error_dir);
    save([error_dir 'ori_errors.mat'], 'ori_errors');
end
%%
[ori_errors] = compute_image_orientation_errors(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512', ...
    'Z:\data\synthetic_lines\real512\results\310941\');
mkdir Z:\data\synthetic_lines\real512\results\310941\errors\
save Z:\data\synthetic_lines\real512\results\310941\errors\ori_errors.mat
%%
px_per_mm = 100/9;
r_max = round([9.6 12.0 16.2 18.6 23.4] * px_per_mm);

load C:\isbe\asymmetry_project\data\misc\ori_maps.mat ori_map_3_4_old
line_map_full = abs(ori_map_3_4_old);
ori_map_full = angle(ori_map_3_4_old);

line_map_half = abs(imresize(ori_map_3_4_old, 0.5, 'bilinear'));
ori_map_half = angle(imresize(ori_map_3_4_old, 0.5, 'bilinear'));
%
profile on;
[f_i1_full f_i2_full] = karssemeijer_radial_projection_multiscale(line_map_full, ori_map_full, 10, r_max, 5, 24, 10, 1);
profsave(profile('info'),'C:\isbe\asymmetry_project\data\misc\profiles\k_full')

profile on;
[f_i1_full2 f_i2_full2] = karssemeijer_radial_projection_multiscale(line_map_full, ori_map_full, 10, r_max, 5, 24, 10, 2);
profsave(profile('info'),'C:\isbe\asymmetry_project\data\misc\profiles\k_full2')

profile on;
[f_i1_half f_i2_half] = karssemeijer_radial_projection_multiscale(line_map_half, ori_map_half, 10, round(r_max/2), 5, 24, 10, 1);
profsave(profile('info'),'C:\isbe\asymmetry_project\data\misc\profiles\k_half')
%%
forest_list = {...
    '243311';...
    '243629';...
    '243630';...
    '286712';...
    '286713';...
    '286714';...
    '286715';...
    '286716';
    '287361';
    '290903';...
    '299400';...
    '299402';...
    '299407';...
    '299430';...
    '299442';...
    '305891';...
    '305894';...
    '305897';...
    '306060';...
    '306076';...
    '306079'};
num_forests = length(forest_list);

mean_ae = zeros(num_forests,1);
median_ae = zeros(num_forests,1);
mean_awae = zeros(num_forests,1);

for ii = 1:num_forests
    try
        errors = u_load(['Z:\data\synthetic_lines\real512\results\' forest_list{ii} '\errors\ori_errors.mat']);
    catch
        errors = u_load(['Z:\data\synthetic_lines\real512\results\' forest_list{ii} '\200\errors\ori_errors.mat']);
    end
        
    ori_errors = sort(abs(errors(:,1)));
    mag_errors = sort(errors(:,2));
    com_errors = sort(prod(abs(errors),2) / mean(mag_errors));
    errors = sortrows(errors, -2);
    
    mean_ae(ii) = mean(ori_errors);
    median_ae(ii) = median(ori_errors);
    mean_awae(ii) = mean(com_errors);
end
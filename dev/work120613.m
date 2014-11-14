for ii = 21:40
    v_mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\' zerostr(ii,2) '_training_v_mask.mat']);
    f_mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\foveal_masks\' zerostr(ii,2) '_training_f_mask.mat']);
    vessel_map = v_mask;
    bg_map = f_mask & ~v_mask;
    
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\orig\vessel_maps\' zerostr(ii,2) '_training_v_map.mat'], 'vessel_map');
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\orig\bg_maps\' zerostr(ii,2) '_training_b_map.mat'], 'bg_map');
    
    figure; 
    subplot(1,2,1); imgray(vessel_map);
    subplot(1,2,2); imgray(bg_map);
end
%%
args.num_samples = 2e4;
args.image_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\';
args.bg_sample_map_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\orig\bg_maps\';
args.vessel_sample_map_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\orig\vessel_maps\';
args.prediction_type = 'detection';
args.debug = 1;
args.decomp_type = 'pixel';
args.win_size = 1;

resample_vessel_data(args);
%%
mkdir C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it1\vessel_maps
mkdir C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it1\bg_maps

for ii = 21:40
    v_mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\' zerostr(ii,2) '_training_v_mask.mat']);
    f_mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\foveal_masks\' zerostr(ii,2) '_training_f_mask.mat']);
    v_prob = load_uint8(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\results\107714\' zerostr(ii,2) '_training_class.mat']);
    
    vessel_map = 1 - v_prob;
    vessel_map(~v_mask) = 0;
    
    bg_map = v_prob;
    bg_map(~f_mask | v_mask) = 0;
    
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it1\vessel_maps\' zerostr(ii,2) '_training_v_map.mat'], 'vessel_map');
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it1\bg_maps\' zerostr(ii,2) '_training_b_map.mat'], 'bg_map');
end
%%
mkdir C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it2\vessel_maps
mkdir C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it2\bg_maps

for ii = 21:40
    v_mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\' zerostr(ii,2) '_training_v_mask.mat']);
    f_mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\foveal_masks\' zerostr(ii,2) '_training_f_mask.mat']);
    v_prob = load_uint8(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\results\107714\' zerostr(ii,2) '_training_class.mat']);
    
    vessel_map = 1 - v_prob + 0.5;
    vessel_map(~v_mask) = 0;
    
    bg_map = v_prob + 0.5;
    bg_map(~f_mask | v_mask) = 0;
    
    figure; 
    subplot(1,2,1); imgray(vessel_map);
    subplot(1,2,2); imgray(bg_map);
    
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it2\vessel_maps\' zerostr(ii,2) '_training_v_map.mat'], 'vessel_map');
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it2\bg_maps\' zerostr(ii,2) '_training_b_map.mat'], 'bg_map');
end
%%
args.num_samples = 2e4;
args.image_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\';
args.bg_sample_map_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it1\bg_maps\';
args.vessel_sample_map_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it1\vessel_maps\';
args.prediction_type = 'detection';
args.debug = 1;
args.decomp_type = 'pixel';
args.win_size = 1;

[training_data training_labels] = resample_vessel_data(args);
%%
d_root = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\';
job = {'108229', '108476', '108025'};

figure; axis([0 1 0 1]); axis equal; hold all;
title('So does the whole resampling work???');

for jj = 1:3
    
    tp_counts = zeros(20, 102);
    fp_counts = zeros(20, 102);
    t_counts = zeros(20, 1);
    f_counts = zeros(20, 1);

    for ii = 1:20

        vessel_prob = ...
            load_uint8([d_root 'images\results\' job{jj} '\' zerostr(ii,2) '_test_class.mat']);

        v_mask = u_load([d_root 'vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
        f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);

        %Compute ROC counts for image
        [~, ~, tp_count fp_count] = calculate_roc_curve(vessel_prob(f_mask), v_mask(f_mask),(-1:100)/100);

        %Increment total counts
        tp_counts(ii,:) = tp_count;
        fp_counts(ii,:) = fp_count;
        t_counts(ii) = sum(v_mask(f_mask));
        f_counts(ii) = sum(~v_mask(f_mask));

    end

    %Compute ROC points for complete set of data
    roc_pts = [sum(fp_counts)' / sum(f_counts) sum(tp_counts)' / sum(t_counts)];

    %Compute AUC for ROC curve
    auc(jj) = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
            0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );

     
    plot(roc_pts(:,1), roc_pts(:,2)); 
    plot(roc_pts(:,1), roc_pts(:,2), 'r.');
    plot(roc_pts(51,1), roc_pts(51,2), 'gx'); 
end
legend(...
    {[ '100 trees, uniform sampling of original images, AUC: ' num2str(auc(1))];...
    [ '100 trees, resampled original images, k = 0, AUC: ' num2str(auc(2))];...
    [ '100 trees, resampled original images, k = 0.5, AUC: ' num2str(auc(3))]});
%%  
tp_counts = zeros(20, 102);
fp_counts = zeros(20, 102);
t_counts = zeros(20, 1);
f_counts = zeros(20, 1);

for ii = 1:20

    vessel_prob1 = ...
        load_uint8([d_root 'images\results\107714\' zerostr(ii,2) '_test_class.mat']);
    vessel_prob2 = ...
        load_uint8([d_root 'images\results\108025\' zerostr(ii,2) '_test_class.mat']);
    vessel_prob = (vessel_prob1 + vessel_prob2) / 2;
    
    v_mask = u_load([d_root 'vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
    f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
    
    vessel_prob1(~f_mask) = 0;
    vessel_prob2(~f_mask) = 0;
    vessel_prob(~f_mask) = 0;
    
%     figure; 
%     axes('units', 'normalized', 'position', [0 .5 .5 .5]); imgray(vessel_prob1); axis off;
%     axes('units', 'normalized', 'position', [.5 .5 .5 .5]); imgray(vessel_prob2); axis off;
%     axes('units', 'normalized', 'position', [0 0 .5 .5]); imgray(vessel_prob); axis off;
%     axes('units', 'normalized', 'position', [.5 0 .5 .5]); imgray(v_mask); axis off;

    

    %Compute ROC counts for image
    [~, ~, tp_count fp_count] = calculate_roc_curve(vessel_prob(f_mask), v_mask(f_mask),(-1:100)/100);

    %Increment total counts
    tp_counts(ii,:) = tp_count;
    fp_counts(ii,:) = fp_count;
    t_counts(ii) = sum(v_mask(f_mask));
    f_counts(ii) = sum(~v_mask(f_mask));

end

%Compute ROC points for complete set of data
roc_pts = [sum(fp_counts)' / sum(f_counts) sum(tp_counts)' / sum(t_counts)];

%Compute AUC for ROC curve
auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
        0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );

figure; axis([0 1 0 1]); axis equal; hold all;
title('So does the whole resampling work???');
plot(roc_pts(:,1), roc_pts(:,2)); 
plot(roc_pts(:,1), roc_pts(:,2), 'r.');
plot(roc_pts(51,1), roc_pts(51,2), 'gx'); 

legend([ '200 trees, uniform and resampling, original images, AUC: ' num2str(auc)]);
%%
for ii = 1:20
    vessel_prob = ...
        load_uint8([d_root 'images\results\107880\' zerostr(ii,2) '_test_class.mat']);
    v_mask = u_load([d_root 'vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
    figure; 
    subplot(1,2,1); imgray(vessel_prob);
    subplot(1,2,2); imgray(v_mask);
end
%%
args.num_samples = 2e4;
args.image_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\';
args.bg_sample_map_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\orig\bg_maps\';
args.vessel_sample_map_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\orig\vessel_maps\';
args.prediction_type = 'detection';
args.debug = 0;
args.decomp_type = 'dt';
args.win_size = 1;

[training_data training_labels] = resample_vessel_data(args);
%%
str_fmt = [];
for ii = 1:48
    str_fmt = [str_fmt '%7.3f '];
end
str_fmt(end) = [];
str_fmt = [str_fmt '\n'];
fid1 = fopen('test_data_file.txt','wt');
fprintf(fid1, '%s \n', ['Sampled data:  ' datestr(now)]);
fprintf(fid1,str_fmt, training_data);
fclose(fid1);
%%
args.num_samples = 2e4;
args.image_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\';
args.sample_map_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\orig\';
args.vessel_mask_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\';
args.prediction_type = 'detection';
args.debug = 0;
args.decomp_type = 'pixel';
args.win_size = 1;

[a b] = resample_retinogram_data(args);
resample_retinogram_data(args);
%%
for ii = 21:40
    v_mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\' zerostr(ii,2) '_training_v_mask.mat']);
    f_mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\foveal_masks\' zerostr(ii,2) '_training_f_mask.mat']);
    v_prob = load_uint8(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\results\108229\' zerostr(ii,2) '_training_class.mat']);
    
    sample_map = 0.5 + abs(v_prob - v_mask);
    sample_map(~f_mask) = 0;
    
    figure; 
    subplot(1,2,1); imgray(v_mask);
    subplot(1,2,2); imgray(sample_map);
    
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it1\' zerostr(ii,2) '_training_s_map.mat'], 'sample_map');
end
%%
args.num_samples = 1e5;
args.image_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\';
args.sample_map_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\resampling_maps\it1\';
args.vessel_mask_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\';
args.prediction_type = 'detection';
args.debug = 0;
args.decomp_type = 'pixel';
args.win_size = 1;

[~, labels] = resample_retinogram_data(args);

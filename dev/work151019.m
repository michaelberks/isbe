%base_dir = 'A:\PROCAS_ALL_VDM\';
%case_list = dir([base_dir 'PROCAS_ALL_*]');

base_dir = 'A:\D_2.2_Datasets\125_375_Single_MLO_Dataset_VOLPARA\CONTROLS\';
case_list = dir([base_dir 'ASSURE_CONTROLS*']);

%%
for i_case = 1:10
    im_name = dir([base_dir case_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    
    if length(im_name) >= 1
        mam_vdm = imread([base_dir case_list(i_case).name '\' im_name(end).name]);
        
        %if ~isempty(strfind(im_name(end).name, 'RML'))
        %    mam_vdm = fliplr(mam_vdm);
        %end
        
        mask = mam_vdm == mam_vdm(1,1);
        mask = ~bwselect(mask, 1, 1);
        mask_inner = imerode(mask, strel('disk', 32));
        dual_tree = compute_dual_tree(mam_vdm, 6);

        figure; imgray(mam_vdm);
        %figure; imgray(mask_inner);
        %
        for i_lev = 5:6%1:6; 

            mask_i = imresize(mask_inner, size(dual_tree{i_lev}(:,:,1)));
            figure;
            for i_band = 1:6; 
                band = dual_tree{i_lev}(:,:,i_band);
                band = complex(real(band), abs(imag(band)));
                glims = [min(abs(band(mask_i))) max(abs(band(mask_i)))];
                subplot(2,3,i_band); imgray(abs(band)); caxis(glims);
                %subplot(2,3,i_band); imgray(complex2rgb(band, [-pi, pi], glims(2)));
            end
        end
    end
end
%%
num_cases = min(1000,length(case_list));
im_sizes = zeros(num_cases,2);
for i_case = 1:num_cases
    im_name = dir([base_dir case_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    
    if length(im_name) >= 1
        mam_vdm = imread([base_dir case_list(i_case).name '\' im_name(end).name]);
        im_sizes(i_case,:) = size(mam_vdm);
    end
end
%%
l6_mean_07 = [];
l6_mean_abs_07 = [];
l6_mean_10 = [];
l6_mean_abs_10 = [];
num_07 = 0;
num_10 = 0;
for i_case = 1:num_cases
    im_name = dir(['A:\PROCAS_ALL_VDM\' case_list(i_case).name '\*LCC*hint_densityMap*.pgm']);
    
    if length(im_name) >= 1
        mam_vdm = imread(['A:\PROCAS_ALL_VDM\' case_list(i_case).name '\' im_name(end).name]);
        dual_tree = compute_dual_tree(mam_vdm, 6);
        
        for i_lev = 1:6
            lev6 = dual_tree{i_lev};
            fold_idx = imag(lev6) < 0;
            lev6(fold_idx) = conj(lev6(fold_idx));
            dual_tree{i_lev} = lev6;
        end
        
        im_size = size(mam_vdm,1);        
        if (im_size < 1000)
            num_07 = num_07 + 1;
            if isempty(l6_mean_07)
                l6_mean_07 = dual_tree;
                l6_mean_abs_07 = cell(6,1);
                for i_lev = 1:6
                    l6_mean_abs_07{i_lev} = abs(dual_tree{i_lev});
                end
            else
                for i_lev = 1:6
                    l6_mean_07{i_lev} = l6_mean_07{i_lev} + dual_tree{i_lev};
                    l6_mean_abs_07{i_lev} = l6_mean_abs_07{i_lev} + abs(dual_tree{i_lev});
                end
            end
        else
            num_10 = num_10 + 1;
            if isempty(l6_mean_10)
                l6_mean_10 = dual_tree;
                l6_mean_abs_10 = cell(6,1);
                for i_lev = 1:6
                    l6_mean_abs_10{i_lev} = abs(dual_tree{i_lev});
                end
            else
                for i_lev = 1:6
                    l6_mean_10{i_lev} = l6_mean_10{i_lev} + dual_tree{i_lev};
                    l6_mean_abs_10{i_lev} = l6_mean_abs_10{i_lev} + abs(dual_tree{i_lev});
                end
            end
        end 
    end
end
%%
for i_lev = 1:6
    figure;
    for i_band = 1:6;
        subplot(2,3,i_band); imgray(complex2rgb(l6_mean_10{i_lev}(:,:,i_band), [0 pi])); 
    end
end
%%
for i_lev = 1:6
    figure;
    for i_band = 1:6;
        subplot(2,3,i_band); imgray(l6_mean_abs_10{i_lev}(:,:,i_band)); 
    end
end
%%
for i_lev = 1:6
    figure;
    for i_band = 1:6;
        band = l6_mean_10{i_lev}(:,:,i_band) ./ l6_mean_abs_10{i_lev}(:,:,i_band);
        subplot(2,3,i_band); imgray(complex2rgb(band, [0 pi])); 
    end
end
%%
all_responses = zeros(1e5, 6);
xy = zeros(1e5, 2);
curr_sample = 0;
for i_case = 1:length(case_list)
    im_name = dir([base_dir case_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    
    if length(im_name) >= 1
        mam_vdm = imread([base_dir case_list(i_case).name '\' im_name(end).name]);
        dual_tree = compute_dual_tree(mam_vdm, 6);
        [r c] = size(dual_tree{6}(:,:,1));
        rows = repmat((1:r)', 1, c);
        cols = repmat(1:c, r, 1);

        num_samples = r*c;
        row_idx = curr_sample+(1:num_samples);
        curr_sample = curr_sample + num_samples;
        for i_band = 1:6; 
            band = dual_tree{i_lev}(:,:,i_band);
            band = complex(real(band), abs(imag(band)));
            all_responses(row_idx,i_band) = band(:);
            xy(row_idx,1) = cols(:);
            xy(row_idx,2) = rows(:);
        end
    end
end
all_responses(curr_sample:end,:) = [];
xy(curr_sample:end,:) = [];
xyr = xy + rand(size(xy)) - 0.5;
%%
train_X = abs(all_responses);
no_dims = 2;
initial_dims = 6;
perplexity = 30;
theta = 0.5;

%Run t-sne
mappedX = fast_tsne(train_X, no_dims, initial_dims, perplexity, theta);
figure; hold all;
plot(mappedX(:,1), mappedX(:,2), 'r.');
%%
figure;
imgray(zeros(100));
plot(mappedX(:,1)+50, mappedX(:,2)+50, 'r.');
%%
%DO MANUALLY
mask = false(100,100,0);
mask(:,:,end+1) = roipoly;
%%
num_masks = size(mask,3);
mask(:,:,num_masks+1) = 0;
for i_m = 1:num_masks
    mask(:,:,num_masks+1) = mask(:,:,num_masks+1) | mask(:,:,i_m);
end
mask(:,:,num_masks+1) = ~mask(:,:,num_masks+1);
num_masks = num_masks + 1;
%%
mapped_idx = sub2ind([100 100], round(mappedX(:,2))+50, round(mappedX(:,1))+50);
figure; axis equal; hold all;
for i_m = 1:num_masks
    mask_i = mask(:,:,i_m);
    plot(mappedX(mask_i(mapped_idx),1), mappedX(mask_i(mapped_idx),2), '.');
end
%%
figure; axis ij equal; hold all;
for i_m = 1:num_masks-1
    mask_i = mask(:,:,i_m);
    plot(xyr(mask_i(mapped_idx),1), xyr(mask_i(mapped_idx),2), '.');
end
%%
figure;
for i_m = 1:num_masks
    mask_i = mask(:,:,i_m);
    subplot(3,4,i_m); axis ij equal; hold all; 
    plot(xyr(mask_i(mapped_idx),1), xyr(mask_i(mapped_idx),2), '.');
end

%%
opts = statset('MaxIter', 1000);
num_k = 10;
train_X = abs(all_responses);
[k_idx, k_centres, sumd] = kmeans(train_X, num_k, 'EmptyAction', 'drop', 'Replicates', 10, 'Options', opts);
%%
figure; axis equal; hold all;
for i_k = 1:num_k    
    plot(mappedX(k_idx == i_k,1), mappedX(k_idx == i_k,2), '.');
end
%%
figure; axis ij equal; hold all;
for i_k = 1:num_k
    if i_k == 6 || i_k == 8
        plot(0,0,'.');
    else
        plot(xyr(k_idx == i_k,1), xyr(k_idx == i_k,2), '.');
    end
end
%%
figure; 
for i_k = 1:num_k
    subplot(2,5,i_k); axis ij equal; hold all;
    plot(xyr(k_idx == i_k,1), xyr(k_idx == i_k,2), '.');   
end
%%
opts = statset('MaxIter', 1000);
num_k = 10;
train_X = [abs(all_responses) angle(all_responses)];
[k_idx, k_centres, sumd] = kmeans(train_X, num_k, 'EmptyAction', 'drop', 'Replicates', 10, 'Options', opts);
%%
figure; axis equal; hold all;
for i_k = 1:num_k    
    plot(mappedX(k_idx == i_k,1), mappedX(k_idx == i_k,2), '.');
end

figure; 
for i_k = 1:num_k
    subplot(3,4,i_k); axis ij equal; hold all;
    plot(xyr(k_idx == i_k,1), xyr(k_idx == i_k,2), '.');   
end
%%
decomposition_args.decomp_type = 'dt';      %Use the dual-tree
decomposition_args.win_size = 1;            %Window size about pixel to sample features
decomposition_args.levels = 6;            %Levels to include
decomposition_args.feature_shape = 'rect';  %Keep this as rect, although there are other things you can play with (see sample_dt_data)
decomposition_args.feature_type = 'conj';   %Again, 'conj' works best, but can use 'mag' for magnitude only, or 'phase' (see convert_complex_representation)
decomposition_args.do_max = 0;              %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0;              %Try and make rotation invariant (best left switched off)
decomposition_args.use_nag = 0;             %Use the NAG toolbox to interpolate (best left swicthed off)
decomposition_args.normalise = 0;

colors = lines(num_k);
for i_case = 1:20
    im_name = dir([base_dir case_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    mam_vdm = imread([base_dir case_list(i_case).name '\' im_name(end).name]);
    mask = mam_vdm == mam_vdm(1,1);
    mask = ~bwselect(mask, 1, 1);
    dual_tree = compute_dual_tree(mam_vdm, 6);

    [rows cols] = find(mask);
    [responses] = compute_filter_responses(mam_vdm, decomposition_args);
    sampled_features = sample_image_features(responses, rows, cols, decomposition_args);
    %
    k_dists = zeros(size(sampled_features,1), num_k);
    for i_k = 1:num_k
        k_dists(:,i_k) = sum(bsxfun(@minus, sampled_features, k_centres(i_k,:)).^2,2);
    end
    [~,assigned_k] = min(k_dists,[],2);
    %
    figure; 
    subplot(1,2,1); imgray(mam_vdm);
    subplot(1,2,2); axis equal ij; hold all;
    for i_k = 1:num_k
        idx = assigned_k == i_k;
        plot(cols(idx), rows(idx), '.', 'markeredgecolor', colors(i_k,:));
    end
    set(gca, 'xlim', [0 size(mam_vdm,2)], 'ylim', [0 size(mam_vdm,1)]);
end


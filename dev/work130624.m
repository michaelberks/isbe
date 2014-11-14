image_dir = 'C:\isbe\nailfold\data\rsa_study\test\images\';
vessel_im = u_load([image_dir '10147c.mat']);
vessel_patch = vessel_im(200:500, 300:500);
figure; imgray(vessel_im); caxis([-20 10]);
figure; imgray(vessel_patch);
%%

rf = u_load('C:\isbe\nailfold\models\vessel\orientation\rf_regression\182263\predictor.mat');
rf.tree_root = 'C:\isbe\nailfold\models\vessel\orientation\rf_regression/';
load('C:\isbe\nailfold\models\vessel\orientation\rf_regression\182263\job_args.mat');

[patch_ori] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', job_args.decomposition_args,...
    'predictor', rf, ...
    'prediction_type', 'rf_regression',...
    'output_type', 'orientation',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
figure; imgray(complex2rgb(patch_ori));
%%
[patch_ori90] = predict_image(...
    'image_in', rot90(vessel_patch),...
    'decomposition_args', job_args.decomposition_args,...
    'predictor', rf, ...
    'prediction_type', 'rf_regression',...
    'output_type', 'orientation',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
figure; imgray(complex2rgb(-rot90(patch_ori90,-1)));
%%
sc_dir = 'C:\isbe\nailfold\data\rsa_study\test\shape_contexts\';
sc_list = dir([sc_dir '*.mat']);
%shape_contexts_ri = [];

shape_contexts_a = [];

n_images = length(sc_list);
n_sampled_so_far = 0;
n_samples = 1e4;
shape_contexts_mp = zeros(n_samples, 120);


for i_sc = 1:n_images
    
    N = n_samples - n_sampled_so_far;
    p = 1 / (n_images - i_sc + 1);
    samples_per_image = sample_from_binomial(N, p, 1);
    display(['Sampling ' num2str(samples_per_image) ' vectors from image ' num2str(i_sc)]);
    
    load([sc_dir sc_list(i_sc).name], 'sc_pi');
    
    total_samples = size(sc_pi, 3);
    samples_per_image = min(samples_per_image, total_samples);
    rand_idx = randperm(total_samples);
    rand_idx(samples_per_image + 1:end) = [];
    
    sc_pi = sc_pi(:,:,rand_idx);
    sc_pi = permute(sc_pi, [3 1 2]);
    
    %pc_pi = permute(pc_pi, [3 1 2]);
    %shape_contexts_ri = [shape_contexts_ri; real(sc_pi(:,:)) imag(sc_pi(:,:))]; %#ok
    shape_contexts_mp(n_sampled_so_far+(1:samples_per_image),:) = [abs(sc_pi(:,:)) angle(sc_pi(:,:))];
    %shape_contexts_a = [shape_contexts_a; pc_pi(:,:)]; %#ok
    
    n_sampled_so_far = n_sampled_so_far + samples_per_image;
end
save('C:\isbe\nailfold\data\rsa_study\models\k_means_samples.mat', 'shape_contexts_mp');
%%
clear; pack;
load('C:\isbe\nailfold\data\rsa_study\models\k_means_samples.mat');

[cluster_idx cluster_centres] = kmeans(shape_contexts_mp, 20, 'options', statset('MaxIter', 1e4), 'replicates', 10, 'emptyaction', 'drop');
for ii = 1:20
    if ~isnan(cluster_centres(ii,1))
        shape_c = reshape(cluster_centres(ii,1:60), 5, 12) .* exp(1i*reshape(cluster_centres(ii,61:120), 5, 12));
        figure; imgray(display_shape_context(shape_c, 10));
    end
end
save('C:\isbe\nailfold\data\rsa_study\models\k_means_mp.mat', 'cluster_centres', 'cluster_idx');

%%
clear; pack;
load('C:\isbe\nailfold\data\rsa_study\models\k_means_samples.mat');
shape_contexts_c = shape_contexts_mp(:,1:60) .* exp(1i * shape_contexts_mp(:,61:120));
clear shape_contexts_mp;
shape_contexts_ri = [real(shape_contexts_c) imag(shape_contexts_c)];
clear shape_contexts_c;

num_clusters = 14;
[cluster_idx cluster_centres] = kmeans(shape_contexts_ri, num_clusters, 'options', statset('MaxIter', 1e4), 'replicates', 10, 'emptyaction', 'drop');
for ii = 1:num_clusters
    if ~isnan(cluster_centres(ii,1))
        shape_c = complex(reshape(cluster_centres(ii,1:60), 5, 12), reshape(cluster_centres(ii,61:120), 5, 12));
        figure; 
        subplot(1,2,1); imgray(display_shape_context(shape_c, 10, 'p'));
        subplot(1,2,2); imgray(display_shape_context(shape_c, 10, 'm'));
    end
end
save('C:\isbe\nailfold\data\rsa_study\models\k_means_ri.mat', 'cluster_centres', 'cluster_idx');

%%

sc_dir = 'C:\isbe\nailfold\data\rsa_study\test\shape_contexts\';
vc_dir = 'C:\isbe\nailfold\data\rsa_study\test\vessel_centres\';
image_dir = 'C:\isbe\nailfold\data\rsa_study\test\images\';
prob_dir = 'C:\isbe\nailfold\data\rsa_study\test/predictions/detection/rf_classification/222836/';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\test/predictions/orientation/rf_regression/222835/';
markup_dir = 'C:\isbe\nailfold\data\rsa_study\markup/tmoore/';

sc_list = dir([sc_dir '*.mat']);

clusters = 1:14;%[5 7 13];
for i_sc = 11:20%:10%length(sc_list)
    im_name = sc_list(i_sc).name(1:6);
    vessel_markup_list = dir([markup_dir '*' im_name '*.txt']);           
        
    load([vc_dir im_name '_vc.mat']);
    load([sc_dir im_name '_sc.mat'], 'sc_pi');
    
    vessel_im = u_load([image_dir im_name '.mat']);
    vessel_prob = u_load([prob_dir im_name '_pred.mat']);
    %vessel_ori = u_load([ori_dir im_name '_pred.mat']);
    vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
    
    include_pts = (vessel_centre_prob>.5) & (vessel_centre_curv>.05);
    vessel_centre_x = vessel_centre_x(include_pts);
    vessel_centre_y = vessel_centre_y(include_pts);
    
    num_pts = size(sc_pi, 3);
    cluster_idx_i = zeros(num_pts, 1);
    
    sc_pi = permute(sc_pi, [3 1 2]);
    sc_pi = [real(sc_pi(:,:)) imag(sc_pi(:,:))];
    
    for i_pt = 1:num_pts;
        c_dists = sum(bsxfun(@minus, cluster_centres, sc_pi(i_pt,:)).^2, 2);
        [~, cluster_idx_i(i_pt)] = min(c_dists);
    end
    
    figure; 
    a1 = subplot(2,1,1); imgray(vessel_im); caxis([100 200]);
    a2 = subplot(2,1,2); imgray(vessel_prob);
    linkaxes([a1 a2]);

    for i_c = clusters
        apex_cluster = cluster_idx_i == i_c;
        
        if i_c > 7
            plot(vessel_centre_x(apex_cluster), vessel_centre_y(apex_cluster), 'x');
        else
            plot(vessel_centre_x(apex_cluster), vessel_centre_y(apex_cluster), '+');
        end
    end
    legend(strcat('cluster ', cellstr(num2str(clusters'))));
    
    for i_v = 1:length(vessel_markup.vessels);

        %Check this is a valid vessel
        anchor_xy = vessel_markup.vessels(i_v).anchor;

        if isempty(anchor_xy); continue; end

        %Check if distal
        is_distal = vessel_markup.vessels(i_v).ncm_vessel_properties.is_distal;

        if is_distal
            num_apices = length(vessel_markup.vessels(i_v).apices);
            for i_a = 1:num_apices
                if isempty(vessel_markup.vessels(i_v).apices(i_a).inner_point)                    
                    %Use the anchor
                    plot(anchor_xy(1), anchor_xy(2), 'ko', 'markersize', 10);
                    
                else
                    %Compute the centre of the apex
                    apex_xy =  ...
                        mean([ vessel_markup.vessels(i_v).apices(i_a).outer_point;...
                          vessel_markup.vessels(i_v).apices(i_a).inner_point]);
                    plot(a1, apex_xy(1), apex_xy(2), 'rx');
                    plot(a2, apex_xy(1), apex_xy(2), 'k*', 'markersize', 10);
                                         
                end
            end
        else
            plot(anchor_xy(1), anchor_xy(2), 'kv', 'markersize', 10);
        end
    end
end
%%
x = 1814;
y = 478;

[sc_temp pc_temp path_map] = shape_context_prob_track_mult(...
                    vessel_prob.*vessel_ori, vessel_prob, vessel_ori,...
                    y, x,...
                    'num_streams', 1e4, 'stopping_prob', 0.2);
sc_temp2 = [real(sc_temp(:))' imag(sc_temp(:))'];
c_dists = sum(bsxfun(@minus, cluster_centres, sc_temp2).^2,2);
[~, c_idx] = min(c_dists);
display(c_idx);

figure; 
subplot(1,2,1); imgray(display_shape_context(sc_temp, 10, 'p'));
subplot(1,2,2); imgray(display_shape_context(sc_temp, 10, 'm'));
figure; imgray(path_map);
%%
vc_dir = 'C:\isbe\nailfold\data\rsa_study\test\vessel_centres\';
prob_dir = 'C:\isbe\nailfold\data\rsa_study\test/predictions/detection/rf_classification/222836/';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\test/predictions/orientation/rf_regression/222835/';
markup_dir = 'C:\isbe\nailfold\data\rsa_study\markup/tmoore/';
all_images = dir([prob_dir '*.mat']);

vessel_prob_smoothing_sigma = 2;
strong_vessel_thresh = 0.25;
g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);

%What is the probability of following a path
%

%We need a l;ist of all the images Tonia has actually made vessel paths
%for
im_list = cell(0,1);
for i_im = 1:length(all_images);
    im_name = all_images(i_im).name(1:6);
    vessel_markup_list = dir([markup_dir '*' im_name '*.txt']);
    if ~isempty(vessel_markup_list)
        vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
        
        centres_marked = false;
        for i_v = 1:length(vessel_markup.vessels);

            if ~isempty(vessel_markup.vessels(i_v).points)
                centres_marked = true;
                break;
            end
        end
        if centres_marked
            im_list{end+1,:} = im_name; %#ok
        end
    end
end
%%        
%Load in prob, ori, markup and centres
for i_im = 8%1:10%12;
    im_name = im_list{i_im};
    vessel_markup_list = dir([markup_dir '*' im_name '*.txt']);           

    try
        load([vc_dir im_name '_vc.mat']);
        vessel_prob = u_load([prob_dir im_name '_pred.mat']);
        %vessel_ori = u_load([ori_dir im_name '_pred.mat']);
        vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
    catch last_err
        display(last_err.message);
        display(['Skipping image ' num2str(i_im)]);
        continue;
    end

    %Smooth prob and ori
    vessel_prob = conv2(g', g, vessel_prob, 'same');
    vessel_ori = conv2(g', g, vessel_ori, 'same');
        
    %Find ends of centres
    vessel_centres = false(size(vessel_prob));
    vessel_centre_idx = sub2ind(size(vessel_prob), vessel_centre_y, vessel_centre_x);
    vessel_centres(vessel_centre_idx) = 1;
    vessel_centre_ends = bwmorph(vessel_centres, 'endpoints');
    [end_y end_x] = find(vessel_centre_ends);

    figure; imgray(vessel_prob);
    plot(vessel_centre_x, vessel_centre_y, 'b.');
    plot(end_x, end_y, 'r.');
    
    for i_v = 1:length(vessel_markup.vessels);

        if ~isempty(vessel_markup.vessels(i_v).anchor)
            plot(vessel_markup.vessels(i_v).anchor(1), vessel_markup.vessels(i_v).anchor(2), 'gx');
        end
        if ~isempty(vessel_markup.vessels(i_v).points)
            plot(vessel_markup.vessels(i_v).points(:,1), vessel_markup.vessels(i_v).points(:,2), 'g');
        end
    end
    
end

[~,~, path_map] = shape_context_prob_track_mult(...
    vessel_ori.*vessel_prob, vessel_prob, vessel_ori,...
    181, 944,...
    'num_streams', 1e4, 'stopping_prob', 0.1);
            
%Workout direction in which to travel

%Determine which ends lie on vessels

%if so

    %Compute how far to the end of the vessel each point lie
    
    %Compute probability of travelling along that path
    



    
    
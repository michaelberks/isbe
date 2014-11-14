rsa_dir = 'rsa_study/';

model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];
fov_mask_dir = ['C:\isbe\nailfold\data\' rsa_dir 'test\fov_masks\'];

markup_names = dir([nailfoldroot 'data/' rsa_dir 'markup/']);
markup_dirs = cell(0,1);
for i_dir = 3:length(markup_names)
    if markup_names(i_dir).isdir
        markup_dirs{end+1,1} = [nailfoldroot 'data/' rsa_dir 'markup/' markup_names(i_dir).name '/'];
    end
end

%
test_list = dir([fov_mask_dir '*_f_mask.mat']);
%
distal_count = zeros(length(test_list),1);

for i_im = 1:10%length(test_list)
    
    im_num = test_list(i_im).name(1:6);
    fov_mask = u_load([fov_mask_dir im_num '_f_mask.mat']);
    figure; imgray(fov_mask);
    [nrows ncols] = size(fov_mask);
    %
    for i_dir = 1:length(markup_dirs)
        vessel_markup_list_i = dir([markup_dirs{i_dir} '*' im_num '*.txt']);
        if ~isempty(vessel_markup_list_i)
            markup_dir = markup_dirs{i_dir};
        	vessel_markup_list = vessel_markup_list_i;
        end
    end
    if isempty(vessel_markup_list)
        distal_count(i_im) = -1;
    else
        
        %Read in markup
        vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
        
        %Loop through each marked vessel apex
        num_vessels = length(vessel_markup.vessels);
        distal_xy = zeros(0,2);
        for i_v = 1:num_vessels

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
                        plot(anchor_xy(:,1), anchor_xy(:,2), 'ro');

                    else
                        %Compute the centre of the apex
                        apex_xy =  ...
                            [ vessel_markup.vessels(i_v).apices(i_a).outer_point;...
                              vessel_markup.vessels(i_v).apices(i_a).inner_point];
                         
                        apex_width = sqrt(sum(diff(apex_xy).^2));
                        
                        distal_count(i_im) = distal_count(i_im) + 1;
%                         distal_y_hist(i_im,:) = distal_y_hist(i_im,:) +...
%                             hist(mean(apex_xy(:,2)) - mid_y, y_bins);
                        
                        distal_xy(end+1,:) = mean(apex_xy); %#ok;
                        if apex_width < 50
                            plot(apex_xy(:,1), apex_xy(:,2), 'g-*')
                        else
                            plot(apex_xy(:,1), apex_xy(:,2), 'r-*')
                        end
                        

                    end
                end
            else
                %Mask out the region around the anchor
                plot(anchor_xy(:,1), anchor_xy(:,2), 'r^');
                
            end
        end
    end
    if distal_count(i_im) > 1
        [p s] = polyfit(distal_xy(:,1), distal_xy(:,2), 1);
        distal_row_x = [1 ncols];
        distal_row_y = polyval(p, distal_row_x);
        plot(distal_row_x, distal_row_y);
    end
    %
    
end
%%
rsa_dir = 'rsa_study/';


centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];
corr_pos_dir = [nailfoldroot 'data/' rsa_dir 'test/corr_pos/'];
create_folder(corr_pos_dir);

markup_names = dir([nailfoldroot 'data/' rsa_dir 'markup/']);
markup_dirs = cell(0,1);
for i_dir = 3:length(markup_names)
    if markup_names(i_dir).isdir
        markup_dirs{end+1,1} = [nailfoldroot 'data/' rsa_dir 'markup/' markup_names(i_dir).name '/'];
    end
end
%
test_list = dir([centre_dir '*_vc.mat']);
%
for i_im = 1:length(test_list)
    
    im_num = test_list(i_im).name(1:6);
    
    load([centre_dir im_num '_vc.mat']);
    vessel_centre_mask = false(nrows, ncols);
    vessel_centre_idx = sub2ind([nrows ncols], vessel_centre_y, vessel_centre_x);
    vessel_centre_mask(vessel_centre_idx) = 1;
    %
    
    centre_labels = bwlabel(vessel_centre_mask, 8);
    centre_labels = centre_labels(vessel_centre_mask);
    
    apex_corr_max_points = false(length(vessel_centre_x),1);
    
    
    for i_lab = 1:max(centre_labels(:))
        label_mask = centre_labels == i_lab;     
        label_corr = vessel_centre_corr(label_mask);

        if any(label_corr > 0.4)

            label_apex_corr_max_points = false(size(label_corr));
            [~,max_idx] = max(label_corr);
            label_apex_corr_max_points(max_idx) = 1;
            apex_corr_max_points(label_mask) = label_apex_corr_max_points;         
        end
        
    end
    max_corr_apex_x = vessel_centre_x(apex_corr_max_points);
    max_corr_apex_y = vessel_centre_y(apex_corr_max_points);
    max_corr_vals = vessel_centre_corr(apex_corr_max_points);
    
    [maxima_corr_pos, maxima_corr_vals] = ...
        apply_local_exclusion([max_corr_apex_x max_corr_apex_y],...
        max_corr_vals, 50);
    
    for i_dir = 1:length(markup_dirs)
        vessel_markup_list_i = dir([markup_dirs{i_dir} '*' im_num '*.txt']);
        if ~isempty(vessel_markup_list_i)
            markup_dir = markup_dirs{i_dir};
        	vessel_markup_list = vessel_markup_list_i;
        end
    end
    
    if ~isempty(vessel_markup_list)
        
        distal_xy = zeros(0,2);
        non_distal_xy = zeros(0,2);
        undefined_xy = zeros(0,2);
        
        %Read in markup
        vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
        
        %Loop through each marked vessel apex
        num_vessels = length(vessel_markup.vessels);
        for i_v = 1:num_vessels

            %Check this is a valid vessel
            anchor_xy = vessel_markup.vessels(i_v).anchor;

            if isempty(anchor_xy); continue; end

            %Check if distal
            is_distal = vessel_markup.vessels(i_v).ncm_vessel_properties.is_distal;

            if is_distal
                
                num_apices = length(vessel_markup.vessels(i_v).apices);
                for i_a = 1:num_apices
                    if isempty(vessel_markup.vessels(i_v).apices(i_a).inner_point)                    
                        %Save the position of the anchor in the undefined holder
                        undefined_xy(end+1,:) = anchor_xy; %#ok

                    else
                        %Compute the centre of the apex and save in distal
                        %holder
                        apex_xy =  ...
                            [ vessel_markup.vessels(i_v).apices(i_a).outer_point;...
                              vessel_markup.vessels(i_v).apices(i_a).inner_point];
                         
                        distal_xy(end+1,:) = mean(apex_xy); %#ok 

                    end
                end
            else
                %Save the position of the anchor in the non-distal holder
                non_distal_xy(end+1,:) = anchor_xy; %#ok
                
            end
        end
    end
    %
    save([corr_pos_dir im_num '_cp.mat'], 'maxima_corr_pos', 'maxima_corr_vals', 'distal_xy', 'non_distal_xy', 'undefined_xy');
end
%%
rsa_dir = 'rsa_study/';
centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];
corr_pos_dir = [nailfoldroot 'data/' rsa_dir 'test/corr_pos/'];
test_list = dir([centre_dir '*_vc.mat']);
%
for i_im = 1:50
    
    im_num = test_list(i_im).name(1:6);
    
    load([centre_dir im_num '_vc.mat']);
    load([corr_pos_dir im_num '_cp.mat']);
    
    if size(distal_xy,1) >= 10 && size(maxima_corr_pos,1) >= 10
        
        [p] = polyfit(distal_xy(:,1), distal_xy(:,2), 1);
        
        dist_to_line = ( p(1)*distal_xy(:,1) - distal_xy(:,2) + p(2) ) / sqrt(p(1)^2 + 1);
        
        [~,min_idx] = min(dist_to_line);
        [~,max_idx] = max(dist_to_line);
        
        pu = p;
        pu(2) = distal_xy(max_idx,2) - pu(1)*distal_xy(max_idx,1);
        pl = p;
        pl(2) = distal_xy(min_idx,2) - pl(1)*distal_xy(min_idx,1);
        
        lx = [1 ncols];
        lyc = polyval(p, lx);
        lyu = polyval(pu, lx);
        lyl = polyval(pl, lx);
        
        
        
        figure; axis equal ij; hold all; 
        plot(distal_xy(:,1), distal_xy(:,2), 'gx')
        plot(non_distal_xy(:,1), non_distal_xy(:,2), 'rv');
        plot(undefined_xy(:,1), undefined_xy(:,2), 'r^');
        plot(maxima_corr_pos(:,1), maxima_corr_pos(:,2), 'mo');
        plot(maxima_corr_pos(1:10,1), maxima_corr_pos(1:10,2), 'bs');
        plot(lx, [1 1]*min(distal_xy(:,2)), 'g--');
        plot(lx, [1 1]*max(distal_xy(:,2)), 'g--');
        plot(lx, lyc, 'c-.');
        plot(lx, lyu, 'c--');
        plot(lx, lyl, 'c--');
    end
end
%%
rsa_dir = 'rsa_study/';
centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];
corr_pos_dir = [nailfoldroot 'data/' rsa_dir 'test/corr_pos/'];
test_list = dir([centre_dir '*_vc.mat']);
%
num_images = length(test_list);
X = zeros(num_images, 30);
y_min = zeros(num_images, 1);
y_max = zeros(num_images, 1);
discard_rows = false(num_images,1);

for i_im = 1:length(test_list)
    
    im_num = test_list(i_im).name(1:6);
    
    load([centre_dir im_num '_vc.mat']);
    load([corr_pos_dir im_num '_cp.mat']);
    
    if size(distal_xy,1) >= 10 && size(maxima_corr_pos,1) >= 10
       
       corr_x = maxima_corr_pos(1:10,1);
       corr_y = maxima_corr_pos(1:10,2);
       med_x = median(corr_x);
       med_y = median(corr_y);
       
       corr_x = corr_x - med_x;
       corr_y = corr_y - med_y;
       
       y_min(i_im) = min(distal_xy(:,2)) - med_y;
       y_max(i_im) = max(distal_xy(:,2)) - med_y;
       
       X(i_im,:) = [corr_x' corr_y' maxima_corr_vals(1:10)'];
        
    else
        discard_rows(i_im) = 1;
    end
end

X(discard_rows,:) = [];
y_min(discard_rows) = [];
y_max(discard_rows) = [];

%%
rf_args.prediction_type = 'rf_regression';
rf_args.n_trees = 100;
rf_args.d = [];
rf_args.w_prior = 0;
rf_args.impure_thresh = 1.0000e-004;
rf_args.split_min = 100;
rf_args.end_cut_min = 25;
rf_args.do_ubound = 0;
rf_args.quiet = 1;
rf_args.do_circular = [];
rf_args.overwrite = 1;
rf_args.minimise_size = 0;
rf_args.split_criterion = 'ssq';
rf_args.var_criterion = 'ssq';

rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_args.decomposition_args = [];

%Train
rf_dir = [corr_pos_dir '/rf/'];
rf_args.tree_dir = [rf_dir 'trees/'];
rf_args.sampling_args.X = X(1:200,:);


rf_args.sampling_args.y = y_min(1:200);
predictor = random_forest_reg_train(rf_args);
y_min_p = random_forest_reg_predict(predictor, X(201:end,:));
figure; plot(y_min(201:end), y_min_p, 'x');
axis equal;

rf_args.sampling_args.y = y_max(1:200);
predictor = random_forest_reg_train(rf_args);
y_max_p = random_forest_reg_predict(predictor, X(201:end,:));
figure; plot(y_max(201:end), y_max_p, 'x');
axis equal;


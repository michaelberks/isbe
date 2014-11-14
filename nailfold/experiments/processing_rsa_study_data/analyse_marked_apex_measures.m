function vessels_stats = analyse_marked_apex_measures()

shape_labels = {'Non-specific', 'Normal', 'Angiogenic','Meandering'};
size_labels = {'Normal', 'Enlarged', 'Giant', 'Irregular', 'Undefined'}; 

image_dir = 'Q:\nailfold\data\rsa_study\images\caps\';
load C:\isbe\nailfold\data\rsa_study\image_id_data.mat

for test_dirc = {'final_test'}%'test_half', 
    
    test_dir = test_dirc{1};
    
    apex_gt_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_gt\'];
    vessel_centre_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\vessel_centres\full_centres\'];

    im_list = dir([apex_gt_dir '*_gt.mat']);

    num_images = length(im_list);
    vessels_stats.(test_dir).category = zeros(num_images,1);
    vessels_stats.(test_dir).gradeable = false(num_images,1);
    vessels_stats.(test_dir).grade = cell(num_images, 1);

    vessels_stats.(test_dir).num_distal_vessels = zeros(num_images,1);
    vessels_stats.(test_dir).num_nondistal_vessels = zeros(num_images,1);
    vessels_stats.(test_dir).num_undefined_vessels = zeros(num_images,1);

    vessels_stats.(test_dir).mean_width = nan(num_images,1);
    vessels_stats.(test_dir).median_width = nan(num_images,1);
    vessels_stats.(test_dir).max_width = nan(num_images,1);
    vessels_stats.(test_dir).std_width = nan(num_images,1);

    vessels_stats.(test_dir).shape_counts = nan(num_images, 4);
    vessels_stats.(test_dir).size_counts = nan(num_images, 5);

    vessels_stats.(test_dir).vessel_density = nan(num_images,1);

    for i_im = 1:num_images
        im_name = im_list(i_im).name(1:6);
        im_idx = strcmp(image_id_data.im_names, im_name);
        switch image_id_data.category{im_idx}
            case 'S'
                vessels_stats.(test_dir).category(i_im) = 1;
            case 'P'
                vessels_stats.(test_dir).category(i_im) = 2;
            case 'HC'
                vessels_stats.(test_dir).category(i_im) = 3;
        end
        
        load([apex_gt_dir  im_name '_gt.mat']);

        vessels_stats.(test_dir).gradeable(i_im) = gradeable;
        vessels_stats.(test_dir).grade{i_im} = majority_grade;

        vessels_stats.(test_dir).num_distal_vessels(i_im) = sum(is_distal);
        vessels_stats.(test_dir).num_nondistal_vessels(i_im) = sum(is_non_distal);
        vessels_stats.(test_dir).num_undefined_vessels(i_im) = sum(is_undefined);

        if vessels_stats.(test_dir).num_distal_vessels(i_im)
            vessels_stats.(test_dir).mean_width(i_im) = mean(apex_widths(is_distal));
            vessels_stats.(test_dir).median_width(i_im) = median(apex_widths(is_distal));
            vessels_stats.(test_dir).max_width(i_im) = max(apex_widths(is_distal));
            vessels_stats.(test_dir).std_width(i_im) = std(apex_widths(is_distal));

            for i_sh = 1:length(shape_labels)
                vessels_stats.(test_dir).shape_counts(i_im, i_sh) = ...
                    sum(strcmpi(apex_shape, shape_labels{i_sh}) & is_distal);
            end
            for i_sz = 1:length(size_labels)
                vessels_stats.(test_dir).size_counts(i_im, i_sz) = ...
                    sum(strcmpi(apex_size, size_labels{i_sz}) & is_distal);
            end
            
            x1 = 1;
            if exist([vessel_centre_dir im_name '_vc.mat'], 'file');
                load([vessel_centre_dir im_name '_vc.mat'], 'ncols');
                x2 = ncols;
            else               
                f = imfinfo([image_dir im_name '.bmp']);
                x2 = f.Width;
            end
            
%             p = polyfit(apex_xy(is_distal,1), apex_xy(is_distal,2),1);
%             y1 = polyval(p, x1);
%             y2 = polyval(p, x2);
%             d1 = sqrt((x1-x2)^2 + (y1-y2)^2);

            d1 = x2 - x1;

            vessels_stats.(test_dir).vessel_density(i_im) = vessels_stats.(test_dir).num_distal_vessels(i_im) / d1;
            
        end
    end
    [vessels_stats.(test_dir).grade_idx vessels_stats.(test_dir).image_grade_labels] = grp2idx(vessels_stats.(test_dir).grade);
    vessels_stats.(test_dir).shape_pct = bsxfun(@rdivide, vessels_stats.test_half.shape_counts, vessels_stats.test_half.num_distal_vessels);
    vessels_stats.(test_dir).size_pct = bsxfun(@rdivide, vessels_stats.test_half.size_counts, vessels_stats.test_half.num_distal_vessels);
end
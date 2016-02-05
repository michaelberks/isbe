create_folder('C:\isbe\nailfold\data\wellcome_study\candidate_patches\');
patch_sz = 128;
patch_sz = patch_sz + 2; %Account for padding
patch_sz2 = (patch_sz - 1)/2;
base_width = 20;

% Get hog size from the hog_args struct
% Set up x,y coordinates for template patch
x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
y = x' + 32;
xy = [x(:) y(:)];

curr_can = 462;
%%
for i_sub = 61:10:111
    display(['Processing subject ' num2str(i_sub) '. ' num2str(curr_can-1) ' candidates saved.']);
    
    for i_seq = seq_nums(i_sub,1):seq_nums(i_sub,2)
        seq_dir = [sequence_dirs{i_seq} 'sequence_data\'];
        cap_dir = [sequence_dirs{i_seq} 'capillary_data\'];
        if ~exist([cap_dir 'apex_candidates.txt'], 'file')
            display([cap_dir ' full mosaic missing']);
            continue;
        end
        vessel_pred = load([cap_dir 'vessels_v_pred.txt']);
        ac = load([cap_dir 'apex_candidates.txt']);

        if isempty(ac)
            continue;
        end
        
        rescores = ac(:,7);
        selected_cans = rescores > 0.5;
        candidate_xy = ac(selected_cans,1:2);
        candidate_widths = ac(selected_cans,3);
        candidate_oris = atan2(ac(selected_cans,5), ac(selected_cans,4))/2;
        num_cans = sum(selected_cans);
        
        nailfold_image = double(imread([seq_dir 'full_mosaic.png']));
        nailfold_image = imresize(nailfold_image, size(vessel_pred), 'lanczos2');
        
%         figure; imgray(nailfold_image); a1 = gca;
        
        for i_can = 1:num_cans
            %Extract scaled/rotated patch                  
            ori_c = candidate_oris(i_can);
            width_c = candidate_widths(i_can);

            %Get scale relative to base width a make rotation matrix
            rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
            scale = width_c / base_width;
            cx = candidate_xy(i_can,1);
            cy = candidate_xy(i_can,2);
            
            %Transform points given scale and angle and translate to
            %candidate position
            xya = xy * rot * scale;
            xa = reshape(xya(:,1) + cx, patch_sz, patch_sz);
            ya = reshape(xya(:,2) + cy, patch_sz, patch_sz);
%             plot(a1, xa, ya, 'b.', 'markersize', 2);
            
            %Sample vessel prob patch
            vessel_feature_patch = interp2(nailfold_image, xa, ya, '*linear', 0);
            vessel_prob_patch = interp2(vessel_pred, xa, ya, '*linear', 0);
            
            save(['C:\isbe\nailfold\data\wellcome_study\candidate_patches\candidate' zerostr(curr_can,5) '.mat'],...
                'vessel_feature_patch', 'vessel_prob_patch', 'cap_dir', 'ori_c', 'width_c', 'base_width', 'cx', 'cy');
            
%             r = rem(curr_can - 1, 24);
%             if ~r
%                 f1 = figure;
%             end           
%             figure(f1); subplot(4,6,r+1); imgray(vessel_feature_patch);
            
            curr_can = curr_can + 1;
        end
    end
end
%%
load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat', 'miccai_selection');
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
im_names = sort(image_id_data.im_names(miccai_selection.validation));

create_folder('C:\isbe\nailfold\data\rsa_study\master_set\candidate_patches\');
patch_sz = 64;
patch_sz = patch_sz + 2; %Account for padding
patch_sz2 = (patch_sz - 1)/2;
base_width = 10;

% Get hog size from the hog_args struct
% Set up x,y coordinates for template patch
x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
y = x';
xy = [x(:) y(:)];

im_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\images\';
vessel_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\predictions\detection\rf_classification\296655\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\predictions\orientation\rf_regression\296621\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\predictions\width\rf_regression\297037\';
rescore_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\apex_maps\set12g_half_296655\miccai_maxima\rescores\';

curr_can = 1;
for i_im = 1:3%length(im_names)
    display(['Processing image ' num2str(i_im) '. ' num2str(curr_can-1) ' candidates saved.']);
    
    nailfold_image = u_load([im_dir im_names{i_im} '.mat']);
    vessel_pred = u_load([vessel_dir im_names{i_im} '_pred.mat']);
    vessel_ori = u_load([ori_dir im_names{i_im} '_pred.mat']);
    vessel_width = u_load([width_dir im_names{i_im} '_pred.mat']);
    candidates = load([rescore_dir im_names{i_im} '_candidates.mat']);

    if isempty(candidates.candidate_xy)
        continue;
    end   
        
    selected_cans = candidates.candidate_rescores > 0.5;
    candidate_xy = candidates.candidate_xy(selected_cans,:);
    
    candidate_oris = interp2(vessel_ori, ...
        candidate_xy(:,1), candidate_xy(:,2));
    candidate_widths = interp2(vessel_width,...
        candidate_xy(:,1), candidate_xy(:,2));
    candidate_oris = angle(candidate_oris) / 2; 
        
    num_cans = sum(selected_cans);
        
    figure; imgray(nailfold_image); a1 = gca;
        
    for i_can = 1:num_cans
        %Extract scaled/rotated patch                  
        ori_c = candidate_oris(i_can);
        width_c = candidate_widths(i_can);

        %Get scale relative to base width a make rotation matrix
        rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
        scale = width_c / base_width;
        cx = candidate_xy(i_can,1);
        cy = candidate_xy(i_can,2);

        %Transform points given scale and angle and translate to
        %candidate position
        xya = xy * rot * scale;
        xa = reshape(xya(:,1) + cx, patch_sz, patch_sz);
        ya = reshape(xya(:,2) + cy, patch_sz, patch_sz);
        plot(a1, xa, ya, 'b.', 'markersize', 2);

        %Sample vessel prob patch
        vessel_feature_patch = interp2(nailfold_image, xa, ya, '*linear', 0);
        vessel_prob_patch = interp2(vessel_pred, xa, ya, '*linear', 0);

        save(['C:\isbe\nailfold\data\rsa_study\master_set\candidate_patches\candidate' zerostr(curr_can,5) '.mat'],...
            'vessel_feature_patch', 'vessel_prob_patch', 'cap_dir', 'ori_c', 'width_c', 'base_width', 'cx', 'cy');

        r = rem(curr_can - 1, 24);
        if ~r
            f1 = figure;
        end           
        figure(f1); subplot(4,6,r+1); imgray(vessel_feature_patch);

        curr_can = curr_can + 1;
    end

end
%work 29/05/2008
dt_list = dir('C:\isbe\dev\background\dual_tree\normal512\*dual_tree*');

p_sum_dt = zeros(60);
map_sum_dt = zeros(60);

for jj = 1:length(dt_list) % = number of regions

    %Load tree
    dual_tree = u_load(['C:\isbe\dev\background\dual_tree\normal512\', dt_list(jj).name]);
    
    % Get data from dual tree:

    %compute num_points;
    sample_mask = ones(size(dual_tree{1}(:,:,1)));
    sample_mask(1:32, :) = 0;
    sample_mask(end-31:end, :) = 0;
    sample_mask(:,1:32) = 0;
    sample_mask(:,end-31:end) = 0;

    [r1_pts c1_pts] = find(sample_mask);

    data_dt = zeros(size(r1_pts,1), 60);

    for level = 1:5 

        r_pts = ceil(r1_pts/2^(level-1));
        c_pts = ceil(c1_pts/2^(level-1));

        idx = sub2ind(size(dual_tree{level}(:,:,1)), r_pts, c_pts);

        for ori = 1:6

            idx_real = (level - 1)*6 + ori;
            idx_imag = (level - 1)*6 + ori + 30;
            temp = dual_tree{level}(:,:,ori);
            data_dt(:,idx_real) = real(temp(idx));
            data_dt(:,idx_imag) = imag(temp(idx));
            clear temp;
        end 
    end

    clear r* c* sample_mask idx*

    %
    [rho p_val] = corr(data_dt);
    
    p_sum_dt = p_sum_dt + p_val;
    map_sum_dt = map_sum_dt + (p_val < 0.01);
    
%     clear sig_cor;
%     
%     [sig_cor(:,1) sig_cor(:,2)] = find(p_val < 0.01);
%     sig_cor(sig_cor(:,1)==sig_cor(:,2),:) = [];
%     sig_cor_list = cell(60);
% 
%     for ii = 1:60
% 
%         temp = sig_cor(sig_cor(:,1) == ii, 2);
%         sig_cor_list{ii}(:,1) = fix((temp-1) ./ 6) + 1;
%         sig_cor_list{ii}(:,2) = rem(temp - 1, 6) + 1;
%         sig_cor_list{ii}(:,3) = sig_cor_list{ii}(:,1) > 5;
%         sig_cor_list{ii}(sig_cor_list{ii}(:,1)>5,1) = ...
%             sig_cor_list{ii}(sig_cor_list{ii}(:,1)>5,1) - 5;
% 
%     end

end

p_sum_dt = p_sum_dt / jj;
%%

%work 29/05/2008
pyr_list = dir('C:\isbe\dev\background\pyramid\normal512\*pyr*');

p_sum_pyr = zeros(20);
map_sum_pyr = zeros(20);

for jj = 1:length(pyr_list) % = number of regions

    %Load tree
    pyr = u_load(['C:\isbe\dev\background\pyramid\normal512\', pyr_list(jj).name]);

    % Get data from dual tree:

    % Set overlap of 64;

    %compute num_points;
    sample_mask = ones(size(pyr{2,1}));
    sample_mask(1:32, :) = 0;
    sample_mask(end-31:end, :) = 0;
    sample_mask(:,1:32) = 0;
    sample_mask(:,end-31:end) = 0;

    [r1_pts c1_pts] = find(sample_mask);

    data_pyr = zeros(size(r1_pts,1), 20);

    for level = 1:5 

        r_pts = ceil(r1_pts/2^(level-1));
        c_pts = ceil(c1_pts/2^(level-1));

        idx = sub2ind(size(pyr{level+1,1}), r_pts, c_pts);

        for ori = 1:4

            idx_real = (level - 1)*4 + ori;
            data_pyr(:,idx_real) = pyr{level+1,ori}(idx);

            clear temp;
        end 
    end

    clear r* c* sample_mask idx*


    %
    [rho p_val] = corr(data_pyr);

    p_sum_pyr = p_sum_pyr + p_val;
    map_sum_pyr = map_sum_pyr + (p_val < 0.01);
    clear data_pyr;

%     clear sig_cor
%     [sig_cor(:,1) sig_cor(:,2)] = find(p_val < 0.01);
%     sig_cor(sig_cor(:,1)==sig_cor(:,2),:) = [];
%     sig_cor_list = cell(60);
% 
%     for ii = 1:60
% 
%         temp = sig_cor(sig_cor(:,1) == ii, 2);
%         sig_cor_list{ii}(:,1) = fix((temp-1) ./ 6) + 1;
%         sig_cor_list{ii}(:,2) = rem(temp - 1, 6) + 1;
%         sig_cor_list{ii}(:,3) = sig_cor_list{ii}(:,1) > 5;
%         sig_cor_list{ii}(sig_cor_list{ii}(:,1)>5,1) = ...
%             sig_cor_list{ii}(sig_cor_list{ii}(:,1)>5,1) - 5;
% 
%     end

end

p_sum_pyr = p_sum_pyr / jj;
%%

pyr_list = dir('C:\isbe\dev\background\pyramid\normal512\*pyr*');
map_sum_pyr_all = zeros(20, 49);

for level = 2:6
    for ori = 1:4
        
        p_sum_pyr = zeros(49);
        map_sum_pyr = zeros(49);

        for jj = 1:length(pyr_list) % = number of regions
            
            display(['Processing mass ', num2str(jj)]);
            %Load tree
            pyr = u_load(['C:\isbe\dev\background\pyramid\normal512\', pyr_list(jj).name]);
    
            % Get data from pyr tree:
            overlap = 2^(8-level);
            sample_mask = ones(size(pyr{level,1}));
            sample_mask(1:overlap, :) = 0;
            sample_mask(end-overlap+1:end, :) = 0;
            sample_mask(:,1:overlap) = 0;
            sample_mask(:,end-overlap+1:end) = 0;

            [r1_pts c1_pts] = find(sample_mask);

            data_pyr = zeros(size(r1_pts,1), 49);
    
            for curr_point = 1:length(r1_pts)

                r = r1_pts(curr_point);
                c = c1_pts(curr_point);
                data_pyr(curr_point,:) = ...
                    reshape(sample_window(pyr{level,ori}, 7, r, c), 1, []);

            end
            
            %
            [rho p_val] = corr(data_pyr);

            p_sum_pyr = p_sum_pyr + p_val;
            map_sum_pyr = map_sum_pyr + (p_val < 0.01);
            clear data_pyr;
        end
        p_sum_pyr = p_sum_pyr / jj;
        map_sum_pyr_all(4*(level-2)+ori, :) = map_sum_pyr(25,:);
        figure; imagesc(reshape(map_sum_pyr(25,:),7,7)); axis image; caxis([0 89]);
        
        display(['Level ', num2str(level), ' Orientation ', num2str(ori), ' complete.']);
    end
end
%%
dt_list = dir('C:\isbe\dev\background\dual_tree\normal512\*dual_tree*');
map_sum_dt_local = zeros(30, 49);

for level = 1:5
    for ori = 1:6
        
        p_sum_dt = zeros(49);
        map_sum_dt = zeros(49);

        for jj = 1:length(dt_list) % = number of regions
            
            display(['Processing mass ', num2str(jj)]);
            %Load tree
            dual_tree = u_load(['C:\isbe\dev\background\dual_tree\normal512\', dt_list(jj).name]);
    
            % Get data from pyr tree:
            overlap = 2^(7-level);
            sample_mask = ones(size(dual_tree{level}(:,:,1)));
            sample_mask(1:overlap, :) = 0;
            sample_mask(end-overlap+1:end, :) = 0;
            sample_mask(:,1:overlap) = 0;
            sample_mask(:,end-overlap+1:end) = 0;

            [r1_pts c1_pts] = find(sample_mask);

            data_dt = zeros(size(r1_pts,1), 49);
    
            for curr_point = 1:length(r1_pts)

                r = r1_pts(curr_point);
                c = c1_pts(curr_point);
                data_dt(curr_point,:) = ...
                    reshape(sample_window(real(dual_tree{level}(:,:,ori)), 7, r, c), 1, []);

            end
            
            %
            [rho p_val] = corr(data_dt);

            p_sum_dt = p_sum_dt + p_val;
            map_sum_dt = map_sum_dt + (p_val < 0.01);
            clear data_pyr;
        end
        p_sum_dt = p_sum_dt / jj;
        map_sum_dt_local(6*(level-1)+ori, :) = map_sum_dt(25,:);
        figure; imagesc(reshape(map_sum_dt(25,:),7,7)); axis image; caxis([0 89]);
        
        display(['Level ', num2str(level), ' Orientation ', num2str(ori), ' complete.']);
    end
end


%%

[P_all, L_all, B_all, mu_all, R_all] = st_pca(data_dt,1);
%%
for ori = 1:6
    [P_level{ori}, L_level{ori}, B_level{ori}, mu_level{ori}, R_levels{ori}]...
        = st_pca(data_dt(:,ori:6:60),1);
    
    figure; bar(L_level{ori}); title(['Orientation = ', num2str(ori)]);
    set(gca,'XTickLabel',cumsum(L_level{ori})/sum(L_level{ori}))
end
%%
for level = 1:5
    [P_ori{level}, L_ori{level}, B_ori{level}, mu_ori{level}, R_ori{level}]...
        = st_pca(data_dt(:,[6*(level-1)+1:6*level 30+(6*(level-1)+1:6*level)]),1);
    
    figure; bar(L_ori{level}); title(['Level = ', num2str(level)]);
    set(gca,'XTickLabel',cumsum(L_ori{level})/sum(L_ori{level}))
end

%%
im_list = dir('C:\isbe\dev\background\images\normal512\o04_*');
%%
for ii = 1:length(im_list)
    
    %for each image
    i1 = double(imread(['C:\isbe\dev\background\images\normal512\', im_list(ii).name]));
    
    [gp p_sizes] = buildGpyr(i1, 5);
    gp = mb_change_pyramid_form(gp, p_sizes, 'g');
    
    CLS_results = cell(3,1);
    for jj = 1:3
        %Do CLS detection
        [CLS_results{jj}] = mb_cls_selection('ImageIn', gp{jj});
    end
    
    save(['C:\isbe\dev\background\cls\normal512\normal_cls_', im_list(ii).name(5:10)], 'CLS_results');    
    clear i1 r c CLS_results cls_mr;
end
%%
copyfile('C:\isbe\dev\background\models\pyramid\normal512_k10_w1_11_w2_0_2_1.mat', 'K:\isbe\normal512_k10_w1_11_w2_0_2_1.mat')

    



        
        
        


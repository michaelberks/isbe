%ACR script

%Load in the data
load('C:\isbe\nailfold\data\2_year_study\results\miccai\auto_stats.mat');
load('C:\isbe\nailfold\data\2_year_study\results\miccai\people_stats.mat');

%Check what images we have each subject at each timepoint
present_mask = squeeze(sum(sum(people_stats.present,2),3));

%As long as no values in the present mask are greater than one, we can use
%the same trick of summing across 2nd and 3rd dimensions to get 2D matrix
%for each measurement
any(present_mask(:) > 1)

%Valid subjects are ones with images at v1, v5 and v6 (0, 1 and 2 years)
valid_subjects = all(present_mask(:,[1 5 6]),2);

%Recategorise
recats = people_stats.super_category;
recats(people_stats.category == 4) = 2;
recats(people_stats.super_category == 2) = 3;
recats(people_stats.category == 5) = 4;

cat_idx = false(length(people_stats.category),4);
for i_cat = 1:4
    cat_idx(:,i_cat) = recats == i_cat & valid_subjects;%(people_stats.category == i_cat)
end
cat_names = {'HC', 'PR', 'SSc',  'UD'};%, 'dSSc', 'lSSc', 
cat_colors = 'rgby';
%%    
selected_features = {...
    'mean_weighted_width',...
    'mean_orientation_entropy',...
   'mean_inter_capillary_distance'};%'mean_max_width',... %'vessel_density2',...
    
    
                
analyse_qseries_apex_measures(auto_stats,...
    'image_names',          [],...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'vessel_centre_dir',    'vessel_centres\',...
    'selected_dir',         'apex_maps\set12g_half_296655\miccai_maxima\selected_apexes',...
    'metrics_dir',          'apex_maps\set12g_half_296655\miccai_maxima\apex_metrics',...
    'do_xls',           0,...
    'do_auto_stats',    0,...
    'do_people_stats',   0,...
    'do_image_plots',   1, ...
    'selected_features', selected_features,... %Do them all for now
    'do_people_plots',  0,...
    'fig_dir', [],...
    'xls_filename', [],...
    'save_dir', []);
%%
im_recats = auto_stats.super_category;
im_recats(auto_stats.category == 4) = 2;
im_recats(auto_stats.super_category == 2) = 3;
im_recats(auto_stats.category == 5) = 4;

grps = {'HC', 'PR', 'SSc', 'UD'};
grp_comps = {'HC,PR', 'HC,SSc', 'HC,UD', 'PR,SSc', 'PR,UD', 'SSc,UD'};
grp_cols = 'rgbm';
feature_names = {'Mean capillary width', 'Mean capillary tortuosity', 'Capillary density'};
xlabels = {'log Width (mm)', 'Tortuosity', 'Density (mm^{-1})'};

fig_dir = 'P:\isbe\conferences_and_symposia\2014\submissions\acr\figs\';

for i_f = 1:3
      
    feature = selected_features(i_f);
    feature_str = feature{1};
    feature_str(feature_str == '_') = ' ';
    feature_str(1) = feature_str(1) - 32;
    
    display(['Tests for ' feature_str]);
    
    feature_scores = auto_stats.(feature{1})(:);    
    switch feature{1}
        case 'mean_inter_capillary_distance'
            feature_scores = 1./ feature_scores;
        case {'mean_weighted_width', 'max_mean_width'}
            feature_scores = log(feature_scores);
    end
    
    valid_scores = isnumeric(feature_scores) & isfinite(feature_scores) & ~isnan(feature_scores);
    
    [p,table,stats] = anova1(feature_scores(valid_scores), im_recats(valid_scores));
    [c,m] = multcompare(stats);
    m(:,3) = m(:,1)+1.96*m(:,2);
    m(:,2) = m(:,1)-1.96*m(:,2);
    
    figure('WindowStyle', 'normal');
    set(gca, 'fontsize', 18);
    hold all;
    title(feature_names{i_f});
    ylabel('f(x)');
    xlabel(xlabels{i_f});
    
    min_val = min(feature_scores(valid_scores));
    max_val = max(feature_scores(valid_scores));
    grid_pts = linspace(min_val, max_val, 100);

    max_y = 0;
    
    leg_text = cell(4,1);
    for i_g = 1:4
        kdist = build_1d_kernel_distribution(feature_scores(valid_scores & im_recats==i_g), grid_pts, 0);
        plot(kdist.x, kdist.D_f, grp_cols(i_g), 'linewidth', 2);
        max_y = max(max_y, max(kdist.D_f));
        leg_text{i_g} = [grps{i_g} ' (n = ' num2str(sum(im_recats==i_g)) ')'];
    end
    if i_f == 3
        legend(leg_text, 'location', 'northwest');
    end
    
    for i_g = 1:4
        plot(m(i_g,1), max_y*(1+i_g/20), [grp_cols(i_g) 'o'], 'markersize', 6,  'markerfacecolor', grp_cols(i_g));
        plot([m(i_g,2) m(i_g,3)], [max_y*(1+i_g/20) max_y*(1+i_g/20)], [grp_cols(i_g) '-'], 'linewidth', 2);
        plot([m(i_g,2) m(i_g,2)], [max_y*(1+i_g/21) max_y*(1+i_g/19)], [grp_cols(i_g) '-'], 'linewidth', 2);
        plot([m(i_g,3) m(i_g,3)], [max_y*(1+i_g/21) max_y*(1+i_g/19)], [grp_cols(i_g) '-'], 'linewidth', 2);
    end
    %exportfig([fig_dir 'f' num2str(i_f) '.png']);
    switch feature{1}
        case 'mean_inter_capillary_distance'
            %m = 1./ m;
        case {'mean_weighted_width', 'max_mean_width'}
            m = exp(m);
    end
    
    
    
    for i_g = 1:4
        display([grps{i_g} ': ' num2str(m(i_g,1),3) ' (' num2str(m(i_g,2),3) ', ' num2str(m(i_g,3),3) ')']);
    end
    disp_str = ['Sig diffs between: '];
    for i_gc = 1:6
        if c(i_gc,6) < 0.05
            disp_str = [disp_str grp_comps{i_gc} '(' num2str(c(i_gc,6),2) '); '];
        end
    end
    display(disp_str);
end
%%
im = [];
for i_f = [3 1 2]
    im_i = imread([fig_dir 'f' num2str(i_f) '.png']);
    if i_f == 3
        r = size(im_i,1);
        c = size(im_i,2);
    else
        im_i = imresize(im_i, [r c]);
    end
    im = cat(2, im, im_i);
end
figure; imgray(im);
imwrite(im, [fig_dir 'all_features.png']);
%%
for feature = selected_features
      
    
    feature_str = feature{1};
    feature_str(feature_str == '_') = ' ';
    feature_str(1) = feature_str(1) - 32;
    
    display(['Tests for ' feature_str]);
    
    feature_scores = people_stats.(feature{1})(:,:,:,[1 5 6]);
    feature_scores(~people_stats.present(:,:,:,[1 5 6])) = 0;
      
    feature_scores = squeeze(sum(sum(feature_scores,2),3));
    
    switch feature{1}
        case 'mean_inter_capillary_distance'
            feature_scores = 1./ feature_scores;
        case {'mean_weighted_width', 'max_mean_width'}
            %feature_scores = log(feature_scores);
    end
    
    valid_scores = all(isnumeric(feature_scores),2) & all(isfinite(feature_scores),2) & all(~isnan(feature_scores),2);

    all_scores = feature_scores(valid_subjects & valid_scores,:);
    all_cats = people_stats.category(valid_subjects & valid_scores);
    
    t = table(all_cats, all_scores(:,1),all_scores(:,2),all_scores(:,3),...
        'VariableNames',{'Categories','t1','t2','t3'});
    time = dataset([0 12 24]','VarNames',{'Time'});
    rm = fitrm(t,'t1-t3~Categories','WithinDesign',time);
    ranovatbl = ranova(rm)
    
%     figure;    
%     ax = zeros(2,1);

    fx(1) = figure('WindowStyle', 'normal');
    ax(1) = gca; hold all;
    fx(2) = figure('WindowStyle', 'normal');
    ax(2) = gca; hold all;
    
    for i_sp = 1:2%4
        %ax(i_sp) = subplot(1,2,i_sp); hold all;
        set(ax(i_sp), 'xlim', [-1 26]);  
        %title(feature_str);   
    end

    set(ax(1), ...
        'xtick', [3 13 23],...
        'xticklabel', {'0yr', '1yr', '2yr'},...
        'fontsize', 18)
    xlabel(ax(1), 'Time points');
    ylabel(ax(1), feature_str);
    
    plot(ax(2), [-1 26], [0 0], 'k--');
    set(ax(2), ...
        'xtick', [3 13 23],...
        'xticklabel', {'1yr - 0yr', '2yr - 1yr', '2yr - 0yr'},...
        'fontsize', 18);
    xlabel(ax(2), 'Time point comparisons');
    ylabel(ax(2), ['Change in ' feature_str]);
    
    for i_cat = 1:4
        cat_idx_i = cat_idx(:,i_cat) & valid_scores;
    
        n_cat = sum(cat_idx_i);
        
        cat_scores = feature_scores(cat_idx_i,:);
        cat_diffs = ([diff(cat_scores,1,2) cat_scores(:,3)-cat_scores(:,1)]);
        
        [~,p12] = ttest(cat_scores(:,1), cat_scores(:,2));
        [~,p13] = ttest(cat_scores(:,1), cat_scores(:,3));
        [~,p23] = ttest(cat_scores(:,2), cat_scores(:,3));   
    
        display([ cat_names{i_cat} ': ' num2str([p12 p13 p23],2)]);
        
        t = table(cat_scores(:,1),cat_scores(:,2),cat_scores(:,3),...
            'VariableNames',{'t1','t2','t3'});
        time = dataset([0 12 24]','VarNames',{'Time'});
        rm = fitrm(t,'t1-t3~1','WithinDesign',time);
        ranovatbl = ranova(rm)

        mu_s = mean(cat_scores);
        mu_d = mean(cat_diffs);
        se_s = std(cat_scores,1)/ sqrt(n_cat);
        se_d = std(cat_diffs,1)/ sqrt(n_cat);
    
        iqr_s = prctile(cat_scores, [25 50 75]);
        iqr_d = prctile(cat_diffs, [25 50 75]);
        
        %Plot1
        %boxplot(ax(1), cat_scores, 'positions', i_cat+[0 10 20],'colors', cat_colors(i_cat), 'widths', 0.5);
    
        %Plot2
        plot(ax(1), [0 10 20]+i_cat, mu_s, [cat_colors(i_cat) 'o'], 'markerfacecolor', cat_colors(i_cat), 'markersize', 10);
        plot(ax(1), [0 10 20]+i_cat, mu_s, cat_colors(i_cat), 'linewidth', 3);
        plot(ax(1), [0 0 0]+i_cat, mu_s(1)+se_s(1)*[-1.96 0 1.96], cat_colors(i_cat), 'linewidth', 2);
        plot(ax(1), [10 10 10]+i_cat, mu_s(2)+se_s(2)*[-1.96 0 1.96], cat_colors(i_cat), 'linewidth', 2);
        plot(ax(1), [20 20 20]+i_cat, mu_s(3)+se_s(3)*[-1.96 0 1.96], cat_colors(i_cat), 'linewidth', 2);
        
        %Plot3
        %boxplot(ax(3), cat_diffs, 'positions', i_cat+[0 10 20],'colors', cat_colors(i_cat), 'widths', 0.5);
    
        %Plot4
        plot(ax(2), [0 10 20]+i_cat, mu_d, [cat_colors(i_cat) 'o'], 'markerfacecolor', cat_colors(i_cat), 'markersize', 10);
        plot(ax(2), [0 0 0]+i_cat, mu_d(1)+se_d(1)*[-1.96 0 1.96], cat_colors(i_cat), 'linewidth', 3);
        plot(ax(2), [10 10 10]+i_cat, mu_d(2)+se_d(2)*[-1.96 0 1.96], cat_colors(i_cat), 'linewidth', 3);
        plot(ax(2), [20 20 20]+i_cat, mu_d(3)+se_d(3)*[-1.96 0 1.96], cat_colors(i_cat), 'linewidth', 3);

    end
    figure(fx(1));
    exportfig(['C:\isbe\nailfold\data\2_year_study\results\figures\time_change1_' feature{1} '.png']);
    figure(fx(2));
    exportfig(['C:\isbe\nailfold\data\2_year_study\results\figures\time_change2_' feature{1} '.png']);
end
%%
figure('WindowStyle', 'normal'); hold all;
leg_caption = cell(4,1);
for i_cat = 1:4
    n_cat = sum(cat_idx(:,i_cat) & valid_scores);
    plot(-1, -1, cat_colors(i_cat), 'linewidth', 3);
    leg_caption{i_cat} = [cat_names{i_cat} ' (n = ' num2str(n_cat) ')'];
end
axis off;
axis([0 1 0 1]);
legend(leg_caption, 'location', 'northeast');
exportfig('C:\isbe\nailfold\data\2_year_study\results\figures\time_change_legend.png');
    
%%
results_dir = 'C:\isbe\nailfold\data\2_year_study\results\';
auto_stats = analyse_qseries_apex_measures(auto_stats,...
    'image_names',          im_names,...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'vessel_centre_dir',    'vessel_centres\',...
    'selected_dir',         'apex_maps\set12g_half_296655\miccai_maxima\selected_apexes',...
    'metrics_dir',          'apex_maps\set12g_half_296655\miccai_maxima\apex_metrics',...
    'do_xls',           0,...
    'do_auto_stats',    1,...
    'do_people_stats',   1,...
    'do_image_plots',   0, ...
    'selected_features', [],... %Do them all for now
    'do_people_plots',  0,...
    'fig_dir', [],...
    'xls_filename', [],...
    'save_dir', results_dir);
%%
extract_apex_measures_set(...
    'image_names',          {'c02V1LD4X3LrgMosaic'},...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'image_dir',            'images',...
    'prob_dir',             'rf_classification/296655',...
    'ori_dir',              'rf_regression/296621',...
    'width_dir',            'rf_regression/297037',...
    'candidates_dir',       'apex_maps\set12g_half_296655\miccai_maxima\rescores',...
    'displacements_dir',    'apex_maps\set12g_half_296655\miccai_maxima\displacements',...
    'selected_dir',         'apex_maps\set12g_half_296655\miccai_maxima\selected_apexes',...
    'metrics_dir',          [],...
    'width_predictor',      [],...[nailfoldroot 'models/apex/width/rf.mat'],...
    'prob_sigma',           2,...
    'ori_sigma',            0,...
    'width_sigma',          2,...
    'all',                  0,...
    'do_aam',               0,...
    'overwrite',            0,...
    'plot', 1);
%%
temp_names = {'l44V1LD4X3LrgMosaic', 'l44V6LD4X3LrgMosaic'};
load([nailfoldroot,'models/apex/final_MAP/miccai_class_MAP/class_map.mat']);
select_vessels_from_candidates( ...
    'image_names',          temp_names,...
    'data_dir', [nailfoldroot 'data/2_year_study/'],...
    'class_map',            class_map,...
    'selected_dir',         'temp',...'apex_maps\set12g_half_296655\miccai_maxima\selected_apexes',...
    'displacement_dir',     'apex_maps\set12g_half_296655\miccai_maxima\displacements',...
    'rescore_dir',          'apex_maps\set12g_half_296655\miccai_maxima\rescores',...
    'strong_vessel_thresh', 0.8,...
    'angle_discard_thresh', 75*pi/180,...
    'do_final_cull',        1,...
    'do_post_merge', 1,...
    'merge_dist_thresh', 60,...
    'merge_connect_thresh', 0.25,...
    'merge_n_pts', 20);
%%
extract_apex_measures_set(...
    'image_names',          temp_names,...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'image_dir',            'images',...
    'prob_dir',             'rf_classification/296655',...
    'ori_dir',              'rf_regression/296621',...
    'width_dir',            'rf_regression/297037',...
    'candidates_dir',       'apex_maps\set12g_half_296655\miccai_maxima\rescores',...
    'displacements_dir',    'apex_maps\set12g_half_296655\miccai_maxima\displacements',...
    'selected_dir',         'temp',...
    'metrics_dir',          'temp',...
    'aam_dir',              'aam',...
    'do_aam',               0,...
    'model_dir',            [nailfoldroot 'data/rsa_study/models/apex_templates/'],...
    'aam_name',             'aam/orig/2/vessel_apex_orig.smd',...
    'width_predictor',      [],...[nailfoldroot 'models/apex/width/rf.mat'],...
    'prob_sigma',           2,...
    'ori_sigma',            0,...
    'width_sigma',          2,...
    'do_distal',            1,...
    'do_nondistal',         1,...
    'plot', 0);
%%
for i_p = [176   179   180   184   186   188   189   191]

    s_cat = people_stats.super_category(i_p);
    if ~s_cat; continue; end

    baseline_im = squeeze(people_stats.present(i_p,:,:,1) & ~people_stats.status(i_p,:,:,1));

    [digit hand] = find(baseline_im);

    if isempty(digit)
        continue;
    elseif length(digit) > 1
        digit = digit(1);
        hand = hand(1);
    end

    followup_visits = find(squeeze(people_stats.present(i_p,digit,hand,5:end) & ~people_stats.status(i_p,digit,hand,5:end)))+1;       
    cat = people_stats.category(i_p);

    switch cat

        case 1
            cat_str = 'c';
        case 2
            cat_str = 'd';
        case 3
            cat_str = 'l';
        case 6
            cat_str = 'n';
        case 4
            cat_str = 'p';

        otherwise
            continue;               
    end        
    people_str = [cat_str zerostr(rem(people_stats.people_ids(i_p),1000),2)];
    visit_str = '156';
    switch hand
        case 1
            hand_str = 'L';
        case 2
             hand_str = 'R';
    end 
    digit_str = ['D' num2str(digit) 'X3LrgMosaic'];

    visit_names = cell(3,1);
    for i_v = 1:3
        if ismember(i_v, [1; followup_visits])
            visit_names{i_v} = [people_str 'V' visit_str(i_v) hand_str digit_str];
        end               
    end
    
    display_automated_markup(... % non-strict mode
        'image_names',         visit_names([1 3]),...
        'data_dir',             [nailfoldroot 'data/2_year_study/'],...
        'image_dir',            'images',...
        'vessel_centre_dir',    'vessel_centres',...
        'metrics_dir',          'apex_maps\set12g_half_296655\miccai_maxima\apex_metrics',...
        'candidates_dir',       'apex_maps\set12g_half_296655\miccai_maxima\rescores',...
        'selected_dir',         'apex_maps\set12g_half_296655\miccai_maxima\selected_apexes',...
        'selected_features', [],...
        'do_xls',           1,...
        'do_make_stats',    0,...
        'do_image_plots',   0, ...
        'do_people_plots',  0,...
        'fig_dir',          [],...
        'um_per_pix',       1.25,...
        'xls_filename',     'auto_stats.xls',...
        'aam_thresh',       -2e4,...
        'plot_rejected',    1,...
        'plot_r', 2,...
        'plot_c', 1);    
                
end
%%% selected_features = {...
%         'mean_orientation_entropy',...
%     'mean_connected_orientation_dispersion',...
%     'mean_weighted_orientation_dispersion',...
%     'dispersion_base_orientation',...
%     'dispersion_connected_orientation',...
%     'dispersion_weighted_orientation'};

selected_features = {...
    'vessel_density1',...
    'vessel_density2',...
    'mean_inter_capillary_distance',...
    'std_inter_capillary_distance',...
    'median_inter_capillary_distance'};

auto_stats = analyse_qseries_apex_measures([],...auto_stats,...
    'image_names',          im_names,...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'vessel_centre_dir',    'vessel_centres\',...
    'selected_dir',         'apex_maps\set12g_half_296655\island_maxima\selected_apexes',...
    'metrics_dir',          'apex_maps\set12g_half_296655\island_maxima\apex_metrics',...
    'do_xls',           0,...
    'do_make_stats',    1,...
    'do_image_plots',   1, ...
    'selected_features', selected_features,... %Do them all for now
    'do_people_plots',  0,...
    'fig_dir', [],...
    'xls_filename', []);
%%

example_names = {'c02V6LD4X3LrgMosaic','d136V1LD4X3LrgMosaic','p04V1LD4X3LrgMosaic','u03V1LD4X3LrgMosaic'};

for i_im = 1:4
    im_name = [example_names{i_im} '.png'];

    n = rgb2gray(imread(['P:\isbe\conferences_and_symposia\2014\acR2104\figs\' im_name]));
    n2 = imread(['C:\isbe\nailfold\data\2_year_study\anonymous_png\' im_name]);

    patch = n(1:101, 1:101);
    corr_map = mb_normxcorr2(patch, n2);

    [~, i] = max(corr_map (:));
    [sr, sc] = ind2sub(size(corr_map ),i);
    
    plot_view = [sc-50 sc+1549 sr-50 sr+349]/2;
    plot_caxis = [min(n(:)) max(n(:))];

    display_automated_markup(... % non-strict mode
        'image_names',         example_names(i_im),...
        'data_dir',             [nailfoldroot 'data/2_year_study/'],...
        'image_dir',            'images',...
        'vessel_centre_dir',    'vessel_centres',...
        'metrics_dir',          'apex_maps\set12g_half_296655\miccai_maxima\apex_metrics',...
        'candidates_dir',       'apex_maps\set12g_half_296655\miccai_maxima\rescores',...
        'selected_dir',         'apex_maps\set12g_half_296655\miccai_maxima\selected_apexes',...
        'selected_features', [],...
        'fig_dir',          [],...
        'um_per_pix',       1.25,...
        'aam_thresh',       inf,...
        'plot_nondistal',    0,...
        'plot_rejected',    0,...
        'plot_r', 1,...
        'plot_c', 1,...
        'plot_view', plot_view,...
        'plot_caxis', plot_caxis); 
    
    title([]);
    axis off;
    set(gcf, 'windowStyle', 'normal', 'position', [1 1 1280 320]);
    exportfig(['P:\isbe\conferences_and_symposia\2014\acR2104\figs\a' im_name]);
end   
%%
extract_apex_measures_set(...
    'image_names',          example_names(3),...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'image_dir',            'images',...
    'prob_dir',             'rf_classification/296655',...
    'ori_dir',              'rf_regression/296621',...
    'width_dir',            'rf_regression/297037',...
    'candidates_dir',       'apex_maps\set12g_half_296655\miccai_maxima\rescores',...
    'displacements_dir',    'apex_maps\set12g_half_296655\miccai_maxima\displacements',...
    'selected_dir',         'apex_maps\set12g_half_296655\miccai_maxima\selected_apexes',...
    'metrics_dir',          [],...
    'width_predictor',      [],...[nailfoldroot 'models/apex/width/rf.mat'],...
    'prob_sigma',           2,...
    'ori_sigma',            0,...
    'width_sigma',          2,...
    'all',                  0,...
    'do_aam',               0,...
    'overwrite',            0,...
    'plot', 1);
%%
people_id = find(...
    people_stats.people_ids == 1002 |...
    people_stats.people_ids == 4004 |...
    people_stats.people_ids == 2136 |...
    people_stats.people_ids == 5003);

display(people_stats.mean_weighted_width(people_id, 4, 1, 1))
%%
im_recats = people_stats.super_category;
im_recats(people_stats.category == 4) = 2;
im_recats(people_stats.super_category == 2) = 3;
im_recats(people_stats.category == 5) = 4;

grps = {'HC', 'PR', 'SSc', 'UD'};
grp_comps = {'HC,PR', 'HC,SSc', 'HC,UD', 'PR,SSc', 'PR,UD', 'SSc,UD'};
grp_cols = 'rgbm';
feature_names = {'Mean capillary width', 'Mean capillary tortuosity', 'Capillary density'};
xlabels = {'log Width (mm)', 'Tortuosity', 'Density (mm^{-1})'};

fig_dir = 'P:\isbe\conferences_and_symposia\2014\submissions\acr\figs\';

for i_f = 1:3
      
    feature = selected_features(i_f);
    feature_str = feature{1};
    feature_str(feature_str == '_') = ' ';
    feature_str(1) = feature_str(1) - 32;
    
    display(['Tests for ' feature_str]);
    
    feature_scores = people_stats.(feature{1})(:,4,1,1);    
    switch feature{1}
        case 'mean_inter_capillary_distance'
            feature_scores = 1./ feature_scores;
        case {'mean_weighted_width', 'max_mean_width'}
            feature_scores = log(feature_scores);
    end
    
    valid_scores = isnumeric(feature_scores) & isfinite(feature_scores) & ~isnan(feature_scores) & ...
        people_stats.present(:,4,1,1);
    
    [p,table,stats] = anova1(feature_scores(valid_scores), im_recats(valid_scores));
    [c,m] = multcompare(stats);
    m(:,3) = m(:,1)+1.96*m(:,2);
    m(:,2) = m(:,1)-1.96*m(:,2);
    
    figure('WindowStyle', 'normal');
    set(gca, 'fontsize', 18);
    hold all;
    title(feature_names{i_f});
    ylabel('f(x)');
    xlabel(xlabels{i_f});
    
    min_val = min(feature_scores(valid_scores));
    max_val = max(feature_scores(valid_scores));
    grid_pts = linspace(min_val, max_val, 100);

    max_y = 0;
    
    leg_text = cell(4,1);
    for i_g = 1:4
        kdist = build_1d_kernel_distribution(feature_scores(valid_scores & im_recats==i_g), grid_pts, 0);
        plot(kdist.x, kdist.D_f, grp_cols(i_g), 'linewidth', 2);
        max_y = max(max_y, max(kdist.D_f));
        leg_text{i_g} = [grps{i_g} ' (n = ' num2str(sum(im_recats==i_g)) ')'];
    end
    if i_f == 3
        legend(leg_text, 'location', 'northwest');
    end
    
    for i_g = 1:4
        plot(m(i_g,1), max_y*(1+i_g/20), [grp_cols(i_g) 'o'], 'markersize', 6,  'markerfacecolor', grp_cols(i_g));
        plot([m(i_g,2) m(i_g,3)], [max_y*(1+i_g/20) max_y*(1+i_g/20)], [grp_cols(i_g) '-'], 'linewidth', 2);
        plot([m(i_g,2) m(i_g,2)], [max_y*(1+i_g/21) max_y*(1+i_g/19)], [grp_cols(i_g) '-'], 'linewidth', 2);
        plot([m(i_g,3) m(i_g,3)], [max_y*(1+i_g/21) max_y*(1+i_g/19)], [grp_cols(i_g) '-'], 'linewidth', 2);
    end
    %exportfig([fig_dir 'f' num2str(i_f) '.png']);
    switch feature{1}
        case 'mean_inter_capillary_distance'
            %m = 1./ m;
        case {'mean_weighted_width', 'max_mean_width'}
            m = exp(m);
    end
    
    
    
    for i_g = 1:4
        display([grps{i_g} ': ' num2str(m(i_g,1),3) ' (' num2str(m(i_g,2),3) ', ' num2str(m(i_g,3),3) '), n = ' num2str(sum(valid_scores&(im_recats==i_g)))]);
    end
    disp_str = ['Sig diffs between: '];
    for i_gc = 1:6
        if c(i_gc,6) < 0.05
            disp_str = [disp_str grp_comps{i_gc} '(' num2str(c(i_gc,6),2) '); '];
        end
    end
    display(disp_str);
end
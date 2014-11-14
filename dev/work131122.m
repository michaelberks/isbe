max_num_vessels = max(im_by_im_counts.all_distal);
max_fp = max(im_by_im_counts.false_positives);
bad = im_by_im_counts.bad_mosaic;
ungradeable = im_by_im_counts.ungradeable;
%%
%----------
figure; hist(im_by_im_counts.all_distal, 0:48);
title('Histogram of the number distal vessels in each image');

figure; hist(im_by_im_counts.valid_distal, 0:48);
title('Histogram of the number distal vessels in each image marked by at least 2 observers');

figure; hist(im_by_im_counts.detected_distal, 0:48);
title('Histogram of the number distal vessels detected in each image');

%-----------
figure; hist(im_by_im_counts.all_distal(~bad), 0:48);
title('Histogram of the number distal vessels in each image');
xlabel('Mosaics detected as badly registered ignored');

figure; hist(im_by_im_counts.valid_distal(~bad), 0:48);
title('Histogram of the number distal vessels in each image marked by at least 2 observers');
xlabel('Mosaics detected as badly registered ignored');

figure; hist(im_by_im_counts.detected_distal(~bad), 0:48);
title('Histogram of the number distal vessels detected in each image');
xlabel('Mosaics detected as badly registered ignored');

%------------
figure; hist(im_by_im_counts.all_distal(~bad & ~ungradeable), 0:48);
title('Histogram of the number distal vessels in each image');
xlabel('Mosaics detected as badly registered or marked by an observer as ungradeable ignored');

figure; hist(im_by_im_counts.valid_distal(~bad & ~ungradeable), 0:48);
title('Histogram of the number distal vessels in each image marked by at least 2 observers');
xlabel('Mosaics detected as badly registered or marked by an observer as ungradeable ignored');

figure; hist(im_by_im_counts.detected_distal(~bad & ~ungradeable), 0:48);
title('Histogram of the number distal vessels detected in each image');
xlabel('Mosaics detected as badly registered or marked by an observer as ungradeable ignored');

%---------------
figure; hist(im_by_im_counts.all_distal - im_by_im_counts.detected_distal, -30:30);
title('Difference between marked and detected counts');

figure; hist(im_by_im_counts.all_distal(~bad) - im_by_im_counts.detected_distal(~bad), -30:30);
title('Difference between marked and detected counts - bad images ignored');

figure; hist(im_by_im_counts.all_distal(~bad & ~ungradeable) - im_by_im_counts.detected_distal(~bad & ~ungradeable), -30:30);
title('Difference between marked and detected counts - bad and ungradeable images ignored');

%-------------
figure; plot(prctile(im_by_im_counts.all_distal(~bad) - im_by_im_counts.detected_distal(~bad), 0:100), 0:100);
title('Cumulative distribution of differences - bad mosaics ignored');
ylabel('%');
xlabel('Difference between marked and detected counts');

figure; plot(prctile(abs(im_by_im_counts.all_distal(~bad) - im_by_im_counts.detected_distal(~bad)), 0:100), 0:100);
title('Cumulative distribution of absolute differences - bad mosaics ignored');
ylabel('%');
xlabel('Difference between marked and detected counts');

figure; plot(im_by_im_counts.all_distal(~bad),...
    im_by_im_counts.all_distal(~bad) - im_by_im_counts.detected_distal(~bad), 'rx');
title('Bland-Altman plot - bad mosaics ignored');
ylabel('%');
xlabel('Marked number of vessels');
ylabel('Difference between marked and detected counts');

figure; plot(prctile(im_by_im_counts.all_distal(~bad) ./ im_by_im_counts.detected_distal(~bad), 0:100), 0:100);
title('Cumulative distribution of ratios - bad mosaics ignored');
ylabel('%');
xlabel('Ratio of marked count / detected count');
%--------------
figure; hist(im_by_im_counts.all_distal/im_by_im_counts.detected_distal, -30:30);
title('Difference between marked and detected counts');

figure; hist(im_by_im_counts.all_distal(~bad) - im_by_im_counts.detected_distal(~bad), -30:30);
title('Difference between marked and detected counts - bad images ignored');

figure; hist(im_by_im_counts.all_distal(~bad & ~ungradeable) - im_by_im_counts.detected_distal(~bad & ~ungradeable), -30:30);
title('Difference between marked and detected counts - bad and ungradeable images ignored');
%--------------
figure; hist(im_by_im_counts.false_positives, 0:max_fp);
title('Histogram of false positives per image');

figure; hist(im_by_im_counts.false_positives(~bad), 0:max_fp);
title('Histogram of false positives per image - bad images ignored');

figure; hist(im_by_im_counts.false_positives(~bad & ~ungradeable), 0:max_fp);
title('Histogram of false positives per image - bad and ungradeable images ignored');
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');

[~, fp_idx] = sort(im_by_im_counts.false_positives, 'descend');
for i_im = fp_idx(1:10)'
    if ~bad(i_im) && ~ungradeable(i_im)
        select_vessels_from_candidates_set('start_i', i_im, 'end_i', i_im,...
            'do_fill_gaps', 0, 'do_distal_sub', 0, 'do_save', 0, 'do_plot', 1, 'do_width', 0,...
            'data_dir', 'C:\isbe\nailfold\data\rsa_study\test_half\',...
            'candidates_dir', 'apex_maps\frog\local_maxima\',...
            'upper_ydist', -70,...
            'lower_ydist', 45);
        load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(i_im).name(1:6) '_apex_clusters.mat']);
        plot_apex_clusters(vessels, gcf, 2);
        legend off
        
        fp = im_by_im_counts.false_positives(i_im);
        td = im_by_im_counts.detected_distal(i_im);
        tm = im_by_im_counts.all_distal(i_im);
        
        xlabel(['Image ' num2str(i_im) ': ' num2str(tm) ' marked vessels, ' num2str(td) ' detected vessels, ' num2str(fp) ' false positives']);
        
    end
end
%%
[~, missing_idx] = sort(im_by_im_counts.all_distal - im_by_im_counts.detected_distal, 'descend');
for i_im = missing_idx(1:10)'
    if ~bad(i_im) && ~ungradeable(i_im)
        select_vessels_from_candidates_batch('start_i', i_im, 'end_i', i_im, 'do_fill_gaps', 0, 'do_save', 0, 'do_plot', 1, 'do_width', 0);
        load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(i_im).name(1:6) '_apex_clusters.mat']);
        plot_apex_clusters(vessels, gcf);
        legend off
        
        fp = im_by_im_counts.false_positives(i_im);
        td = im_by_im_counts.detected_distal(i_im);
        tm = im_by_im_counts.all_distal(i_im);
        
        xlabel([num2str(tm) ' marked vessels, ' num2str(td) ' detected vessels, ' num2str(fp) ' false positives']);
        
    end
end
%%
[~, missing_idx] = sort(im_by_im_counts.valid_distal - im_by_im_counts.detected_distal, 'descend');
for i_im = missing_idx(1:6)'
    if ~bad(i_im) && ~ungradeable(i_im)
        select_vessels_from_candidates_batch('start_i', i_im, 'end_i', i_im, 'do_fill_gaps', 0, 'do_save', 0, 'do_plot', 1, 'do_width', 0);
        load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(i_im).name(1:6) '_apex_clusters.mat']);
        plot_apex_clusters(vessels, gcf);
        legend off
        
        fp = im_by_im_counts.false_positives(i_im);
        td = im_by_im_counts.detected_distal(i_im);
        tm = im_by_im_counts.valid_distal(i_im);
        
        xlabel([num2str(tm) ' (double) marked vessels, ' num2str(td) ' detected vessels, ' num2str(fp) ' false positives']);
        
    end
end
%%
[~, missing_idx] = sort(im_by_im_counts.all_distal - im_by_im_counts.detected_distal);
for i_im = missing_idx(1:10)'
    if ~bad(i_im) && ~ungradeable(i_im)
        select_vessels_from_candidates_batch('start_i', i_im, 'end_i', i_im, 'do_fill_gaps', 0, 'do_save', 0, 'do_plot', 1, 'do_width', 0);
        load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(i_im).name(1:6) '_apex_clusters.mat']);
        plot_apex_clusters(vessels, gcf);
        legend off
        
        fp = im_by_im_counts.false_positives(i_im);
        td = im_by_im_counts.detected_distal(i_im);
        tm = im_by_im_counts.all_distal(i_im);
        
        xlabel([num2str(tm) ' marked vessels, ' num2str(td) ' detected vessels, ' num2str(fp) ' false positives']);
        
    end
end
%%
%%
[~, best_idx] = sort(im_by_im_counts.detected_distal - im_by_im_counts.false_positives, 'descend');
for i_im = best_idx(1:10)'
    if ~bad(i_im) && ~ungradeable(i_im)
        select_vessels_from_candidates_set('start_i', i_im, 'end_i', i_im,...
            'do_fill_gaps', 0, 'do_distal_sub', 0, 'do_save', 0, 'do_plot', 1, 'do_width', 1,...
            'data_dir', 'C:\isbe\nailfold\data\rsa_study\test_half\',...
            'candidates_dir', 'apex_maps\frog\local_maxima\',...
            'upper_ydist', -70,...
            'lower_ydist', 45);
        load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(i_im).name(1:6) '_apex_clusters.mat']);
        plot_apex_clusters(vessels, gcf, 2);
        legend off
        
        fp = im_by_im_counts.false_positives(i_im);
        td = im_by_im_counts.detected_distal(i_im);
        tm = im_by_im_counts.all_distal(i_im);
        
        xlabel([num2str(tm) ' marked vessels, ' num2str(td) ' detected vessels, ' num2str(fp) ' false positives']);
        
    end
end
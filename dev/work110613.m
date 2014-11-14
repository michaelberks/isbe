num_bins = 30;
[g2] = k_maps_hist('g2d', 'g2d', 'num_bins', num_bins);
[rfp] = k_maps_hist('rf_prob', 'rf_thin', 'line_stem',...
    'C:\isbe\asymmetry_project\data\orientation_maps\', 'thresh', 0.35,...
    'num_bins', num_bins, 'do_vm', 1);
%[rft] = k_maps_hist('rf_thin', [], 'thresh', 0.35, 'num_bins', num_bins);
%
rf = rfp;
%
bins = linspace(-pi/2, pi/2, 30);
%%
figure;
subplot(2,2,1);
bar(bins, mean(rf.n_hist_line(:,:,1)));
title('RF abnormals (line pixels)'); %axis([-pi/2 pi/2 0 1000]);
subplot(2,2,2);
bar(bins, mean(rf.n_hist_line(:,:,2)));
title('RF normals (line pixels)'); %axis([-pi/2 pi/2 0 1000]);
subplot(2,2,3);
bar(bins, mean(g2.n_hist_line(:,:,1)));
title('G2 abnormals (line pixels)'); %axis([-pi/2 pi/2 0 1000]);
subplot(2,2,4);
bar(bins, mean(g2.n_hist_line(:,:,2)));
title('G2 normals (line pixels)'); %axis([-pi/2 pi/2 0 1000]);

figure;
subplot(2,2,1);
bar(bins, mean(rf.n_hist_all(:,:,1)));
title('RF abnormals (all pixels)'); %axis([-pi/2 pi/2 0 2500]);
subplot(2,2,2);
bar(bins, mean(rf.n_hist_all(:,:,2)));
title('RF normals (all pixels)'); %axis([-pi/2 pi/2 0 2500]);
subplot(2,2,3);
bar(bins, mean(g2.n_hist_all(:,:,1)));
title('G2 abnormals (all pixels)'); %axis([-pi/2 pi/2 0 2500]);
subplot(2,2,4);
bar(bins, mean(g2.n_hist_all(:,:,2)));
title('G2 normals (all pixels)'); %axis([-pi/2 pi/2 0 2500]);


%
figure;
subplot(2,2,1);
bar(bins, mean(rf.n_hist_all(:,:,1)) - mean(rf.n_hist_all(:,:,2))); 
%axis([-pi/2 pi/2 0 250]);
title('RF - difference between abnormal and normal (all pixels)');

subplot(2,2,2); 
bar(bins, mean(g2.n_hist_all(:,:,1)) - mean(g2.n_hist_all(:,:,2))); 
%axis([-pi/2 pi/2 0 250]);
title('G2 - difference between abnormal and normal (all pixels)');

subplot(2,2,3); 
bar(bins, mean(rf.n_hist_line(:,:,1)) - mean(rf.n_hist_line(:,:,2))); 
%axis([-pi/2 pi/2 0 250]);
title('RF - difference between abnormal and normal (line pixels)');

subplot(2,2,4); 
bar(bins, mean(g2.n_hist_line(:,:,1)) - mean(g2.n_hist_line(:,:,2))); 
%axis([-pi/2 pi/2 0 250]);
title('G2 - difference between abnormal and normal (line pixels)');
%
figure; 
subplot(2,2,1);
bar(bins, mean(rf.n_hist_all(:,:,1)) - mean(g2.n_hist_all(:,:,1))); 
%axis([-pi/2 pi/2 0 250]);
title('Difference between RF and G2 abnormals (all pixels)');

subplot(2,2,2);
bar(bins, mean(rf.n_hist_all(:,:,2)) - mean(g2.n_hist_all(:,:,2))); 
%axis([-pi/2 pi/2 0 250]);
title('Difference between RF and G2 normals (all pixels)');

subplot(2,2,3);
bar(bins, mean(rf.n_hist_line(:,:,1)) - mean(g2.n_hist_line(:,:,1))); 
%axis([-pi/2 pi/2 -200 300]);
title('Difference between RF and G2 abnormals (line pixels)');

subplot(2,2,4);
bar(bins, mean(rf.n_hist_line(:,:,2)) - mean(g2.n_hist_line(:,:,2))); 
%axis([-pi/2 pi/2 -200 300]);
title('Difference between RF and G2 normals (line pixels)');
%%
%
figure; hold all;
plot(sort(rf.f1_all(:,1)), (0:145) / 145, 'x-');
plot(sort(rf.f1_all(:,2)), (0:145) / 145, 'x-');
plot(sort(g2.f1_all(:,1)), (0:145) / 145, 'x-');
plot(sort(g2.f1_all(:,2)), (0:145) / 145, 'x-');
legend({...
    'RF abnormals (all pixels)',...
    'RF normals (all pixels)',...
    'G2 abnormals (all pixels)',...
    'G2 normals (all pixels)'});
title('CDF of f1 scores');
%
figure; hold all;
plot(sort(rf.f1_line(:,1)), (0:145) / 145, 'x-');
plot(sort(rf.f1_line(:,2)), (0:145) / 145, 'x-');
plot(sort(g2.f1_line(:,1)), (0:145) / 145, 'x-');
plot(sort(g2.f1_line(:,2)), (0:145) / 145, 'x-');
legend({...
    'RF abnormals (line pixels)',...
    'RF normals (line pixels)',...
    'G2 abnormals (line pixels)',...
    'G2 normals (line pixels)'});
title('CDF of f1 scores');
%%
[g2_roc_line g2_auc_line] = calculate_roc_curve(...
    g2.f1_line(:), [true(146,1); false(146,1)], sort(g2.f1_line(:)));
[rf_roc_line rf_auc_line] = calculate_roc_curve(...
    rf.f1_line(:), [true(146,1); false(146,1)], sort(rf.f1_line(:)));
[g2_roc_all g2_auc_all] = calculate_roc_curve(...
    g2.f1_all(:), [true(146,1); false(146,1)], sort(g2.f1_all(:)));
[rf_roc_all rf_auc_all] = calculate_roc_curve(...
    rf.f1_all(:), [true(146,1); false(146,1)], sort(rf.f1_all(:)));

figure; hold all;
plot(rf_roc_line(:,1), rf_roc_line(:,2), '-');
plot(rf_roc_all(:,1), rf_roc_all(:,2), '-');
plot(g2_roc_line(:,1), g2_roc_line(:,2), '-');
plot(g2_roc_all(:,1), g2_roc_all(:,2), '-');
axis([0 1 0 1]); axis equal;
legend({...
    'RF (line pixels)',...
    'RF (all pixels)',...
    'G2(line pixels)',...
    'G2 (all pixels)'});
title('ROC curve of f1 scores');
%
[g2_roc_line g2_auc_line_n] = calculate_roc_curve(...
    g2.n_line(:), [true(146,1); false(146,1)], sort(g2.n_line(:)));
[rf_roc_line rf_auc_line_n] = calculate_roc_curve(...
    rf.n_line(:), [true(146,1); false(146,1)], sort(rf.n_line(:)));
[g2_roc_all g2_auc_all_n] = calculate_roc_curve(...
    g2.n_all(:), [true(146,1); false(146,1)], sort(g2.n_all(:)));
[rf_roc_all rf_auc_all_n] = calculate_roc_curve(...
    rf.n_all(:), [true(146,1); false(146,1)], sort(rf.n_all(:)));

figure; hold all;
plot(rf_roc_line(:,1), rf_roc_line(:,2), '-');
plot(rf_roc_all(:,1), rf_roc_all(:,2), '-');
plot(g2_roc_line(:,1), g2_roc_line(:,2), '-');
plot(g2_roc_all(:,1), g2_roc_all(:,2), '-');
axis([0 1 0 1]); axis equal;
legend({...
    'RF (line pixels)',...
    'RF (all pixels)',...
    'G2(line pixels)',...
    'G2 (all pixels)'});
title('ROC curve of n counts');
%%
figure; hold all;
plot(rf.n_line(:,1), rf.f1_line(:,1), 'x');
plot(g2.n_line(:,1), g2.f1_line(:,1), 'x');

plot(rf.n_line(:,2), rf.f1_line(:,2), 'o');
plot(g2.n_line(:,2), g2.f1_line(:,2), 'o');

legend({...
    'RF abnormal lines',...
    'G2 abnormal lines',...
    'RF normal lines',...
    'G2 normal lines'});
%%
rf.n_line_1 = [sum(rf.n_hist_line(:,11:20,1),2) sum(rf.n_hist_line(:,11:20,2),2)];
rf.N_line_1 = [sum(rf.N_hist_line(:,11:20,1),2) sum(rf.N_hist_line(:,11:20,2),2)];
g2.n_line_1 = [sum(g2.n_hist_line(:,11:20,1),2) sum(g2.n_hist_line(:,11:20,2),2)];
g2.N_line_1 = [sum(g2.N_hist_line(:,11:20,1),2) sum(g2.N_hist_line(:,11:20,2),2)];

rf.f1_line_1 = (rf.n_line_1 - (rf.N_line_1 .* rf.p_line)) ./...
    sqrt(rf.N_line_1 .* rf.p_line .* (1 - rf.p_line));
g2.f1_line_1 = (g2.n_line_1 - (g2.N_line_1 .* g2.p_line)) ./...
    sqrt(g2.N_line_1 .* g2.p_line .* (1 - g2.p_line));

[g2_roc_line_1 g2_auc_line_1] = calculate_roc_curve(...
    g2.f1_line_1(:), [true(146,1); false(146,1)], sort(g2.f1_line_1(:)));
[rf_roc_line_1 rf_auc_line_1] = calculate_roc_curve(...
    rf.f1_line_1(:), [true(146,1); false(146,1)], sort(rf.f1_line_1(:)));

figure; hold all;
plot(rf_roc_line_1(:,1), rf_roc_line_1(:,2), '-');
plot(g2_roc_line_1(:,1), g2_roc_line_1(:,2), '-');

axis([0 1 0 1]); axis equal;
legend({...
    'RF (line pixels)',...
    'G2(line pixels)'});
title('ROC curve of f1 scores');
%%
rf.n_line_2 = [sum(rf.n_hist_line(:,[1:10 21:30],1),2) sum(rf.n_hist_line(:,[1:10 21:30],2),2)];
rf.N_line_2 = [sum(rf.N_hist_line(:,[1:10 21:30],1),2) sum(rf.N_hist_line(:,[1:10 21:30],2),2)];
g2.n_line_2 = [sum(g2.n_hist_line(:,[1:10 21:30],1),2) sum(g2.n_hist_line(:,[1:10 21:30],2),2)];
g2.N_line_2 = [sum(g2.N_hist_line(:,[1:10 21:30],1),2) sum(g2.N_hist_line(:,[1:10 21:30],2),2)];

rf.f1_line_2 = (rf.n_line_2 - (rf.N_line_2 .* rf.p_line)) ./...
    sqrt(rf.N_line_2 .* rf.p_line .* (1 - rf.p_line));
g2.f1_line_2 = (g2.n_line_2 - (g2.N_line_2 .* g2.p_line)) ./...
    sqrt(g2.N_line_2 .* g2.p_line .* (1 - g2.p_line));

[g2_roc_line_2 g2_auc_line_2] = calculate_roc_curve(...
    g2.f1_line_2(:), [true(146,1); false(146,1)], sort(g2.f1_line_2(:)));
[rf_roc_line_2 rf_auc_line_2] = calculate_roc_curve(...
    rf.f1_line_2(:), [true(146,1); false(146,1)], sort(rf.f1_line_2(:)));

figure; hold all;
plot(rf_roc_line_2(:,1), rf_roc_line_2(:,2), '-');
plot(g2_roc_line_2(:,1), g2_roc_line_2(:,2), '-');

axis([0 1 0 1]); axis equal;
legend({...
    'RF (line pixels)',...
    'G2(line pixels)'});
title('ROC curve of f1 scores');
%%
for ii = 1:10
    figure; 
    subplot(2,1,1); plot(g2.n_hist_line(ii,:,1)');
    subplot(2,1,2); plot(g2.n_hist_line(ii,:,2)');
end
%%
[g2_roc_line g2_auc_line] = calculate_roc_curve(...
    g2.f1_line(:), [true(146,1); false(146,1)], sort(g2.f1_line(:)));

[g2_roc_vm g2_auc_vm] = calculate_roc_curve(...
    g2.f1_line_vm(:), [true(146,1); false(146,1)], sort(g2.f1_line_vm(:)));


figure; hold all;

plot(g2_roc_line(:,1), g2_roc_line(:,2), '-');
plot(g2_roc_vm(:,1), g2_roc_vm(:,2), '-');
axis([0 1 0 1]); axis equal;
legend({...
    'G2',...
    'G2 vm'});
title('ROC curve of f1 scores');

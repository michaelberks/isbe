load C:\isbe\dev\files\small_files.mat
load C:\isbe\dev\files\large_files.mat
load C:\isbe\dev\files\r_files50.mat
load C:\isbe\dev\files\r_files51.mat
load C:\isbe\dev\files\u_files.mat
%%
%
% Get loo errors for scale models
%
[er_small.com er_small.ind] = model_errors_loo(small_files, [], 'C:\isbe\dev\models\loo_small\model');
save C:\isbe\dev\scale\er_small er_small;
[er_large.com er_large.ind] = model_errors_loo(large_files, [], 'C:\isbe\dev\models\loo_large\model');
save C:\isbe\dev\scale\er_large er_large;
[er_r51.com er_r51.ind] = model_errors_loo(r_files51, [], 'C:\isbe\dev\models\loo_r51\model');
save C:\isbe\dev\scale\er_r51 er_r51;
[er_r50.com er_r50.ind] = model_errors_loo(r_files50, [], 'C:\isbe\dev\models\loo_r50\model');
save C:\isbe\dev\scale\er_r50 er_r50;
[er_u.com er_u.ind] = model_errors_loo(u_files1, [], 'C:\isbe\dev\models\loo_u1\model');
save C:\isbe\dev\scale\er_u er_u;
%%
load C:\isbe\dev\scale\er_small.mat
load C:\isbe\dev\scale\er_large.mat
load C:\isbe\dev\scale\er_r50.mat
load C:\isbe\dev\scale\er_r51
load C:\isbe\dev\scale\er_u.mat
%%
%
% Compute error matrices comparing error for each mass
ers_combined = nan;
ers_combined = ers_combined(ones(179,3));
ers_shape = ers_combined;
ers_tex = ers_combined;
ers_scale = ers_combined;
ers_indie = ers_combined;

ers_indie(idx_small,1) = er_small.ind.total;
ers_shape(idx_small,1) = er_small.ind.shape;
ers_tex(idx_small,1) = er_small.ind.tex;
ers_scale(idx_small,1) = er_small.com.scale;
ers_combined(idx_small,1) = er_small.com.total;

ers_indie(idx_large, 1) = er_large.ind.total;
ers_shape(idx_large, 1) = er_large.ind.shape;
ers_tex(idx_large, 1) = er_large.ind.tex;
ers_scale(idx_large, 1) = er_large.com.scale;
ers_combined(idx_large, 1) = er_large.com.total;

ers_indie(idx_r50,2) = er_r50.ind.total;
ers_shape(idx_r50,2) = er_r50.ind.shape;
ers_tex(idx_r50,2) = er_r50.ind.tex;
ers_scale(idx_r50,2) = er_r50.com.scale;
ers_combined(idx_r50,2) = er_r50.com.total;

ers_indie(idx_r51,2) = er_r51.ind.total;
ers_shape(idx_r51,2) = er_r51.ind.shape;
ers_tex(idx_r51,2) = er_r51.ind.tex;
ers_scale(idx_r51,2) = er_r51.com.scale;
ers_combined(idx_r51,2) = er_r51.com.total;

ers_indie(idx_u1, 3) = er_u.ind.total;
ers_shape(idx_u1, 3) = er_u.ind.shape;
ers_tex(idx_u1, 3) = er_u.ind.tex;
ers_scale(idx_u1, 3) = er_u.com.scale;
ers_combined(idx_u1, 3) = er_u.com.total;

ers_indie(setdiff(1:179,idx_u1), :) = [];
ers_shape(setdiff(1:179,idx_u1), :) = [];
ers_tex(setdiff(1:179,idx_u1), :) = [];
ers_scale(setdiff(1:179,idx_u1), :) = [];
ers_combined(setdiff(1:179,idx_u1), :) = [];
%%
% for ii = 1:50
%     ers_indie(idx_small(ii),1) = er_small.ind.total(ii);
%     ers_shape(idx_small(ii),1) = er_small.ind.shape(ii);
%     ers_tex(idx_small(ii),1) = er_small.ind.tex(ii);
%     ers_scale(idx_small(ii),1) = er_small.com.scale(ii);
%     ers_combined(idx_small(ii),1) = er_small.com.total(ii);
% 
%     ers_indie(idx_r50(ii),3) = er_r50.ind.total(ii);
%     ers_shape(idx_r50(ii),3) = er_r50.ind.shape(ii);
%     ers_tex(idx_r50(ii),3) = er_r50.ind.tex(ii);
%     ers_scale(idx_r50(ii),3) = er_r50.com.scale(ii);
%     ers_combined(idx_r50(ii),3) = er_r50.com.total(ii);
% end
% for ii = 1:51
%     ers_indie(idx_large(ii), 2) = er_large.ind.total(ii);
%     ers_shape(idx_large(ii), 2) = er_large.ind.shape(ii);
%     ers_tex(idx_large(ii), 2) = er_large.ind.tex(ii);
%     ers_scale(idx_large(ii), 2) = er_large.com.scale(ii);
%     ers_combined(idx_large(ii), 2) = er_large.com.total(ii);
%     
%     ers_indie(idx_r51(ii),4) = er_r51.ind.total(ii);
%     ers_shape(idx_r51(ii),4) = er_r51.ind.shape(ii);
%     ers_tex(idx_r51(ii),4) = er_r51.ind.tex(ii);
%     ers_scale(idx_r51(ii),4) = er_r51.com.scale(ii);
%     ers_combined(idx_r51(ii),4) = er_r51.com.total(ii);
% end
% for ii = 1:101
%     ers_indie(idx_u1(ii), 5) = er_u.ind.total(ii);
%     ers_shape(idx_u1(ii), 5) = er_u.ind.shape(ii);
%     ers_tex(idx_u1(ii), 5) = er_u.ind.tex(ii);
%     ers_scale(idx_u1(ii), 5) = er_u.com.scale(ii);
%     ers_combined(idx_u1(ii), 5) = er_u.com.total(ii);
% end

ers_indie = sortrows(ers_indie, 5);
ers_indie(102:179,:) = [];
ers_shape = sortrows(ers_shape, 5);
ers_shape(102:179,:) = [];
ers_tex = sortrows(ers_tex, 5);
ers_tex(102:179,:) = [];
ers_scale = sortrows(ers_scale, 5);
ers_scale(102:179,:) = [];
ers_combined = sortrows(ers_combined, 5);
ers_combined(102:179,:) = [];

ers_indie2 = ers_indie(~isnan(ers_indie(:,1)), [1 5]);
ers_shape2 = ers_shape(~isnan(ers_shape(:,1)), [1 5]);
ers_tex2 = ers_tex(~isnan(ers_tex(:,1)), [1 5]);
ers_scale2 = ers_scale(~isnan(ers_scale(:,1)), [1 5]);
ers_combined2 = ers_combined(~isnan(ers_combined(:,1)), [1 5]);

ers_indie3 = ers_indie(~isnan(ers_indie(:,2)), [2 5]);
ers_shape3 = ers_shape(~isnan(ers_shape(:,2)), [2 5]);
ers_tex3 = ers_tex(~isnan(ers_tex(:,2)), [2 5]);
ers_scale3 = ers_scale(~isnan(ers_scale(:,2)), [2 5]);
ers_combined3 = ers_combined(~isnan(ers_combined(:,2)), [2 5]);

ers_indie4 = ers_indie(~isnan(ers_indie(:,3)), [3 5]);
ers_shape4 = ers_shape(~isnan(ers_shape(:,3)), [3 5]);
ers_tex4 = ers_tex(~isnan(ers_tex(:,3)), [3 5]);
ers_scale4 = ers_scale(~isnan(ers_scale(:,3)), [3 5]);
ers_combined4 = ers_combined(~isnan(ers_combined(:,3)), [3 5]);

ers_indie5 = ers_indie(~isnan(ers_indie(:,4)), [4 5]);
ers_shape5 = ers_shape(~isnan(ers_shape(:,4)), [4 5]);
ers_tex5 = ers_tex(~isnan(ers_tex(:,4)), [4 5]);
ers_scale5 = ers_scale(~isnan(ers_scale(:,4)), [4 5]);
ers_combined5 = ers_combined(~isnan(ers_combined(:,4)), [4 5]);

%save the error matrices
save C:\isbe\dev\scale\errors ers*
%%
% Plot and save errors calculated above
figure('windowstyle', 'docked'); hold on; plot(ers_indie, 'x');
legend({'small model', 'large model', 'random 50 model', 'random 51 model', 'full model'}); 
title('Total individual model errors');
saveas(gcf, 'C:\isbe\dev\scale\figures\1.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_shape, 'x');
legend({'small model', 'large model', 'random 50 model', 'random 51 model', 'full model'}); 
title('Shape individual model errors');
saveas(gcf, 'C:\isbe\dev\scale\figures\2.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_tex, 'x');
legend({'small model', 'large model', 'random 50 model', 'random 51 model', 'full model'}); 
title('Texture individual model errors');
saveas(gcf, 'C:\isbe\dev\scale\figures\3.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_scale, 'x');
legend({'small model', 'large model', 'random 50 model', 'random 51 model', 'full model'}); 
title('Scale combined model errors');
saveas(gcf, 'C:\isbe\dev\scale\figures\4.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_combined, 'x');
legend({'small model', 'large model', 'random 50 model', 'random 51 model', 'full model'}); 
title('Total combined model errors');
saveas(gcf, 'C:\isbe\dev\scale\figures\5.fig');

figure('windowstyle', 'docked'); hold on; plot(ers_indie2, 'x');
legend({'small model', 'full model'}); 
title('Total individual model errors - small');
saveas(gcf, 'C:\isbe\dev\scale\figures\6.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_shape2, 'x');
legend({'small model', 'full model'}); 
title('Shape individual model errors - small');
saveas(gcf, 'C:\isbe\dev\scale\figures\7.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_tex2, 'x');
legend({'small model','full model'}); 
title('Texture individual model errors - small');
saveas(gcf, 'C:\isbe\dev\scale\figures\8.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_scale2, 'x');
legend({'small model', 'full model'}); 
title('Scale combined model errors - small');
saveas(gcf, 'C:\isbe\dev\scale\figures\9.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_combined2, 'x');
legend({'small model', 'full model'}); 
title('Total combined model errors - small');
saveas(gcf, 'C:\isbe\dev\scale\figures\10.fig');

figure('windowstyle', 'docked'); hold on; plot(ers_indie3, 'x');
legend({'large model', 'full model'}); 
title('Total individual model errors - large');
saveas(gcf, 'C:\isbe\dev\scale\figures\11.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_shape3, 'x');
legend({'large model', 'full model'}); 
title('Shape individual model errors - large');
saveas(gcf, 'C:\isbe\dev\scale\figures\12.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_tex3, 'x');
legend({'large model', 'full model'}); 
title('Texture individual model errors - large');
saveas(gcf, 'C:\isbe\dev\scale\figures\13.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_scale3, 'x');
legend({'large model', 'full model'}); 
title('Scale combined model errors - large');
saveas(gcf, 'C:\isbe\dev\scale\figures\14.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_combined3, 'x');
legend({'large model', 'full model'}); 
title('Total combined model errors - large');
saveas(gcf, 'C:\isbe\dev\scale\figures\15.fig');

figure('windowstyle', 'docked'); hold on; plot(ers_indie4, 'x');
legend({'random 50 model', 'full model'}); 
title('Total individual model errors - r50');
saveas(gcf, 'C:\isbe\dev\scale\figures\16.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_shape4, 'x');
legend({'random 50 model', 'full model'}); 
title('Shape individual model errors - r50');
saveas(gcf, 'C:\isbe\dev\scale\figures\17.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_tex4, 'x');
legend({'random 50 model','full model'}); 
title('Texture individual model errors - r50');
saveas(gcf, 'C:\isbe\dev\scale\figures\18.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_scale4, 'x');
legend({'random 50 model', 'full model'}); 
title('Scale combined model errors - r50');
saveas(gcf, 'C:\isbe\dev\scale\figures\19.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_combined4, 'x');
legend({'random 50 model','full model'}); 
title('Total combined model errors - r50');
saveas(gcf, 'C:\isbe\dev\scale\figures\20.fig');

figure('windowstyle', 'docked'); hold on; plot(ers_indie5, 'x');
legend({'random 51 model', 'full model'}); 
title('Total individual model errors - r51');
saveas(gcf, 'C:\isbe\dev\scale\figures\21.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_shape5, 'x');
legend({'random 51 model', 'full model'}); 
title('Shape individual model errors - r51');
saveas(gcf, 'C:\isbe\dev\scale\figures\22.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_tex5, 'x');
legend({'random 51 model', 'full model'}); 
title('Texture individual model errors - r51');
saveas(gcf, 'C:\isbe\dev\scale\figures\23.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_scale5, 'x');
legend({'random 51 model', 'full model'}); 
title('Scale combined model errors - r51');
saveas(gcf, 'C:\isbe\dev\scale\figures\24.fig');
figure('windowstyle', 'docked'); hold on; plot(ers_combined5, 'x');
legend({'random 51 model', 'full model'}); 
title('Total combined model errors - r51');
saveas(gcf, 'C:\isbe\dev\scale\figures\25.fig');
%%

ers_combined_a = zeros(101, 3);
ers_shape_a = zeros(101, 3);
ers_scale_a = zeros(101, 3);
ers_tex_a = zeros(101, 3);
ers_indie_a = zeros(101, 3);

ind_small = ~isnan(ers_combined(:,1));
ind_large = ~isnan(ers_combined(:,2));
ind_r50 = ~isnan(ers_combined(:,3));
ind_r51 = ~isnan(ers_combined(:,4));

ers_combined_a(ind_small, 1) = ers_combined(ind_small, 1);
ers_shape_a(ind_small, 1) = ers_shape(ind_small, 1);
ers_scale_a(ind_small, 1) = ers_scale(ind_small, 1);
ers_tex_a(ind_small, 1) = ers_tex(ind_small, 1);
ers_indie_a(ind_small, 1) = ers_indie(ind_small, 1);

ers_combined_a(~isnan(ers_combined(:,2)), 1) = ...
    ers_combined(~isnan(ers_combined(:,2)), 2);
ers_shape_a(~isnan(ers_shape(:,2)), 1) = ...
    ers_shape(~isnan(ers_shape(:,2)), 2);
ers_scale_a(~isnan(ers_scale(:,2)), 1) = ...
    ers_scale(~isnan(ers_scale(:,2)), 2);
ers_tex_a(~isnan(ers_tex(:,2)), 1) = ...
    ers_tex(~isnan(ers_tex(:,2)), 2);
ers_indie_a(~isnan(ers_indie(:,2)), 1) = ...
    ers_indie(~isnan(ers_indie(:,2)), 2);

ers_combined_a(~isnan(ers_combined(:,3)), 2) = ...
    ers_combined(~isnan(ers_combined(:,3)), 3);
ers_shape_a(~isnan(ers_shape(:,3)), 2) = ...
    ers_shape(~isnan(ers_shape(:,3)), 3);
ers_scale_a(~isnan(ers_scale(:,3)), 2) = ...
    ers_scale(~isnan(ers_scale(:,3)), 3);
ers_tex_a(~isnan(ers_tex(:,3)), 2) = ...
    ers_tex(~isnan(ers_tex(:,3)), 3);
ers_indie_a(~isnan(ers_indie(:,3)), 2) = ...
    ers_indie(~isnan(ers_indie(:,3)), 3);

ers_combined_a(~isnan(ers_combined(:,4)), 2) = ...
    ers_combined(~isnan(ers_combined(:,4)), 4);
ers_shape_a(~isnan(ers_shape(:,4)), 2) = ...
    ers_shape(~isnan(ers_shape(:,4)), 4);
ers_scale_a(~isnan(ers_scale(:,4)), 2) = ...
    ers_scale(~isnan(ers_scale(:,4)), 4);
ers_tex_a(~isnan(ers_tex(:,4)), 2) = ...
    ers_tex(~isnan(ers_tex(:,4)), 4);
ers_indie_a(~isnan(ers_indie(:,4)), 2) = ...
    ers_indie(~isnan(ers_indie(:,4)), 4);

ers_combined_a(~isnan(ers_combined(:,5)), 3) = ...
    ers_combined(~isnan(ers_combined(:,5)), 5);
ers_shape_a(~isnan(ers_shape(:,5)), 3) = ...
    ers_shape(~isnan(ers_shape(:,5)), 5);
ers_scale_a(~isnan(ers_scale(:,5)), 3) = ...
    ers_scale(~isnan(ers_scale(:,5)), 5);
ers_tex_a(~isnan(ers_tex(:,5)), 3) = ...
    ers_tex(~isnan(ers_tex(:,5)), 5);
ers_indie_a(~isnan(ers_indie(:,5)), 3) = ...
    ers_indie(~isnan(ers_indie(:,5)), 5);
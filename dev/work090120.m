%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Displaying dual-tree feature vector data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1) 2D Histograms across each pair dimensions
%
% Useful for comparing fit of GMM to data as in code work090109. Below we
% simply look at phase-phase pairs to compare with 2
counts12 = hist3(data(:,2*[1 2]), [100, 100]);
counts13 = hist3(data(:,2*[1 3]), [100, 100]);
counts14 = hist3(data(:,2*[1 4]), [100, 100]);

counts23 = hist3(data(:,2*[2 3]), [100, 100]);
counts24 = hist3(data(:,2*[2 4]), [100, 100]);

counts34 = hist3(data(:,2*[3 4]), [100, 100]);
%
figure; imagesc(counts12); axis image;
title('Un-weighted histograms of phase-phase pairs');
xlabel('Level 2 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\counts_1_2.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(counts13); axis image;
title('Un-weighted histograms of phase-phase pairs');
xlabel('Level 3 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\counts_1_3.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(counts14); axis image;
title('Un-weighted histograms of phase-phase pairs');
xlabel('Level 4 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\counts_1_4.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(counts23); axis image;
title('Un-weighted histograms of phase-phase pairs');
xlabel('Level 3 phase'); ylabel('Level 2: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\counts_2_3.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(counts24); axis image;
title('Un-weighted histograms of phase-phase pairs');
xlabel('Level 4 phase'); ylabel('Level 2: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\counts_2_4.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(counts34); axis image;
title('Un-weighted histograms of phase-phase pairs');
xlabel('Level 4 phase'); ylabel('Level 3: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\counts_3_4.eps';
% saveas(gcf, f_name, 'psc2');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2) 2D Histograms of phase v phase weighted by a function of magnitudes
%
% Picks out relations of phase for 'important' structures

w_counts12 = mb_dual_tree_weighted_hist(data, {100,100}, [1 2], [1 2]);
w_counts13 = mb_dual_tree_weighted_hist(data, {100,100}, [1 3], [1 3]);
w_counts14 = mb_dual_tree_weighted_hist(data, {100,100}, [1 4], [1 4]);

w_counts23 = mb_dual_tree_weighted_hist(data, {100,100}, [2 3], [2 3]);
w_counts24 = mb_dual_tree_weighted_hist(data, {100,100}, [2 4], [2 4]);

w_counts34 = mb_dual_tree_weighted_hist(data, {100,100}, [3 4], [3 4]);

figure; imagesc(w_counts12); axis image;
title('Magnitude weighted histograms of phase-phase pairs');
xlabel('Level 2 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);

figure; imagesc(w_counts13); axis image;
title('Magnitude weighted histograms of phase-phase pairs');
xlabel('Level 3 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);

figure; imagesc(w_counts14); axis image;
title('Magnitude weighted histograms of phase-phase pairs');
xlabel('Level 4 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);

figure; imagesc(w_counts23); axis image;
title('Magnitude weighted histograms of phase-phase pairs');
xlabel('Level 3 phase'); ylabel('Level 2: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);

figure; imagesc(w_counts24); axis image;
title('Magnitude weighted histograms of phase-phase pairs');
xlabel('Level 4 phase'); ylabel('Level 2: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);

figure; imagesc(w_counts34); axis image;
title('Magnitude weighted histograms of phase-phase pairs');
xlabel('Level 4 phase'); ylabel('Level 3: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
%%
% Lets also look at the log-scale of the weighted counts because some
% counts will be particularly high, skewing the scale
wl_counts12 = log2(w_counts12);
wl_counts13 = log2(w_counts13);
wl_counts14 = log2(w_counts14);

wl_counts23 = log2(w_counts23);
wl_counts24 = log2(w_counts24);

wl_counts34 = log2(w_counts34);
%
figure; imagesc(wl_counts12); axis image;
title('Magnitude weighted histograms of phase-phase pairs - log scaled');
xlabel('Level 2 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\weighted_log_counts_1_2.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(wl_counts13); axis image;
title('Magnitude weighted histograms of phase-phase pairs - log scaled');
xlabel('Level 3 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\weighted_log_counts_1_3.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(wl_counts14); axis image;
title('Magnitude weighted histograms of phase-phase pairs - log scaled');
xlabel('Level 4 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\weighted_log_counts_1_4.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(wl_counts23); axis image;
title('Magnitude weighted histograms of phase-phase pairs - log scaled');
xlabel('Level 3 phase'); ylabel('Level 2: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\weighted_log_counts_2_3.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(wl_counts24); axis image;
title('Magnitude weighted histograms of phase-phase pairs - log scaled');
xlabel('Level 4 phase'); ylabel('Level 2: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\weighted_log_counts_2_4.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(wl_counts34); axis image;
title('Magnitude weighted histograms of phase-phase pairs - log scaled');
xlabel('Level 4 phase'); ylabel('Level 3: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\weighted_log_counts_3_4.eps';
% saveas(gcf, f_name, 'psc2');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3)
% How about looking at the weighted counts, normalised by the unweighted
% counts (i.e. the average magnitude sum for each bin)
a_counts12 = w_counts12 ./ counts12;
a_counts12(isinf(a_counts12)) = 0;
a_counts13 = w_counts13 ./ counts13;
a_counts13(isinf(a_counts13)) = 0;
a_counts14 = w_counts14 ./ counts14;
a_counts14(isinf(a_counts14)) = 0;

a_counts23 = w_counts23 ./ counts23;
a_counts23(isinf(a_counts23)) = 0;
a_counts24 = w_counts24 ./ counts12;
a_counts24(isinf(a_counts24)) = 0;

a_counts34 = w_counts34 ./ counts34;
a_counts34(isinf(a_counts34)) = 0;
%%
figure; imagesc((a_counts12)); axis image; caxis(prctile(a_counts12(:), [0 98]));
title('Average magnitude per bin for phase-phase pairs - log scaled');
xlabel('Level 2 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\average_log_counts_1_2.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc((a_counts13)); axis image; caxis(prctile(a_counts12(:), [0 98]));
title('Average magnitude per bin for phase-phase pairs - log scaled');
xlabel('Level 3 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\average_log_counts_1_3.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc((a_counts14)); axis image; caxis(prctile(a_counts12(:), [0 98]));
title('Average magnitude per bin for phase-phase pairs - log scaled');
xlabel('Level 4 phase'); ylabel('Level 1: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\average_log_counts_1_4.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc((a_counts23)); axis image; caxis(prctile(a_counts12(:), [0 98]));
title('Average magnitude per bin for phase-phase pairs - log scaled');
xlabel('Level 3 phase'); ylabel('Level 2: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\average_log_counts_2_3.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc((a_counts24)); axis image; caxis(prctile(a_counts12(:), [0 98]));
title('Average magnitude per bin for phase-phase pairs - log scaled');
xlabel('Level 4 phase'); ylabel('Level 2: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\average_log_counts_2_4.eps';
% saveas(gcf, f_name, 'psc2');

figure; imagesc(log2(a_counts34)); axis image;
title('Average magnitude per bin for phase-phase pairs - log scaled');
xlabel('Level 4 phase'); ylabel('Level 3: phase');
set(gca, 'xtick',0:10:100, 'xticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
set(gca, 'ytick', 0:10:100, 'yticklabel', round(100*linspace(-pi/2, pi/2, 11))/100);
% f_name = 'C:\isbe\dev\background\location\figures\average_log_counts_3_4.eps';
% saveas(gcf, f_name, 'psc2');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4)
% Take a PCA transformation of the data and look pairwise at the histograms
% of parameters along each mode
%[px, lx, bx, mx] = st_pca(data, 1); clear data; pack

c1 = prctile(bx(:,1), [5 95]);
lims1 = linspace(c1(1), c1(2), 102);
c2 = prctile(bx(:,2), [5 95]);
lims2 = linspace(c2(1), c2(2), 102);
%%
lims1 = linspace(-2, 2, 102);

num_dims = size(px, 1);
for dim1 = 1:num_dims
    for dim2 = dim1+1:num_dims
        
       
        
        [bin_counts bin_centres] = hist3(bx(:,[dim1 dim2]), {lims1, lims1});
        
        figure;
        
        imagesc(bin_counts(2:end-1, 2:end-1)); axis image;
        set(gca, 'xtick',10:10:100, 'xticklabel', round(100*bin_centres{2}(10:10:end))/100);
        set(gca, 'ytick', 10:10:100, 'yticklabel', round(100*bin_centres{1}(10:10:end))/100);
        xlabel(['Mode ', num2str(dim2)]);
        ylabel(['Mode ', num2str(dim1)]);
        title('Pairwise histogram counts of parameters along the principal components');
        
        zerox = interp1(bin_centres{2}, 1:102, 0);
        zeroy = interp1(bin_centres{1}, 1:102, 0);
        hold on;
        plot(1:100, zeroy, 'm:', zerox, 1:100, 'm:');
        
        f_name = ['C:\isbe\dev\background\location\figures\pca_counts_', num2str(dim1), '_', num2str(dim2), '.eps'];
        saveas(gcf, f_name, 'psc2');
        
    
    end
end

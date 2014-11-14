cluster_list = dir('C:\isbe\dev\background\models\dual_tree\ilp_mag\*normal512*');

n = length(cluster_list);

for ii = 1:n
    models(ii).clusters = u_load(['C:\isbe\dev\background\models\dual_tree\ilp_mag\', cluster_list(ii).name]); %#ok
end

%
markers = '+x*ds^hv><';
colors = 'rgbkmcy';
%%
figure('name', 'Cluster centres for each model');

%Work through the dimensions of data
for d = 1:7
    subplot(2,4,d); hold on;
    %figure; hold on;
    %Plot each cluster centre in each model
    for m = 1:1
        for k = 1:10
            %plot d, d+1, d+2 dimensions of k-th cluster in m-th model
            dim1 = models(m).clusters.Means(k,d);
            dim2 = models(m).clusters.Means(k,d+1);
            plot(dim1, dim2, [colors(m) markers(k)], 'MarkerSize', 20-(2*m));
        end
    end
    title(['Dims ', num2str(d), ' and ', num2str(d+1)]);  
end
%%
figure('name', 'Cluster centres for each model');

%Work through the dimensions of data
for d = 1:8
    subplot(2,4,d); hold on;
    
    %Plot each cluster centre in each model
    for m = 1:n
        for k = 1:8
            %plot d, d+1, d+2 dimensions of k-th cluster in m-th model
            dim1 = models(m).clusters.Means(k,d);
            dim2 = models(m).clusters.Means(k,d+1);
            dim3 = models(m).clusters.Means(k,d+2);
            plot3(dim1, dim2, dim3, [colors(m) markers(k)]);
        end
    end
      
end
%%
%
m_cluster_list = dir('C:\isbe\dev\background\models\dual_tree\*mass_2*');

n_m = length(m_cluster_list);

for ii = 1:n_m
    m_models(ii).clusters = u_load(['C:\isbe\dev\background\models\dual_tree\', m_cluster_list(ii).name]); %#ok
end

%
markers = '+x*ds^hv';
colors = 'rgbkmcy';
%%
figure('name', 'Cluster centres for each model');

%Work through the dimensions of data
for d = 1:9
    subplot(2,4,d); hold on;
    %figure; hold on;
    %Plot each cluster centre in each model
    for m = 1:n_m
        for k = 1:8
            %plot d, d+1, d+2 dimensions of k-th cluster in m-th model
            dim1 = m_models(m).clusters.Means(k,d);
            dim2 = m_models(m).clusters.Means(k,d+1);
            plot(dim1, dim2, [colors(m) markers(k)], 'MarkerSize', 20-(2*m));
        end
    end
    title(['Dims ', num2str(d), ' and ', num2str(d+1)]);  
end
%%
%
figure;
%Work through the dimensions of data
for d = 1:6
    subplot(2,3,d); hold on;
    %figure; hold on;
    
    [bin_counts bin_centres] = hist3(data(:,[d d+1]), [100 100]);
    imagesc(bin_counts); axis image;
    set(gca, 'xtick',10:10:100, 'xticklabel', bin_centres{1}(10:10:end));
    set(gca, 'ytick', 10:10:100, 'yticklabel', bin_centres{2}(10:10:end));
    title(['Dims ', num2str(d), ' and ', num2str(d+1)]);  
end
%%
figure;

for d = 1:8
    subplot(2,4,d);
    if rem(d, 2)
        hist(data(:,d), linspace(0, 20, 100));
    else
        hist(data(:,d), linspace(-20, 20, 200));
    end
end
%%
cluster_model = models(1).clusters;
markers = '+x*ds^hv';

for d = 1:4
    
    [bin_counts bin_centres] = hist3(data(:,[2*d 2*d-1]), {linspace(-10, 10, 100), linspace(-10, 10, 100)});
    figure; imagesc(bin_counts); axis equal;
    set(gca, 'xtick',10:10:100, 'xticklabel', bin_centres{1}(10:10:end));
    set(gca, 'ytick', 10:10:100, 'yticklabel', bin_centres{2}(10:10:end));
    
    xlabel(['Dimension ', num2str(2*d-1)]);
    ylabel(['Dimension ', num2str(2*d)]);
    
    title(['Dims ', num2str(2*d - 1), ' and ', num2str(2*d)]);
    hold on;
    
    for k = 1:8
        dim2 = interp1(bin_centres{1}, 1:100, cluster_model.Means(k,2*d));
        dim1 = interp1(bin_centres{2}, 1:100, cluster_model.Means(k,2*d-1));
        plot(dim1, dim2, ['m' markers(k)], 'MarkerSize', 15);
        
        [transform, radii] = eig(cluster_model.CovMats{k}(2*d-1:2*d, 2*d-1:2*d));
        [x y] = ellipse(2*sqrt(radii(1,1)), 2*sqrt(radii(2,2)), cluster_model.Means(k,2*d-1), cluster_model.Means(k,2*d), transform(:,1));
        x = interp1(bin_centres{2}, 1:100, x);
        y = interp1(bin_centres{1}, 1:100, y);
        plot(x,y,'m');
    end
    
end
%%
%For data in phase magnitude form, lets look at some histogram images of
%data density across sets of dimensions

% 1) Magnitude vs Phase in each level
for lev = 1:4
    
    [bin_counts bin_centres] = hist3(data(:,[2*lev-1 2*lev]), {linspace(0, 10, 100), linspace(-pi/2, pi/2, 100)});
    figure; imagesc(bin_counts); axis image;
    set(gca, 'xtick',10:10:100, 'xticklabel', bin_centres{2}(10:10:end));
    set(gca, 'ytick', 10:10:100, 'yticklabel', bin_centres{1}(10:10:end));
    
    xlabel(['Level ', num2str(lev), ': Phase']);
    ylabel(['Dimension ', num2str(lev), ': Magnitude']);
    
    hold on;
    
end
%%
%2) Magnitude in L vs magnitude in L+1
for lev = 1:3
    
    [bin_counts bin_centres] = hist3(data(:,[2*lev-1 2*lev+1]), {linspace(0, 10, 100), linspace(0, 10, 100)});
    figure; imagesc(bin_counts); axis image;
    set(gca, 'xtick',10:10:100, 'xticklabel', bin_centres{2}(10:10:end));
    set(gca, 'ytick', 10:10:100, 'yticklabel', bin_centres{1}(10:10:end));

    xlabel(['Level ', num2str(lev+1), ': Magnitude']);
    ylabel(['Dimension ', num2str(lev), ': Magnitude']);
    
    hold on;
    
end
%%
c_means = models(1).clusters.Means;
[num_k num_dims] = size(c_means);

covar_mats = zeros(num_dims,num_dims,6);
for k = 1:num_k
    covar_mats(:,:,k) = models(1).clusters.CovMats{k};
end
c_probs = models(1).clusters.ClusterProbs;

for dim1 = num_dims:-1:1
    for dim2 = 1:dim1 - 1
        
        lev1 = ceil(dim1/2);
        lev2 = ceil(dim2/2);
        if rem(dim1, 2)
            %dim1 is mag
            lims1 = linspace(0, 3, 100);
            ylab = ['Level ', num2str(lev1), ': Magnitude'];
        else
            %dim1 is phase
            lims1 = linspace(-pi/2, pi/2, 100);
            ylab = ['Level ', num2str(lev1), ': Phase'];
        end
        if rem(dim2, 2)
            %dim2 is mag
            lims2 = linspace(0, 3, 100);
            xlab = ['Level ', num2str(lev2), ': Magnitude'];
        else
            %dim2 is phase
            lims2 = linspace(-pi/2, pi/2, 100);
            xlab = ['Level ', num2str(lev2), ': Phase'];
        end
            
        [bin_counts bin_centres] = hist3(data(:,[dim1 dim2]), {lims1, lims2});
        dims_gmm = gmdistribution(c_means(:,[dim2 dim1]), covar_mats([dim2 dim1],[dim2 dim1],:),c_probs);
        [x y] = meshgrid(bin_centres{2}, bin_centres{1});        
        dims_pdf = reshape(pdf(dims_gmm, [x(:) y(:)]), size(x));
        
        figure; 
        
        subplot(1,2,1); imagesc(bin_counts(1:end-1, 1:end-1)); axis image;
        set(gca, 'xtick',10:10:100, 'xticklabel', round(100*bin_centres{2}(10:10:end))/100);
        set(gca, 'ytick', 10:10:100, 'yticklabel', round(100*bin_centres{1}(10:10:end))/100);
        xlabel(xlab);
        ylabel(ylab);
    
        subplot(1,2,2); imagesc(dims_pdf); axis image;
        set(gca, 'xtick',10:10:100, 'xticklabel', round(100*bin_centres{2}(10:10:end))/100);
        set(gca, 'ytick', 10:10:100, 'yticklabel', round(100*bin_centres{1}(10:10:end))/100);
        xlabel(xlab);
        ylabel(ylab);
        f_name = ['C:\isbe\dev\background\location\figures\counts_vs_model_', num2str(dim1), '_', num2str(dim2), '.eps'];
        saveas(gcf, f_name, 'psc2');
%         for k = 1:8
%             c1 = interp1(bin_centres{1}, 1:100, cluster_model.Means(k,dim1));
%             c2 = interp1(bin_centres{2}, 1:100, cluster_model.Means(k,dim2));
%             plot(c2, c1, 'mx', 'MarkerSize', 15);
% 
%             [transform, radii] = eig(cluster_model.CovMats{k}([dim2 dim1], [dim2 dim1]));
%             [x y] = ellipse(2*sqrt(radii(1,1)), 2*sqrt(radii(2,2)), cluster_model.Means(k,dim2), cluster_model.Means(k,dim1), transform(:,1));
%             x = interp1(bin_centres{2}, 1:100, x);
%             y = interp1(bin_centres{1}, 1:100, y);
%             plot(x,y,'m');
%         end
    
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Now let's look at the posterior probs for each cluster for points in
% normal mammogramphic patches

%Load a model and conevrt into Matlab's gmm object
models(1).clusters = u_load('C:\isbe\dev\background\models\dual_tree\ilp_mag\normal512_k10_size_512_uid_1_std.mat');
c_means = models(1).clusters.Means;
[num_k num_dims] = size(c_means);
covar_mats = zeros(num_dims,num_dims,6);
for k = 1:num_k
    covar_mats(:,:,k) = models(1).clusters.CovMats{k};
end
c_probs = models(1).clusters.ClusterProbs;
gmm1 = gmdistribution(c_means, covar_mats, c_probs);
%%
%Randomly select a normal patch
normal_list = dir('C:\isbe\dev\background\dual_tree\normal512\*.mat');
im_list = dir('C:\isbe\dev\background\images\normal512\*.bmp');
ridx = randsample(1:length(normal_list), 1);
normal_dt = u_load(['C:\isbe\dev\background\dual_tree\normal512\', normal_list(ridx).name]);
normal_im = imread(['C:\isbe\dev\background\images\normal512\', im_list(ridx).name]);
normal_ilp = mb_dual_tree_transform(normal_dt);
%%
%Convert the ILP construct into a list of feature vectors and divide the
%magnitude dimensions by the standard devs from the model
[r c] = size(normal_ilp{1}(:,:,1));
[cols rows] = meshgrid(1:c, 1:r);
normal_data = mb_get_dt_vector([rows(:) cols(:)], normal_dt, 4);
normal_data(:,1:2:num_dims) = normal_data(:,1:2:num_dims) ./ repmat(stds(1:2:num_dims), r*c, 1);

% Compute posterior probabilities for each cluster in the model
post_probs = posterior(gmm2, normal_data(:,1:6));

%%
%Look at the posterior probs for each cluster
% 10 clusters for k-means form the following groups:
% {6}
% {9}
% {2, 4, 5, 8}
% {1, 3, 7, 10}
clust_im1 = reshape(post_probs(:,6), r, c);
clust_im2 = reshape(post_probs(:,9), r, c);
clust_im3 = reshape(sum(post_probs(:,[2 4 5 8]), 2), r, c);
clust_im4 = reshape(sum(post_probs(:,[1 3 7 10]), 2), r, c);

figure; imagesc(normal_im); axis image; colormap(gray(256));
figure; imagesc(clust_im1); axis image;
figure; imagesc(clust_im2); axis image;
figure; imagesc(clust_im3); axis image;
figure; imagesc(clust_im4); axis image;

%%
figure;
subplot(2,2,1);
imagesc(normal_im); axis image; colormap(gray(256));
title('Patch of real normal mammogram');
subplot(2,2,2);
imagesc(log10(clust_im1+clust_im2)); axis image; colorbar;
title('Posterior probabilities for cluster group most closely related to structure - log10 scaled');
subplot(2,2,3);
imagesc(clust_im3); axis image; colorbar;
title('Posterior probabilities for 1 of 2 cluster groups associated with general texture');
subplot(2,2,4);
imagesc(clust_im4); axis image; colorbar;
title('Posterior probabilities for 1 of 2 cluster groups associated with general texture');
f_name = 'C:\isbe\dev\background\location\figures\posterior_probs.eps';
saveas(gcf, f_name, 'psc2');
    

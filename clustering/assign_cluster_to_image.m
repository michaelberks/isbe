function [cluster_image] = assign_cluster_to_image(cluster_model, image_in)
% Given a cluster model of sample window data, assign each pixel in an
% image to it's nearest cluster

%get the number of clusters
k = size(cluster_model.Means, 1);

%get the window size
window_size = sqrt(size(cluster_model.Means, 2));

%Pre-compute the inverses of the cluster covariance matrices to save time
%later
covari = cell(1,40);
for cl = 1:k
    covari{cl} = pinv(cluster_model.CovMats{cl});
end

[rows cols] = size(image_in);
overlap = (window_size - 1) / 2;

%Pre-allocate cluster image to store assigned cluster
cluster_image = zeros(size(image_in));

%for each pixel in the image, ignoring edge pixels
for row = overlap+1:rows-overlap
    for col = overlap+1:cols-overlap;
        
        %extract sample window
        samp_window = sample_window(image_in, window_size, row, col);
        %pre-allocate vector to store Mahalonobis distances
        mahals = zeros(k,1);
        for cl = 1:k
            %Compute mahalanobis distance to each cluster
            x_minus_mu = cluster_model.Means(cl,:) - samp_window(:)';            
            mahals(cl) = x_minus_mu * covari{cl} * x_minus_mu';
        end
        % Save the nearest cluster to the cluster image;
        [dummy min_idx] = min(mahals);
        cluster_image(row, col) = min_idx;
    end
end

%Display both images
% figure; imagesc(image_in); colormap(jet(256)); axis image;
% figure; imagesc(cluster_image); colormap(jet(256)); axis image;


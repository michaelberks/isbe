%_------------------------------------------------------------------------
% Script showing how to cluster the markup of several observers for a set
% of nailfold images
%-------------------------------------------------------------------------
%%
%-------------------------------------------------------------------------
%% 1) For better ground truth, make up apex clusters using all available
% markers

rsa_dir = 'rsa_study/';
nailfoldroot = 'C:\isbe\nailfold\';

%Set paths to directories containing the images, image markup, and where
%you want to save the clustered markup
image_dir = [nailfoldroot 'data/' rsa_dir 'test/images/'];
markup_dir = [nailfoldroot 'data/' rsa_dir 'markup/'];
cluster_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_clusters/'];
post_cluster_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_clusters_merged/'];

mkdir(cluster_dir);
mkdir(post_cluster_dir);

%Create a cell array of marker names (that match the folder names in the
%markup dir - it may be easier just to explicitly list these)
markers = dir(markup_dir);
markers = markers(3:end);
discard_markers = struct2cell(markers);
markers(~cell2mat(discard_markers(4,:))) = []; 
    
%Get list of the images to process
 %I've converted the image to .mat files on my computer, you'll probably want to change the extension to .png
test_list = dir([image_dir '*.mat']);

%% Loop through the set of images computing the megered clusters
%Set default parameters
min_apex_dist = 20; %The minimum distance allowed between separate clusters (the larger out of this value or the apex width is used for each cluster)
include_nondistal = 1; %Record the position of non-distal vessels, and merge these with vessel marked as distal by a different observer if necessary

markers_per_image = zeros(length(test_list),1);
for i_im = 1:length(test_list)
    
    im_num = test_list(i_im).name(1:6);      
    [vessels] = cluster_vessel_apices(im_num, markup_dir, markers, min_apex_dist, include_nondistal);
    save([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    markers_per_image(i_im) = length(vessels.markers);
end
%--------------------------------------------------------------------------
%% Now try and post merge cluster based on the vessel detection images

vessel_prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/257273/'];
dist_thresh = 100; %What's the largest separation allowed between cluster centres
patch_sz2 = 50; %What's the halfwidth of the patch in which we look for connectivity
connect_thresh = 0.5; %What's the threshold for deciding 2 clusters are connected?
n_connect_pts = 20; %How many points (spaced between 0 and 1) to test the connectivity of pairs
for i_im = 1:length(test_list);
    
    if markers_per_image(i_im) > 1
        im_num = test_list(i_im).name(1:6);  
        vessel_prob = u_load([vessel_prob_dir im_num '_pred.mat']);
        load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
        
        [vessels] = post_merge_apexes(vessels, vessel_prob, dist_thresh, patch_sz2,...
            connect_thresh, n_connect_pts);
        
        save([post_cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    end
end
%-------------------------------------------------------------------------
%% Visually inspect the clusters
for i_im = 1:20
    
    im_num = test_list(i_im).name(1:6);    
    if exist([post_cluster_dir im_num '_apex_clusters.mat'], 'file')
        load([post_cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    else
        load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    end
    nailfold = u_load(['C:\isbe\nailfold\data\rsa_study\test\images\' im_num '.mat']);
    
    figure; imgray(nailfold);
    plot_apex_clusters(vessels, gcf);
end


model = model_2_1;

for ii = 1:10
    figure;
    imagesc(reshape(model_2_1.Means(ii,:), 11, 11));
    colormap(jet(256)); axis image;
end

%%

for level = 2:5
    for ori = 1:2
        figure;
        imagesc(cluster_image{level, ori});
        colormap(jet(256)); axis image;
    end   
end

%%
syn_args.SamplePyramid = u_load('C:\isbe\dev\background\pyramid\lines\line001_pyramid.mat');
syn_args.ClusterImage = u_load('C:\isbe\dev\background\misc\lines\line001_pyramidcluster_image.mat');
syn_args.FilledImage = logical(padarray(zeros(128), [64 64], 1));
syn_args.ModelDir = 'C:\isbe\dev\background\results\lines\';
syn_args.ModelName = 'lines_model';
syn_args.CutOffLevel = 5;
syn_args.Plot = 1;
syn_args.PyrLevels = 5;
syn_args.PyrOrientations = 2;
syn_args.ForceCluster = 1;
      
syn_args.SaveFile = [mberksroot, 'background\syn\test_lines_model'];
[synthesised_image, new_pyr] = mb_gmm_pyr_synthesis(syn_args);
figure; imagesc(synthesised_image); colormap(gray(256)); axis image;

%%

for level = 2:5
    for ori = 1:2
        figure;
        imagesc(new_pyr{level, ori});
        colormap(jet(256)); axis image;
    end   
end

%%
for level = 2:5
    for ori = 1:2
        figure;
        imagesc(syn_args.SamplePyramid{level, ori});
        colormap(jet(256)); axis image;
    end   
end
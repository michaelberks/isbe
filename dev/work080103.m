for ii = 1:7
    for jj = 1:1
        figure; imagesc(pyr10{ii,jj}); axis image; colormap(jet(256));
        if ii == 1 || ii == 7
            break;
        end
    end
end

%%
pyramid = u_load('C:/isbe/dev/background/pyramid/normal_2/normal001.bmp_pyramid.mat');
[levels orientations] = meshgrid(2:3, 1:5);

%%
filled_image = ones(size(im));
filled_image(50:249, 50:249) = 0;

syn_args.PathToTextureGMM = 'C:/isbe/dev/background/results/pyramid/normal_2_model_2_4_a.mat';
syn_args.SeededImage =  im;
syn_args.FilledImage = filled_image;
syn_args.Uint8Mode = 0;

syn_args.SynthesisMode = 'pixel-wise';
[syn_image clust_image] = mb_gmm_tex_synthesis(syn_args);

%%
c_model = u_load('C:/isbe/dev/background/results/pyramid/normal_2_model_3_4.mat');
clims = [-1 1];
for ii = 21:40
    mean_window = reshape(c_model.Means(ii,:), 15, 15);
    figure; imagesc(mean_window, clims); axis image;
end
%%
im = pyramid{levels(4), orientations(4)};
cluster_image = zeros(size(im));
clear pyramid;
%%
covari = cell(1,40);
for cl = 1:40
    covari{cl} = pinv(covars{cl});
end
[rows cols] = size(im);
overlap = 7;
for row = overlap+1:rows-overlap
    for col = overlap+1:cols-overlap;
        
        samp_window = sample_window(im, 15, row, col);
        mahals = zeros(40,1);
        for cl = 1:40
            x_minus_mu = means(cl,:) - samp_window(:)';
            
            pinv_C = pinv(covars{cl});
            mahals(cl) = x_minus_mu * pinv_C * x_minus_mu';
        end
        cluster_image(row, col) = max(mahals);
    end
end
figure; imagesc(syn_image); colormap(jet(256)); axis image;
figure; imagesc(cluster_image); colormap(jet(256)); axis image;
            
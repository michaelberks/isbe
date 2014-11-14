i1 = double(imread('C:/isbe/dev/background/images/normal_2/normal001.bmp'));
%%
filled_image = ones(size(i1));
filled_image(100:199, 100:199) = 0;

syn_args.PathToPyramidGMM = [mberksroot, '/background/results/pyramid/normal_2_model'];
syn_args.SampleImage = i1;
syn_args.FilledImage = filled_image;
          
%           'SynthesisMode', 'simple', ...
% 		  'MovieFilename', [], ...
% 		  'NormaliseClusterProbs', (1==1), ...
% 		  'SaveFrequency', 1000,...
%           'SaveFile', 'C:\isbe\dev\background\syn\temp\syn_image',...
%           'PyrLevels', 5,...
%           'PyrOrientations', 5,...
%           'CutOffLevel', 2

profile on;
[synthesised_image] = mb_gmm_pyr_synthesis(syn_args);
profile viewer;

%%
[p_rows p_cols] = find(~filled_image);
%p_idx = sub2ind(size(args.SampleImage), p_rows, p_cols);

area_idx = cell(7, 1);

for level = 6:-1:2
    new_rows = ceil(p_rows/2^(level-2));
    new_cols = ceil(p_cols/2^(level-2));
    area_idx{level} = sub2ind(size(pyramid1{level,1}), new_rows, new_cols);
end

%%
model_means = zeros(27,40);
for ii = 1:40
    model_means(:,ii) = diag(model.CovMats{ii});
end
%%
[synthesised_image, pyramid2] = mb_gmm_pyr_synthesis(syn_args);
pyramid2 = test(pyramid2, 5, 5, p_sizes);
figure; imagesc(synthesised_image); axis image; colormap(gray(256));
for ii = 2:3
    for jj = 1:5
        figure; imagesc(pyramid3{ii,jj}); colormap(gray(256));
    end
end
%%
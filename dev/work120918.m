warning('off', 'ASYM:unexpectedArgument');

%Arguments to extract features
decomposition_args.decomp_type = 'dt'; %Use the dual-tree
decomposition_args.win_size = 3; %Window size about pixel to sample features
decomposition_args.levels = 1:5; %Levels to include
decomposition_args.feature_shape = 'rect'; %Keep this as rect, although there are other things you can play with (see sample_dt_data)
decomposition_args.feature_type = 'conj'; %Again, 'conj' works best, but can use 'mag' for magnitude only, or 'phase' (see convert_complex_representation)
decomposition_args.do_max = 0; %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0; %Try and make rotation invariant (best left switched off)
decomposition_args.use_nag = 0; %Use the NAG toolbox to interpolate (best left swicthed off)
%%
im = imread('C:\isbe\asymmetry_project\data\fibre\CCMdatabase_Analysed\T090173-0000-02-0038.bmp');
r_map = (im(:,:,1) > 0.9) & (im(:,:,2) < 0.1) & (im(:,:,3) < 0.1);
g_map = (im(:,:,1) < 0.1) & (im(:,:,2) > 0.9) & (im(:,:,3) < 0.1);
b_map = (im(:,:,1) < 0.1) & (im(:,:,2) < 0.1) & (im(:,:,3) > 0.9);
junction_map = bwmorph(g_map, 'thin', inf);
fibre_map = r_map | b_map | junction_map;
red_fibre_map = r_map | junction_map;
blue_fibre_map = b_map | junction_map;

figure; imgray(im);
figure; imgray(fibre_map);
figure; imgray(red_fibre_map);
figure; imgray(blue_fibre_map);
figure; imgray(junction_map);
%%

rows = randi(size(im,1), 10, 1);
cols = randi(size(im,2), 10, 1);

[responses] = compute_filter_responses(im, decomposition_args);
[features] = sample_image_features(responses, rows, cols, decomposition_args);
%%
fibre_list = dir([asymmetryroot 'data\fibre\CCMdatabase_Analysed\*.bmp']);

for ii = 1:20

    %Load in the original fibre image and the annotated copy
    orig_im = imread([asymmetryroot 'data\fibre\CCMdatabase\' fibre_list(ii).name]);
    if size(orig_im,3) == 3
        orig_im = rgb2gray(orig_im);
    end
    gabor = compute_gabor_responses(orig_im, 2, 8);
    max_gabor = max(abs(gabor), [], 4);
    
    gaussian = imfilter(double(orig_im), fspecial('gaussian', 12, 2));
    
    figure;
    subplot(2,2,1); imgray(orig_im);
    subplot(2,2,2); imgray(max_gabor);
    subplot(2,2,3); imgray(gaussian);
    subplot(2,2,4); imgray(max_gabor - gaussian);
end
%%
y = 0:0.05:5;
k = 2;

for x = 1:4
    figure; hold all; title(['x = ' num2str(x)]);
    leg_text = cell(0,1);
    for a = 1:4
        g = (x + a*y) ./ (1 + exp(-2*k*(x - a*y)));
         plot(y, g);
         leg_text{a} = ['\alpha = '  num2str(a)];
    end
    legend(leg_text);
end
%%
[x y] = meshgrid(-5:0.05:5, 0:0.05:5);
k = 1;
a = 1;
g = (x + a*y) ./ (1 + exp(-2*k*(x - a*y)));

figure; surf(x, y, g);
xlabel('x');
ylabel('y');

gx = x ./ (1 + exp(-2*k*x));
figure; plot(x(1,:), gx);
%%
fibre_list = dir([asymmetryroot 'data\fibre\CCMdatabase_Analysed\*.bmp']);
a = 0.5;
k = 1;
for ii = 19%1:20

    %Load in the original fibre image and the annotated copy
    orig_im = imread([asymmetryroot 'data\fibre\CCMdatabase\' fibre_list(ii).name]);
    if size(orig_im,3) == 3
        orig_im = rgb2gray(orig_im);
    end
    gabor = real(compute_gabor_responses(orig_im, 2, 8));
    [~, max_idx] = max(abs(gabor), [], 4);
    max_gabor = zeros(size(orig_im));
    for ori = 1:8
        mask = max_idx == ori;
        band = gabor(:,:,1,ori);
        max_gabor(mask) = band(mask);
    end
        
    gaussian = compute_gaussian_responses(orig_im, 2);
    
    r_pos = max_gabor + a*gaussian;
    r_neg = max_gabor - a*gaussian;
    
    dual_model = r_pos ./ (1 + exp(-2*k*r_neg));
    hack = max_gabor ./ (1 + exp(-2*k*max_gabor));
    
    figure;
    ax(1) = subplot(2,3,1); imgray(gaussian); title('original image');
    ax(2) = subplot(2,3,2); imgray(max_gabor); title('Gabor response - real part maximised across orientation');
    ax(3) = subplot(2,3,3); imgray(r_pos); title('\Gamma_{p}');
    ax(4) = subplot(2,3,4); imgray(r_neg); title('\Gamma_{n}');
    ax(5) = subplot(2,3,5); imgray(dual_model); title('Dual model');
    ax(6) = subplot(2,3,6); imgray(hack); title('Gabor through sigmoid');
    linkaxes(ax);
end

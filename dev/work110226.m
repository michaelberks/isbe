nd = 16;
num_levels = 4;
win_size = 1;

D = 6*num_levels*win_size^2;
nd2 = 2^nd;
a_dims = ceil(rand(nd,1)*D);
b_dims = ceil(rand(nd,1)*(D-1));
b_dims(b_dims >= a_dims) = b_dims(b_dims >= a_dims)+1;

brief_sums = sparse(nd2,1);
brief_counts = sparse(nd2,1);
%
for ii = 200:200
    for jj = 1:14
        display(['Processing block (' num2str(ii) ', ' num2str(jj) ')']);
        tic;
        X = u_load(['G:\asymmetry_project\data\synthetic_data\real512_dt\' zerostr(ii,3) '\X_' zerostr(jj,2) '.mat']);
        y = u_load(['G:\asymmetry_project\data\synthetic_data\real512_dt\' zerostr(ii,3) '\y_' zerostr(jj,2) '.mat']);
        toc;
        X(:,:,:,num_levels+1:end) = [];

        if win_size == 1
            X = X(:,5,:,:);
        end
    
        X = convert_dt_representation(X, 'feature_type', 'mag', 'win_size', win_size);
        tic;
        brief_codes = bin_row_to_dec(X(:,a_dims) < X(:,b_dims)+2);
        toc;
        tic;
        brief_sums = brief_sums + sparse(brief_codes, 1, y, nd2,1);
        brief_counts = brief_counts + sparse(brief_codes, 1, 1, nd2,1);
        toc;
        
%         for kk = 1:size(X,1);
%             brief_sums(brief_codes(kk)) = brief_sums(brief_codes(kk)) + y(kk);
%             brief_counts(brief_codes(kk)) = brief_counts(brief_codes(kk)) + 1;
%         end
    end
end
%%
px_per_mm = 100 / 9;
test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\';
prob_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\';

%
%Compute line and orientation maps for analytic g2d methods
param_dir = 'g2d_scales_16';

mkdir([prob_dir param_dir '\lines\']);
mkdir([prob_dir param_dir '\orientations\']);

for ii = 1:100
    display(['Processing image ', num2str(ii) ' of 100']);
    %Load in test image
    test_im = u_load([test_dir '\image' zerostr(ii,3) '.mat']);
    
    [line_orientation, line_map] = karssemeijer_line_detection(...
        test_im,...
        'line_scales', px_per_mm*[0.1 0.17 0.29],...
        'grad_scale', 10,...
        'grad_ori_thresh', pi/6,...
        'grad_strength_thresh', 25,...
        'line_strength_thresh', 0,...
        'binary_map', 0);
    
    if ii < 11
        figure; colormap(gray(256));
        subplot(1,2,1); imagesc(test_im); axis image;
        subplot(1,2,2); imagesc(line_map); axis image;
    end
        
    save([prob_dir param_dir '\lines\image' zerostr(ii,3) '_class.mat'], 'line_map');
    save([prob_dir param_dir '\orientations\image' zerostr(ii,3) '_class.mat'], 'line_orientation');

end

    
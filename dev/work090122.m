%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will build a group a dual-tree for simple black and white
% images of lines, curved lines and blobs (filled ellipses)

%%%
%%
% mkdir C:\isbe\dev\background\images\lecb
% mkdir C:\isbe\dev\background\dual_tree\lecb
%
% 1) Build lines
for ii = 21:40
    line_256 = ones(256);
    half_length = randsample(32:64, 1);
    width = randsample(1:4, 1);
    
    line_256(128:127+width, 129-half_length:128+half_length) = 2;
    
    rnd_angle = 360*rand(1);
    temp = imrotate(line_256, rnd_angle, 'crop');
    line_128 = temp(65:192, 65:192);
    
%     figure; imagesc(line_128); colormap(gray); axis image;
    dt_line_128 = dtwavexfm2(line_128, 5, 'near_sym_b','qshift_b');
    
    im_name = ['C:\isbe\dev\background\images\lecb\line_', zerostr(ii,2)];
    dt_name = ['C:\isbe\dev\background\dual_tree\lecb\dual_tree_line_', zerostr(ii,2)];
    save(im_name, 'line_128');
    save(dt_name, 'dt_line_128');
    
end

%%
% 2) Build edges

edge_256 = [ones(128,256); 2*ones(128,256)];

for ii = 1:20;
    
    rnd_angle = 360*rand(1);
    temp = imrotate(edge_256, rnd_angle, 'crop');
    edge_128 = temp(65:192, 65:192);
    
%     figure; imagesc(edge_128); colormap(gray); axis image;
%     dt_edge_128 = dtwavexfm2(edge_128, 5, 'near_sym_b','qshift_b');
%     
%     im_name = ['C:\isbe\dev\background\images\lecb\edge_', zerostr(ii,2)];
%     dt_name = ['C:\isbe\dev\background\dual_tree\lecb\dual_tree_edge_', zerostr(ii,2)];
%     save(im_name, 'edge_128');
%     save(dt_name, 'dt_edge_128');
    
end

%%
% 3) Build curves (as open ellipses)
[xx yy] = meshgrid(1:256, 1:256);

for ii = 21:40
    
    v_axis = randsample(16:32, 1);
    h_axis = randsample(v_axis:48, 1);
    width = randsample(2:4, 1);
    
    ellipse1 = ((xx - 128) / h_axis).^2 + ((yy - 128) / v_axis).^2 < 1;
    ellipse2 = ((xx - 128) / (h_axis-width)).^2 + ((yy - 128) / (v_axis-width)).^2 >= 1;
    
    ellipse_256 = (ellipse1 & ellipse2) + 1;
    if rand(1) > 0.5
        ellipse_256(1:128,:) = 1;
    else
        ellipse_256(129:end,:) = 1;
    end
    
    rnd_angle = 360*rand(1);
    temp = imrotate(ellipse_256, rnd_angle, 'crop');
    ellipse_128 = temp(65:192, 65:192);
    
%     figure; imagesc(ellipse_128); colormap(gray); axis image;
    dt_ellipse_128 = dtwavexfm2(ellipse_128, 5, 'near_sym_b','qshift_b');
    
    im_name = ['C:\isbe\dev\background\images\lecb\ellipse_', zerostr(ii,2)];
    dt_name = ['C:\isbe\dev\background\dual_tree\lecb\dual_tree_ellipse_', zerostr(ii,2)];
    save(im_name, 'ellipse_128');
    save(dt_name, 'dt_ellipse_128');
%     
end

%%
% Build blobs
[xx yy] = meshgrid(1:256, 1:256);

for ii = 1:20
    
    v_axis = randsample(8:16, 1);
    h_axis = randsample(8:16, 1);
    
    blob_256 = (((xx - 128) / h_axis).^2 + ((yy - 128) / v_axis).^2 < 1) + 1;
    
    rnd_angle = 360*rand(1);
    temp = imrotate(blob_256, rnd_angle, 'crop');
    blob_128 = temp(65:192, 65:192);
    
    %figure; imagesc(blob_128); colormap(gray); axis image;
    dt_blob_128 = dtwavexfm2(blob_128, 5, 'near_sym_b','qshift_b');
    
    im_name = ['C:\isbe\dev\background\images\lecb\blob_', zerostr(ii,2)];
    dt_name = ['C:\isbe\dev\background\dual_tree\lecb\dual_tree_blob_', zerostr(ii,2)];
    save(im_name, 'blob_128');
    save(dt_name, 'dt_blob_128');
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we can extract feature vectors and histogram/cluster them as with
% real data
data = mb_hist_dual_tree('C:\isbe\dev\background\dual_tree\lecb\', 128, 4);
%%
num_dims = size(data, 2);
for dim1 = 3:num_dims
    for dim2 = dim1 + 1:num_dims
        
        lev1 = ceil(dim1/2);
        lev2 = ceil(dim2/2);
        if rem(dim1, 2)
            %dim1 is mag
            nbins1 = 101;
            lims1 = 2:101;
            ylab = ['Level ', num2str(lev1), ': Magnitude'];
        else
            %dim1 is phase
            nbins1 = 100;
            lims1 = 1:100;
            ylab = ['Level ', num2str(lev1), ': Phase'];
        end
        if rem(dim2, 2)
            %dim2 is mag
            nbins2 = 101;
            lims2 = 2:101;
            xlab = ['Level ', num2str(lev2), ': Magnitude'];
        else
            %dim2 is phase
            nbins2 = 100;
            lims2 = 1:100;
            xlab = ['Level ', num2str(lev2), ': Phase'];
        end
            
        [bin_counts bin_centres] = hist3(data(:,[dim1 dim2]), [nbins1 nbins2]); % {lims1, lims2}
        [bin_counts2 bin_centres2] = hist3(data2(:,[dim1 dim2]), [nbins1 nbins2]); % {lims1, lims2}
        
        figure; 
        subplot(1,2,1); imagesc(log10(bin_counts(lims1, lims2))); axis image;
        set(gca, 'xtick',10:10:100, 'xticklabel', bin_centres{2}(10:10:end));
        set(gca, 'ytick', 10:10:100, 'yticklabel', bin_centres{1}(10:10:end));
        xlabel(xlab);
        ylabel(ylab);
        
        subplot(1,2,2); imagesc(log10(bin_counts2(lims1, lims2))); axis image;
        set(gca, 'xtick',10:10:100, 'xticklabel', bin_centres2{2}(10:10:end));
        set(gca, 'ytick', 10:10:100, 'yticklabel', bin_centres2{1}(10:10:end));
        xlabel(xlab);
        ylabel(ylab);
    
        

    
    end
end

%%
% Why do we get points with ILP phase close to -pi/2 - i.e. suggesting we
% have a negative line, when all lines should be positive?
colors = 'rgbm';
for ii = 1:80
    im = u_load(['C:\isbe\dev\background\images\lecb\', im_list(ii).name]);
    dt = u_load(['C:\isbe\dev\background\dual_tree\lecb\', dt_list(ii).name]);
    ilp = mb_dual_tree_transform(dt);
    
    if ~rem(ii-1,4)
        figure;
    end
    
    subplot(2,2, rem(ii-1, 4)+1); imagesc(im); axis image; colormap(gray); hold on;
    
    for level = 1:4
        
        [y x] = find(angle(max(ilp{level}, [], 3)) > 1.0 & angle(max(ilp{level}, [], 3)) < 1.3); 
        dim = 2^level;
        st = (dim+1)/2;
        
        plot(dim*(x-1) + st, dim*(y-1) + st, [colors(level), 'x']);
    end
end

%%
[px, lx, bx, mx] = st_pca(data, 1); clear data; pack
[px2, lx2, bx2, mx2] = st_pca(data2, 1); clear data; pack

for ii = 1:8; [prctile(bx(:,ii), [2 98]) prctile(bx2(:,ii), [2 98])], end
%%
lims1 = linspace(-2, 6, 102);

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
        
%         f_name = ['C:\isbe\dev\background\location\figures\pca_counts_', num2str(dim1), '_', num2str(dim2), '.eps'];
%         saveas(gcf, f_name, 'psc2');
        
    
    end
end

        
    
    
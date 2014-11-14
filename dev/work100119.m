norm_list = dir('C:\isbe\dev\background\images\normal512\*.bmp');
mkdir C:\isbe\dev\background\images\normal_smooth512\

for ii = 20:40;%1:length(norm_list)

    norm_bg = double(imread(['C:\isbe\dev\background\images\normal512\', norm_list(ii).name]));

    norm_dt = dtwavexfm2(norm_bg, 7);
    %%
    norm_dt2 = norm_dt;
    norm_dt3 = norm_dt;

    for lev = 1:6
        for ori = 1:6
            subband = norm_dt2{lev}(:,:,ori);
            [dummy sort_idx] = sort(abs(subband(:)));
            sort_idx(end-round(2*end/5):end) = [];
            subband(sort_idx) = 0;
            norm_dt2{lev}(:,:,ori) = subband;
        end
        norm_dt3{lev}(:) = 0;
    end
    norm_bg2 = dtwaveifm2(norm_dt2);

    norm_noise = norm_bg - norm_bg2;
    norm_coarse = dtwaveifm2(norm_dt3);

    %bg = norm_coarse(193:end-192,193:end-192)+norm_noise(193:end-192,193:end-192);
    bg = norm_coarse + norm_noise;
    
%     if rand > .5
%         bg(17:end-16, 60:68) = bg(17:end-16, 60:68) + 16;
%     else
%         bg(17:end-16, 60:68) = bg(17:end-16, 60:68) + 8;
%     end
%    figure; imagesc(bg); axis image; colormap(gray(256));
%     figure; 
%     subplot(2,2,1); imagesc(norm_bg2); axis image; colormap(gray(256));
%     subplot(2,2,2); imagesc(norm_noise); axis image; colormap(gray(256));
%     subplot(2,2,3); imagesc(norm_coarse); axis image; colormap(gray(256));
%     subplot(2,2,4); imagesc(bg); axis image; colormap(gray(256));
%     
   save(['C:\isbe\dev\background\images\normal_smooth512\bg', zerostr(ii,3)], 'bg');
    
end
%%
for ii = 1:30
    load(['C:\isbe\dev\background\images\normal_smooth128\bg', zerostr(ii,3)]);
    figure; imagesc(bg); axis image; colormap(gray(256));
end
    
%%
norm_mod_dt = dtwavexfm2(norm_mod, 5);
%%
for lev = 1:5
    figure;
    for ori = 1:6
        subplot(2,3,ori);
        imagesc(abs(norm_mod_dt{lev}(:,:,ori))); axis image; colorbar;
    end
    figure;
    for ori = 1:6
        subplot(2,3,ori);
        imagesc(abs(norm_dt{lev}(:,:,ori))); axis image; colorbar;
    end
end
%%
x = repmat(1:256, 256, 1);
y = repmat((1:256)', 1, 256);
 
orientations = (1:30:180)*pi/180;
for o = 1:length(orientations)
    a = sin(orientations(o));
    b = cos(orientations(o));
    c = -128*(a + b);
    dx = abs(a*x + b*y + c);
    
    for halfwidth = 2:4:32
        
        sigma2 = (halfwidth^2) / log(2);
        ymax = 1/sqrt(2*pi*sigma2);

        scaling = 1 / ymax;
        pos_line = scaling*exp(-(dx.^2 / sigma2)) / sqrt(2*pi*sigma2);
        figure; imagesc(pos_line); axis image; colormap gray; colorbar
    end
end
%%
norm_list = dir('C:\isbe\dev\background\images\normal512\*.bmp');
mkdir C:\isbe\dev\background\images\normal_128\
idx = 1;
for ii = 1:length(norm_list);

    norm_bg = double(imread(['C:\isbe\dev\background\images\normal512\', norm_list(ii).name]));
    
    for rr = [1 129 257 385]
        for cc = [1 385]
            bg = norm_bg(rr:rr+127, cc:cc+127);
            save(['C:\isbe\dev\background\images\normal_128\bg', zerostr(idx,3)], 'bg');
            idx = idx+1;
        end
    end   
end
%%
norm_list = dir('C:\isbe\dev\background\images\normal512\*.bmp');
mkdir C:\isbe\dev\background\images\normal_512a\
for ii = 1:length(norm_list);

    bg = double(imread(['C:\isbe\dev\background\images\normal512\', norm_list(ii).name]));
    save(['C:\isbe\dev\background\images\normal_512a\bg', zerostr(ii,3)], 'bg');

end
%%
%mb_make_new_function('generate_bar', 0, {'type', 'intensity', 
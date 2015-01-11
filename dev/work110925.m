load(['C:\isbe\nailfold\images\ncm\Visit1\Lhand\Digit4\X300\mb_reg\frames' zerostr(1, 2) '.mat'], 'frames_reg');
for ii = 1:16    
    if ii == 1
        imwrite(uint8(255*(frames_reg(20:end-9,10:end-9,ii)-100)/155), gray(256),...
            'P:\isbe\nailfold\figures\frames.gif',...
            'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.2);
    else
        imwrite(uint8(255*(frames_reg(20:end-9,10:end-9,ii)-100)/155), gray(256),...
            'P:\isbe\nailfold\figures\frames.gif',...
            'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end
%%
tile = imread('C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1FComp.bmp');
imwrite(uint8(255*(double(tile(20:end-9,10:end-9))-100)/155), 'C:\isbe\nailfold\presentations\figures\frame.bmp');

%%
tile = double(tile(20:end-9,10:end-9));

[nailfold_strength, nailfold_orientation] = ...
                gaussian_clover_line(tile, [2 4 8]);
nailfold_nms = mb_non_maximal_supp(nailfold_strength, nailfold_orientation);

nailfold_bao = bwareaopen(nailfold_nms > 0, 10);

[y_b x_b] = find(nailfold_bao);

figure; imagesc(tile); axis image; colormap(gray(256)); hold on;
plot(x_b, y_b, 'y.', 'markersize', 2);

%%
f = figure(...
    'windowstyle', 'normal',...
    'PaperPositionMode','auto');

imagesc(tile); axis image off; colormap(gray(256)); hold on;
print('-dtiff', '-noui', '-painters', f, '-r300', 'C:\isbe\nailfold\presentations\figures\frame_no_lines.tif');
plot(x_b, y_b, 'y.', 'markersize', 2);
print('-dtiff', '-noui', '-painters', f, '-r300', 'C:\isbe\nailfold\presentations\figures\frame_lines.tif');
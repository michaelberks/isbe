function read_mass_data2

bg = rgb2gray(imread('C:\isbe\project\shapes\background.bmp'));
for ii = 1:9;
    sh = imread(strcat('C:\isbe\project\shapes\shape00', num2str(ii), '.bmp'));
    [start_r start_c] = find(sh, 1);
    sb = bwtraceboundary(sh, [start_r, start_c], 'E');
    data(ii).shape_border = [sb(:,2), sb(:,1)];
    %px = bg(sh);
    pl = regionprops(bwlabel(sh, 4), 'PixelList');
    pxl = pl.PixelList;
    data(ii).shape_texture = pxl;%[pxl, double(px)];
    data(ii).shape_ROI = bg;
    
    %p = floor((length(shape_border(:,1)) / 150));
    %shape_vec = shape_border(1:p:end, :);
    %shape_mat(ii, :) = [shape_vec(1:150,1)', shape_vec(1:150,2)'];
end
for ii = 10:20;
    sh = imread(strcat('C:\isbe\project\shapes\shape0', num2str(ii), '.bmp'));
    [start_r start_c] = find(sh, 1);
    sb = bwtraceboundary(sh, [start_r, start_c], 'E');
    data(ii).shape_border = [sb(:,2), sb(:,1)];
    %px = bg(sh);
    pl = regionprops(bwlabel(sh, 4), 'PixelList');
    pxl = pl.PixelList;
    data(ii).shape_texture = pxl;%[pxl, double(px)];
    data(ii).shape_ROI = bg;
end
mass_input_data = data;
save C:\isbe\project\shapes\mass_data mass_input_data
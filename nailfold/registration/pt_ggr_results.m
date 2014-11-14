clc;
clear;
% close all;

datroot = 'U:\projects\nailfold\capture';
imgpath = fullfile(datroot, '\2012_10_18\Left\Digit4\x300\seq2');
ggrpath = fullfile(imgpath, 'ggr\Stage0001');

nframes = 30;

for iframe = 1:nframes
    d = dir(fullfile(imgpath, sprintf('*%04d.png', iframe)));
    imgfile = fullfile(imgpath, d(1).name);
    ptsfile = fullfile(ggrpath, [d(1).name,'.pts']);

    img = imread(imgfile);
    pts = read_pts(ptsfile);
    
    x = pts(1,1):2:pts(end,1);
    y = pts(1,2):2:pts(end,2);
    [xx,yy] = meshgrid(x,y);

    if (iframe == 1)
        intimg = interp2(double(img), xx, yy, '*linear');
    else
        intimg = intimg + interp2(double(img), xx, yy, '*linear');
    end
end

intimg = intimg / nframes;
figure(1); clf; colormap(gray(256));
imagesc(intimg);


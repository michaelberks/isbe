imgpath = 'U:\projects\nailfold\synthesis\_movies\zoom\ncm_video';

d = dir(fullfile(imgpath,'frame*.png'));
nImages = length(d);

outpath = fullfile(imgpath, 'looped');
if ~exist(outpath, 'dir')
    mkdir(outpath);
else
    delete(fullfile(outpath, 'frame*.png'));
end

% Blend some proportion of the frames
overlap = 0.1;
nOverlap = floor(overlap*nImages);

% Blend first few images
alpha = (1:nOverlap) / nOverlap;
for i = 1:nOverlap
    img1 = imread(fullfile(imgpath, d(i).name));
    img2 = imread(fullfile(imgpath, d(nImages-nOverlap+i).name));
    
    img =    alpha(i)  * img1 + ...
          (1-alpha(i)) * img2;
      
    imwrite(uint8(img), fullfile(outpath, d(i).name));   
end
nImages = nImages - nOverlap;

% Copy the rest
for i = nOverlap+1:nImages
    copyfile(fullfile(imgpath, d(i).name), ...
             fullfile(outpath, d(i).name));
end

clc; clear all;

imgpath = 'U:\data\satellite\ukcities\';

% do google roadmaps
imgdir = dir([imgpath,'googleroads/*.png']);
cols =	[	255,253,139; % minor (B) roads
			255,195,69; % major (A) roads
			149,209,123; % dual carriageway
			119,134,222; % motorway
			153,179,204 % water
		];
for i = 1:length(imgdir)
	imgname = imgdir(i).name;
	img1 = imread([imgpath,'googleroads/',imgname]);
	
	map = false(size(img1(:,:,1)));
	for c = 1:size(cols,1)
		map = map |	  (	img1(:,:,1)==cols(c,1) & ...
						img1(:,:,2)==cols(c,2) & ...
						img1(:,:,3)==cols(c,3) );
	end
		
	se = strel('disk',4);
	map = imclose(map,se);
	
	imwrite(uint8(map*255),[imgpath,'masks/googleroads/',imgname]);
end

return

% do overlays
imgdir = dir([imgpath,'overlays/*.png']);
for i = 1:length(imgdir)
	imgname = imgdir(i).name;
	img0 = imread([imgpath,imgname]);
	img1 = imread([imgpath,'overlays/',imgname]);
	diff = any(img0~=img1,3);
	imwrite(uint8(diff*255),[imgpath,'masks/overlays/',imgname]);
end

if 0
	figure(1); clf; 
		image(img0); axis('image');
	figure(3); clf; colormap(gray(256));
		image(uint8(diff*255)); axis('image');
end


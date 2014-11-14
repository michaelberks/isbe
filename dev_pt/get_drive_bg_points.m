% script to sample background points from retinal training images

% define image path and get list of PNG images
imgpath	= 'U:\data\retinal\drive\training\';
imgdir	= dir([imgpath,'images\*.png']);
mskdir	= dir([imgpath,'mask\*.png']);
lindir	= dir([imgpath,'1st_manual\*.png']);

% path where we'll store the coefficients
outpath = [asymmetryroot,'data\synthetic_data\drive\dt\training\'];

sampling_args = struct(...
    'feature_type', 'raw',...
	'levels',1:6,...
    'win_size', 3);

%% generate training data
fprintf('Generating training data...');
for i = 1:length(imgdir)
	% load image and take green channel
	img = imread([imgpath,'images/',imgdir(i).name]);
	img = img(:,:,2);
	
	% get manual markup
	imgmask = logical(imread([imgpath,'mask/',mskdir(i).name]));
	linemask = logical(imread([imgpath,'1st_manual/',lindir(i).name]));
	
	% define points >T from a line as background
	T = 0; % threshold distance
	bgmask = double(~linemask);
	bgmask(~linemask & imgmask) = inf;
	bgmask_dt = dtrans(bgmask);
	bgmask = (bgmask_dt>T);
	
	% pick points at random
	Nbg = sum(bgmask(:)); % number of background points
	inds = find(bgmask);
	inds = inds(randperm(length(inds))); % shuffle indices
	inds = inds(1:min(Nbg,10000)); % sample up to 10k points
	
	% compute DTCWT coefficients for whole image
	dt = compute_dual_tree(img,6,0);
	
	% convert indices to (row,col) pairs
	[rinds,cinds] = ind2sub(size(img),inds(:));

	% get these coefficients
	Xbg = sample_dt_data(dt,rinds,cinds,sampling_args);
	
	% save to file
	filename = sprintf('Xbg_%s.mat',imgdir(i).name(1:2));
	save([outpath,filename],'Xbg');
end
fprintf('done\n');


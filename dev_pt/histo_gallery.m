% clear;

outpath = 'U:\projects\mammography\data\histograms\';
outdir = dir(outpath);
outdir = outdir([outdir(:).isdir]);
outpath = fullfile(outpath,outdir(end).name);

load(fullfile(outpath,'hist_clusters_050.mat'));
load(fullfile(outpath,'hist_samples.mat'));
load('A:\data\vq_histograms\normals_hists_set.mat',...
	'xy_set','r_max','scl');

% preload all N images
trnpath = 'A:\data\orientation_maps\rf\2004_screening_processed\normals';
trndir = dir(fullfile(trnpath,'*_class.mat'));
if ~exist('trn_images','var') || isempty(trn_images)
	for i = 1:length(train_imgs)
		filename = fullfile(trnpath,trndir(train_imgs(i)).name);
		[p,basename,e] = fileparts(trndir(train_imgs(i)).name);
		basename = basename(1:6);

		% load orientation map and mask 
		orimap = load_uint8(filename,true);
		load(['A:\data\masks\2004_screening\normals\',basename,'_mask.mat']); 

		% make sure they're the same height (they aren't always) and
		% subsample
		imsz = min([size(orimap); size(mask)]);
		orimap = orimap(1:scl:imsz(1),1:scl:imsz(2));
		
		% flip left view, and compute line strength and orienation
		if basename(4)=='L'
			orimap = conj(orimap(:,end:-1:1)); % flip l-r and angle
		end

		trn_images{i} = orimap;
	end
end

% reshape train_inds to column order
train_inds = train_inds';

% for each cluster...
for c = 1:size(centres,2)
	% get index of point and image for all associated histograms
	h_inds = find((idx==c));
	[pt,img] = ind2sub(size(train_inds),h_inds);
	
	% get associated histograms themselves
	hists_c = hists(:,:,h_inds);
	hists_c = reshape(hists_c,[numel(hists_c(:,:,1)),size(hists_c,3)]);

	% find distance of each example to cluster centre
	d = zeros(1,size(hists_c,2));
	for i = 1:size(hists_c,2)
		d(i) = hist_dist(hists_c(:,i),centres(:,c));
	end

	% sort into ascending order
	[d,d_inds] = sort(d);

	% make folder for images
	imgpath = fullfile(outpath,sprintf('centre_%03d',c));
	if ~exist(imgpath,'dir'), mkdir(imgpath); end
	
	% get indices of image and point in ascending order
	img = img(d_inds); pt = pt(d_inds);
	for i = 1:length(d)		
		trn_img = train_imgs(img(i));
		trn_pt = train_inds(pt(i),img(i));

		x = xy_set{trn_img}(trn_pt,1);
		y = xy_set{trn_img}(trn_pt,2);
		patch = trn_images{img(i)}(y-r_max:y+r_max,x-r_max:x+r_max);
								
		patch_rgb = ori2rgb(patch);
		imwrite(uint8(255*patch_rgb),...
				fullfile(imgpath,sprintf('sample_%03d.png',i)));
	end
end



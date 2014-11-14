clear;

%% define log-polar map
% logpolar map with ncells_r radial cells and ncells_t angular cells, 
% applied for a disc with inner radius r_min and outer radius r_max
r_min = 10; r_max = 50; ncells_r = 3; ncells_t = 6;
[lp_map,lp_dist] = logpolar_map(r_min,r_max,ncells_r,ncells_t);

% update lp_dist for multiple orientations
nbins_t = 3; % number of bins we want for orientation
[xx,yy] = meshgrid(1:nbins_t,1:nbins_t);
t_dist = min(abs(xx-yy),nbins_t-abs(xx-yy));
lp_dist =	kron(ones(nbins_t,nbins_t),lp_dist) + ...
			kron(t_dist,ones(size(lp_dist)));

%% get histograms
dataroot	= 'a:\data\';
for subgroup = {'normals','abnormals'}
	subgroup = subgroup{1};

	datapath	= [dataroot,'orientation_maps\rf\2004_screening_processed\',subgroup,'\'];
	ext			= '_class.mat';

	% subsampling scale (4 = quarter size)
	scl = 4;

	% get list of all images in datapath
	datadir	= dir([datapath,'*',ext]);

	%% read images
	n_images = min(inf,length(datadir)); % inf = all images
	n_per_pair = 1000;

	hists_set = cell(1,n_images);
	xy_set = cell(1,n_images);
	tb = timebar('limit',n_images,'title','Sampling histograms');
	for iimg = 1:n_images
		[p,basename,e] = fileparts(datadir(iimg).name);
		basename = basename(1:6);

		% load orientation matp
		filename = [datapath,basename,ext];
		orimap = load_uint8(filename,true);
		
		% load mask
		load([dataroot,'masks\2004_screening\',subgroup,'\',basename,'_mask.mat']); 

		% make sure they're the same height (they aren't always)
		imsz = min([size(orimap); size(mask)]);
		orimap = orimap(1:imsz(1),1:imsz(2));
		mask = mask(1:imsz(1),1:imsz(2));

		% subsample
		orimap = orimap(1:scl:end,1:scl:end);
		mask = mask(1:scl:end,1:scl:end);
		
		% flip left view, and compute line strength and orienation
		if basename(4)=='L'
			orimap = conj(orimap(:,end:-1:1)); % flip l-r and angle
			mask = mask(:,end:-1:1); % flip l-r
		end
		
		% get angle
		ori_ang = angle(orimap);
		
		% avoid the edges
		mask([1:r_max,end-r_max+1:end],:) = 0;
		mask(:,[1:r_max,end-r_max+1:end]) = 0;

		% ditch values that are outside the mask
		orimap(~mask) = NaN;
		ori_ang(~mask) = NaN;

		% randomize all points within mask
		[y,x] = find(mask);
		inds = [x,y]; 
		inds = inds(randperm(length(x)),:);

		% pick points and compute histograms
		hists = zeros(ncells_r*ncells_t,nbins_t,n_per_pair);
		for i_pt = 1:n_per_pair
			x = inds(i_pt,1); y = inds(i_pt,2);
			ori_ref = ori_ang(y-r_max:y+r_max,x-r_max:x+r_max);
			hists(:,:,i_pt) = ...
				compute_lp_hist(ori_ref,lp_map,nbins_t,ncells_r*ncells_t);
		end
		
		hists_set{iimg} = hists;
		xy_set{iimg} = inds(1:n_per_pair,:);
		
		timebar(tb,'advance');
	end	
	timebar(tb,'close');

	% save histograms
% 	outpath = 'U:\projects\mammography\data\histograms\';
	outpath = fullfile(dataroot,'vq_histograms');
	filename = sprintf('%s_hists_set.mat',subgroup);
	save(fullfile(outpath,filename),...
		'r_min','r_max','ncells_r','ncells_t','nbins_t',...
		'scl','xy_set','hists_set','lp_dist');
end



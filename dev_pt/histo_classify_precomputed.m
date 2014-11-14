% clear;

% load cluster centres
load(fullfile(outpath,'hist_clusters_050.mat'),...
	'centres','idx');
K = size(centres,2);

% get indices of training images
clear train_imgs;
load(fullfile(outpath,'hist_samples.mat'),'train_imgs');

% define file stuff
dataroot	= 'a:\data\';

for subgroups = {'normals','abnormals'}
	subgroup = subgroups{1};
	
	% get list of all images in datapath
	datapath	= [dataroot,'orientation_maps\rf\2004_screening_processed\',subgroup,'\'];
	ext			= '_class.mat';
	datadir	= dir([datapath,'*',ext]);

	% load precomputed histograms
	hists_set = {};
	load(fullfile(dataroot,['vq_histograms\',subgroup,'_hists_set.mat']),...
		'hists_set');

	% number of histograms per pair
	n_per_pair = min(inf,size(hists_set{1},3));
	
	tb = timebar('title',sprintf('Matching histograms (%s)',subgroup));
	vq_dists = [];
	for Lind = 1:length(datadir)
		[p,basename,e] = fileparts(datadir(Lind).name);
		basename = basename(1:6);

		% only need to deal with each pair once
		if (basename(4)=='R'), continue; end
		basenameL = basename;
		basenameR = basename; basenameR(4) = 'R';
		
		% find index of right image
		Rind = find(strcmp([basenameR,ext],{datadir(:).name}));
		
		% if one of the pair is a training image then skip
		if strcmp(subgroup,'normals') && ...
		   (any(train_imgs==Lind) || any(train_imgs==Rind))
			fprintf('Excluding %s and %s\n',basenameL,basenameR);
			continue;
		end
		
		% pick points and compute histograms
		distL = zeros(K,n_per_pair);
		hinds = randperm(size(hists_set{Lind},3));
		for i_pt = 1:n_per_pair
			for i_centre = 1:K
				distL(i_centre,i_pt) = hist_dist(hists_set{Lind}(:,:,hinds(i_pt)),...
												 centres(:,i_centre),'L1');
			end
		end
		[ignore,matchesL] = min(distL);
		vq_histL = hist(matchesL,1:K);

		% do same for right image
		distR = zeros(K,n_per_pair);
		hinds = randperm(size(hists_set{Rind},3));
		for i_pt = 1:n_per_pair
			for i_centre = 1:K
				distR(i_centre,i_pt) = hist_dist(hists_set{Rind}(:,:,hinds(i_pt)),...
												 centres(:,i_centre),'L1');
			end
		end
		[ignore,matchesR] = min(distR);
		vq_histR = hist(matchesR,1:K);

		% compute distance between vector quantized histograms
		vq_dists(end+1) = hist_dist(vq_histL,vq_histR);
		timebar(tb,'advance');
	end	
	timebar(tb,'close');

	% save distances
	varname = sprintf('vq_dists_%s',subgroup);
	eval([varname,' = vq_dists;']);
	save(fullfile(outpath,[varname,'.mat']),varname);
end % for subgroup



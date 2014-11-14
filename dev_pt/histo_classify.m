% clear;

% load cluster centres
% outpath = 'U:\projects\mammography\data\histograms\';
load(fullfile(outpath,['hist_clusters_',imgview,'_050.mat']),...
	'centres','idx');

K = size(centres,2);

% define file stuff
dataroot	= 'a:\data\';

for subgroups = {'normals','abnormals'}
	subgroup = subgroups{1};
	
	datapath	= [dataroot,'orientation_maps\rf\2004_screening_processed\',subgroup,'\'];
	ext			= '_class.mat';

	% scl = 1; sclpath = '';
	% scl = 2; sclpath = 'halfsize\';
	scl = 4; sclpath = 'qtrsize\';

	datadir		= dir([datapath,'*L',imgview,ext]);

%% define log-polar map
	if ~exist('lp_dist','var')
		% logpolar map with 2 radial bins and 4 angular bins, applied for a disc
		% with inner radius 5 and outer radius r_max (50 at the moment)
		r_min = 10; r_max = 50; nbins_r = 3; nbins_t = 6;
		[lp_map,lp_dist] = logpolar_map(r_min,r_max,nbins_r,nbins_t);
		nbins_lp = nbins_r*nbins_t;
		r_rng = round(logspace(log10(r_min),log10(r_max),nbins_r+1));
		t_rng = linspace(-pi,pi,nbins_t+1); t_rng(end) = [];

		% update lp_dist for multiple orientations
		nbins_t = 3; % number of bins we want for orientation
		[xx,yy] = meshgrid(1:nbins_t,1:nbins_t);
		t_dist = min(abs(xx-yy),nbins_t-abs(xx-yy));
		lp_dist =	kron(ones(nbins_t,nbins_t),lp_dist) + ...
					kron(t_dist,ones(size(lp_dist)));
	end

%% read images
	n_images = min(inf,length(datadir));
	n_per_pair = 500;
	img_dists = nan(1,n_images);

	hists = {};
	for iimg = 1:n_images
		% choose an image at random
		[p,basename,e] = fileparts(datadir(iimg).name);
		basename = basename(1:6);
		basenameL = basename; basenameL(4) = 'L';
		basenameR = basename; basenameR(4) = 'R';
		
		% skip any images that were used for picking cluster centres
		if	strcmp(subgroup,'normals') && any(strcmp(basenameL,train_imgs))
			fprintf('Excluding %s\n',basename);
			continue;			
		end

		% load images and masks
		lml_name = [datapath,basenameL,ext];
			lml = load_uint8(lml_name,true);
			load([dataroot,'masks\2004_screening\',subgroup,'\',basenameL,'_mask.mat']); 
			maskL = mask;
		rml_name = [datapath,basenameR,ext];
			rml = load_uint8(rml_name,true);
			load([dataroot,'masks\2004_screening\',subgroup,'\',basenameR,'_mask.mat']);
			maskR = mask;

		% make sure they're the same height (they aren't always)
		imsz = min([size(lml); size(rml); size(maskL); size(maskR)]);
		lml = lml(1:imsz(1),1:imsz(2));
		rml = rml(1:imsz(1),1:imsz(2));
		maskL = maskL(1:imsz(1),1:imsz(2));
		maskR = maskR(1:imsz(1),1:imsz(2));

		% flip left view, and compute line strength and orienation
		lml_map = abs(lml(1:scl:end,end:-scl:1,1)); % flipped l-r
			lml = angle(conj(lml(1:scl:end,end:-scl:1,1))); % flip angle, too
			maskL = maskL(1:scl:end,end:-scl:1); 
		rml_map = abs(rml(1:scl:end,1:scl:end,1));
			rml = angle(rml(1:scl:end,1:scl:end,1));
			maskR = maskR(1:scl:end,1:scl:end); 

		% ditch values that are outside of the mask
		lml(~maskL) = NaN; rml(~maskR) = NaN;

		% trim off edges
		trim = 20;
		lml = lml(trim+1:end-trim,trim+1:end-trim);
			lml_map = lml_map(trim+1:end-trim,trim+1:end-trim);
			maskL = maskL(trim+1:end-trim,trim+1:end-trim);
		rml = rml(trim+1:end-trim,trim+1:end-trim);
			rml_map = rml_map(trim+1:end-trim,trim+1:end-trim);
			maskR = maskR(trim+1:end-trim,trim+1:end-trim);

		% avoid the edges
		maskL([1:r_max,end-r_max+1:end],:) = 0;
		maskL(:,[1:r_max,end-r_max+1:end]) = 0;
		maskR([1:r_max,end-r_max+1:end],:) = 0;
		maskR(:,[1:r_max,end-r_max+1:end]) = 0;

		% build a colour image
		lml_rgb = (lml+pi/2)/pi; lml_rgb(:,:,2) = 1; lml_rgb(:,:,3) = lml_map; lml_rgb = hsv2rgb(lml_rgb);
		rml_rgb = (rml+pi/2)/pi; rml_rgb(:,:,2) = 1; rml_rgb(:,:,3) = rml_map; rml_rgb = hsv2rgb(rml_rgb);

	% 	% display left and right
	% 	figure(1); clf;
	% 	subplot(2,2,1); image(lml_rgb); title(''); axis('off','image');
	% 		hold on;
	% 		for r = r_rng
	% 			rectangle(	'position',[r_max-r,r_max-r,2*r,2*r],...
	% 						'curvature',[1,1],'edgecolor','y');
	% 		end
	% 		t_x = [r_min*cos(t_rng); r_max*cos(t_rng)];
	% 		t_y = [r_min*sin(t_rng); r_max*sin(t_rng)];
	% 		plot(r_max+t_x,r_max+t_y,'y-');
	% 	subplot(2,2,2); image(rml_rgb); title(''); axis('off','image');

		% choose points for sampling
		[y,x] = find(maskL); indsL = [x,y]; indsL = indsL(randperm(length(x)),:);
		[y,x] = find(maskR); indsR = [x,y]; indsR = indsR(randperm(length(x)),:);

		% pick points and compute histograms
		matchesL = zeros(1,K); matchesR = zeros(1,K);
		for i_pt = 1:n_per_pair
			% compute histogram for random point in left image and 
			% match to cluster centre
			pt_ind = ceil(rand*length(indsL));
			x = indsL(pt_ind,1); y = indsL(pt_ind,2);
			ori_ref = lml(y-r_max:y+r_max,x-r_max:x+r_max);
			histL = compute_lp_hist(ori_ref,lp_map,nbins_t,nbins_lp);
			for i_centre = 1:K
				dst(i_centre) = hist_dist(histL(:),centres(:,i_centre),'L1');
			end
			[ignore,ind] = min(dst);
			matchesL(ind) = matchesL(ind)+1;

			% do same for right image
			pt_ind = ceil(rand*length(indsR));
			x = indsR(pt_ind,1); y = indsR(pt_ind,2);
			ori_ref = rml(y-r_max:y+r_max,x-r_max:x+r_max);
			histR = compute_lp_hist(ori_ref,lp_map,nbins_t,nbins_lp);
			for i_centre = 1:K
				dst(i_centre) = hist_dist(histR(:),centres(:,i_centre),'L1');
			end
			[ignore,ind] = min(dst);
			matchesR(ind) = matchesR(ind)+1;

			% save these for later
			hists{iimg,i_pt,1} = uint16(histL);
			hists{iimg,i_pt,2} = uint16(histR);
		end

		% compute distance between vector quantized histograms
		img_dists(iimg) = hist_dist(matchesL,matchesR);
	end	

	% save distances
	varname = sprintf('dists_%s_%s',subgroup,imgview);
	eval([varname,' = img_dists;']);
	save(fullfile(outpath,[varname,'.mat']),varname);

	% save sampled histograms
	varname = sprintf('hists_%s_%s',subgroup,imgview);
	eval([varname,' = hists;']);
	save(fullfile(outpath,[varname,'.mat']),varname);
	
end % for subgroup



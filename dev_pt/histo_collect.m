% clear;

dataroot	= 'a:\data\';
subgroup	= 'normals';
datapath	= [dataroot,'orientation_maps\rf\2004_screening_processed\',subgroup,'\'];
ext			= '_class.mat';

% scl = 1; sclpath = '';
% scl = 2; sclpath = 'halfsize\';
scl = 4; sclpath = 'qtrsize\';

datadir	= dir([datapath,'*L',imgview,ext]); % choose from left side only
datadir	= datadir(randperm(length(datadir))); % randomize

%% read images
n_images = min(5,length(datadir));
n_per_pair = 50;

train_imgs = cell(1,n_images);
hists = zeros(0,0,0); % need to make it 3D from the start
for iimg = 1:n_images
	% choose an image at random
	[p,basename,e] = fileparts(datadir(iimg).name);
	basename = basename(1:6);
	basenameL = basename; basenameL(4) = 'L';
	basenameR = basename; basenameR(4) = 'R';

	train_imgs{iimg} = basename;

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
	for i_pt = 1:n_per_pair
		img_ind = (rand>0.5)+1; % 1 or 2

		% compute histogram for random point
		switch img_ind,
			case 1,
				pt_ind = ceil(rand*length(indsL));
				x = indsL(pt_ind,1); y = indsL(pt_ind,2);
				ori_ref = lml(y-r_max:y+r_max,x-r_max:x+r_max);
			case 2,
				pt_ind = ceil(rand*length(indsR));
				x = indsR(pt_ind,1); y = indsR(pt_ind,2);
				ori_ref = rml(y-r_max:y+r_max,x-r_max:x+r_max);
		end
		hists(:,:,end+1) = compute_lp_hist(ori_ref,lp_map,nbins_t,nbins_lp);
	end
end	

%% plot 3D histograms on a simplex
if numel(hists(:,:,1))==3
	% reshape and normalize histograms
	hists3 = reshape(hists,[3,size(hists,3)]);
	hists3 = hists3 ./ (ones(3,1)*sum(hists3,1));
	
	% define basis for points
	theta = [0 1 2]*2*pi/3 + pi/2;
	basis = [cos(theta); sin(theta)];
	
	% project to image points
	pts = basis*hists3;

	% plot
	figure(3); clf; hold on; box on;
		plot(basis(1,[1,2,3,1]),basis(2,[1,2,3,1]),'-','color',0.9*[1,1,1]);
		plot(pts(1,:),pts(2,:),'b.');
		axis('equal',[-1,1,-0.6,1.1]);
% 		title(['Distance: ',distname]);
	
% 	figpath = 'U:\projects\mammography\figs\matches\';
% 	graph(3); exportfig([figpath,basename,'_',distname,'_hist3d']);
end

% save data
% outpath = 'U:\projects\mammography\data\histograms\';
save(fullfile(outpath,['hist_samples_',imgview,'.mat']),...
	'hists','lp_dist','train_imgs');



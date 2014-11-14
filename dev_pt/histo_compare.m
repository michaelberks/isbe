% clear;

dataroot	= 'a:\data\';
subgroup	= 'normals';
% subgroup	= 'abnormals';
datapath	= [dataroot,'mammograms\2004_screening\',subgroup,'\'];
datapath	= [dataroot,'orientation_maps\rf\2004_screening_processed\',subgroup,'\'];
	ext		= '_class.mat';

% scl = 1; sclpath = '';
% scl = 2; sclpath = 'halfsize\';
scl = 4; sclpath = 'qtrsize\';

% datapath	= [datapath,'png\',sclpath];
% 	datadir		= dir([datapath,'*LCC.png']);
datadir		= dir([datapath,'*LCC',ext]);

% choose an image at random
ip = ceil(rand*length(datadir));
[p,basename,e] = fileparts(datadir(ip).name);
basename = basename(1:6);
basename = '044LCC';
display(basename);

% load left and right images
if ~exist('lml','var')
	lml_name	= [datapath,basename(1:3),'LCC',ext];
	lml = load_uint8(lml_name,true);
	rml_name	= [datapath,basename(1:3),'RCC',ext];
	rml = load_uint8(rml_name,true);

	% make sure they're the same height (they aren't always)
	rml = rml(1:size(lml,1),:);

	% flip left view, and compute line strength and orienation
	lml_map = abs(lml(1:scl:end,end:-scl:1,1)); % flipped l-r
	lml = angle(conj(lml(1:scl:end,end:-scl:1,1))); % flip angle, too
	rml_map = abs(rml(1:scl:end,1:scl:end,1));
	rml = angle(rml(1:scl:end,1:scl:end,1));
	
	% load left and right masks
	load([dataroot,'masks\2004_screening\',subgroup,'\',basename(1:3),'LCC_mask.mat']); 
		maskl = mask(1:scl:end,end:-scl:1); lml(~maskl) = NaN;
	load([dataroot,'masks\2004_screening\',subgroup,'\',basename(1:3),'RCC_mask.mat']);
		maskr = mask(1:scl:end,1:scl:end); rml(~maskr) = NaN;
		maskr = maskr(1:size(maskl,1),:);

	% trim off edges
	trim = 20;
	lml = lml(trim+1:end-trim,trim+1:end-trim);
	lml_map = lml_map(trim+1:end-trim,trim+1:end-trim);
	maskl = maskl(trim+1:end-trim,trim+1:end-trim);
	rml = rml(trim+1:end-trim,trim+1:end-trim);
	rml_map = rml_map(trim+1:end-trim,trim+1:end-trim);
	maskr = maskr(trim+1:end-trim,trim+1:end-trim);

	% build a colour image
	lml_rgb = (lml+pi/2)/pi; lml_rgb(:,:,2) = 1; lml_rgb(:,:,3) = lml_map; lml_rgb = hsv2rgb(lml_rgb);
	rml_rgb = (rml+pi/2)/pi; rml_rgb(:,:,2) = 1; rml_rgb(:,:,3) = rml_map; rml_rgb = hsv2rgb(rml_rgb);
end

% display left and right
figure(1); clf;
	subplot(2,2,1); image(lml_rgb); title(''); axis('off','image');
	subplot(2,2,2); image(rml_rgb); title(''); axis('off','image');

%% define log-polar map
% logpolar map with 2 radial bins and 4 angular bins, applied for a disc
% with inner radius 5 and outer radius r_max (50 at the moment)
r_min = 5; r_max = 50; nbins_r = 1; nbins_t = 1;
[lp_map,lp_dist] = logpolar_map(r_min,r_max,nbins_r,nbins_t);
nbins_lp = max(lp_map(:));

% draw it on the image
r_rng = round(logspace(log10(r_min),log10(r_max),nbins_r+1));
t_rng = linspace(-pi,pi,nbins_t+1); t_rng(end) = [];
figure(1); subplot(2,2,1); hold on;
	for r = r_rng
		rectangle(	'position',[r_max-r,r_max-r,2*r,2*r],...
					'curvature',[1,1],'edgecolor','y');
	end
	t_x = [r_min*cos(t_rng); r_max*cos(t_rng)];
	t_y = [r_min*sin(t_rng); r_max*sin(t_rng)];
	plot(r_max+t_x,r_max+t_y,'y-');
	axis('equal');

% update lp_dist for multiple orientations
nbins_t = 3; % number of bins we want for orientation
[xx,yy] = meshgrid(1:nbins_t,1:nbins_t);
t_dist = min(abs(xx-yy),nbins_t-abs(xx-yy));
lp_dist =	kron(ones(nbins_t,nbins_t),lp_dist) + ...
			kron(t_dist,ones(size(lp_dist)));

% get reference point for histogram computation
% subplot(2,2,1); [x,y,button] = ginput(1);

% how big a step when computing histogram signatures
r_step = round(r_max/2);

%% define reference positions
% if button==1
% 	xy_ref = round([x y]);
% else
	% ignore selected point and define a grid
	[xx,yy] = meshgrid( r_max+1:r_step:size(lml,2)-r_max-1,...
						r_max+1:r_step:size(lml,1)-r_max-1 );
	xy_ref = [xx(:) yy(:)];
% end
for i_ref = size(xy_ref,1):-1:1
	if ~maskl(xy_ref(i_ref,2),xy_ref(i_ref,1)), xy_ref(i_ref,:) = []; end
end

%% define target positions
[xx,yy] = meshgrid( r_max+1:r_step:size(rml,2)-r_max-1,...
					r_max+1:r_step:size(rml,1)-r_max-1 );
xy_tgt = [xx(:) yy(:)];
for i_tgt = size(xy_tgt,1):-1:1
	if ~maskr(xy_tgt(i_tgt,2),xy_tgt(i_tgt,1)), xy_tgt(i_tgt,:) = []; end
end

%% compute histograms for reference and target positions
tic;
hists_ref = []; hists_tgt = [];
for i_ref = 1:size(xy_ref,1)
	x = xy_ref(i_ref,1); y = xy_ref(i_ref,2);
	ori_ref = lml(y-r_max:y+r_max,x-r_max:x+r_max);
	hists_ref(:,:,i_ref) = compute_lp_hist(ori_ref,lp_map,nbins_t,nbins_lp);
end
for i_tgt = 1:size(xy_tgt,1)
	x = xy_tgt(i_tgt,1); y = xy_tgt(i_tgt,2);
	ori_tgt = rml(y-r_max:y+r_max,x-r_max:x+r_max);
	hists_tgt(:,:,i_tgt) = compute_lp_hist(ori_tgt,lp_map,nbins_t,nbins_lp);
end
toc;

%% compute distances between histograms
% distname = 'bhattacharyya';
% distname = 'L1';
% distname = 'L2';
% distname = 'chisquared';
% distname = 'jeffrey';
% distname = 'emd'; dist_params = {lp_dist};
distname = 'emdn'; dist_params = {lp_dist};
% distname = 'match';

tic;
dists = nan(size(xy_ref,1),size(xy_tgt,1));
for i = 1:size(dists,1)
	for j = 1:size(dists,2)
		dists(i,j) = hist_dist(hists_ref(:,:,i),hists_tgt(:,:,j),distname);
	end
end
toc;

%% display distance (or similarity) image
% create a match (i.e. similarity = 1-dist) image
matches = zeros(size(lml));
matches(sub2ind(size(matches),xy_tgt(:,2),xy_tgt(:,1))) = ...
	1 - dists(end,:);

% smooth the matches map (i.e. Parzen windows approximation) and display
matches(isnan(matches)) = 0;
g = linspace(-3,3,2*r_max+1); g = exp(-g.*g);
matches_smooth = conv2(g,g,matches,'same');
figure(1); colormap([0 0 0; jet(255)]);
	subplot(2,2,3); imagesc(matches); axis('off','image');
	subplot(2,2,4); imagesc(matches_smooth); axis('off','image');
	
	
%% show matches between sides

[ignore,match_inds] = min(dists,[],2);
if ~all(isnan(match_inds))
	match_pts = [xy_ref xy_tgt(match_inds,:)];
	match_pts(:,3) = match_pts(:,3)+size(lml,2);
	
	figure(2); clf; hold on; colormap(gray(256));
		imagesc(uint8(255*[lml_map rml_map])); axis('off','image');
		plot(match_pts(:,[1,3])',match_pts(:,[2,4])','g.-');
		for r = r_rng
			rectangle(	'position',[r_max-r,r_max-r,2*r,2*r],...
						'curvature',[1,1],'edgecolor','y');
		end
		t_x = [r_min*cos(t_rng); r_max*cos(t_rng)];
		t_y = [r_min*sin(t_rng); r_max*sin(t_rng)];
		plot(r_max+t_x,r_max+t_y,'y-');
	
	figpath = 'U:\projects\mammography\figs\matches\';
	exportfig([figpath,basename,'_',distname,'_matches']);
end

%% if histograms are 3D then plot on the simplex
if prod(size(hists_ref(:,:,1)))==3
	% reshape histograms
	hists_ref = reshape(hists_ref,[3,size(hists_ref,3)]);
	hists_tgt = reshape(hists_tgt,[3,size(hists_tgt,3)]);
	
	% normalize
	hists_ref = hists_ref ./ (ones(3,1)*sum(hists_ref,1));
	hists_tgt = hists_tgt ./ (ones(3,1)*sum(hists_tgt,1));
	
	% define basis for points
	theta = [0 1 2]*2*pi/3 + pi/2;
	basis = [cos(theta); sin(theta)];
	
	% project to image points
	ref_pts = basis*hists_ref;
	tgt_pts = basis*hists_tgt;

	% plot
	figure(3); clf; hold on; box on;
		plot(basis(1,[1,2,3,1]),basis(2,[1,2,3,1]),'-','color',0.9*[1,1,1]);
		plot([ref_pts(1,:);tgt_pts(1,match_inds)],...
			 [ref_pts(2,:);tgt_pts(2,match_inds)],'g-');
		plot(ref_pts(1,:),ref_pts(2,:),'r.');
		plot(tgt_pts(1,:),tgt_pts(2,:),'b.');
		axis('equal',[-1,1,-0.6,1.1]);
		title(['Distance: ',distname]);
	
	figpath = 'U:\projects\mammography\figs\matches\';
	graph(3); exportfig([figpath,basename,'_',distname,'_hist3d']);
end


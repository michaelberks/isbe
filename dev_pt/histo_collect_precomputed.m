% clear;

dataroot	= 'a:\data\vq_histograms\';
load(fullfile(dataroot,'normals_hists_set.mat'));

%% sample point histograms
n_images = min(5,length(hists_set));
n_per_pair = min(100,size(hists_set{1},3));

train_imgs = randperm(length(hists_set));
train_imgs = train_imgs(1:n_images);

hists = [];
train_inds = [];
for iimg = train_imgs
	inds = randperm(size(hists_set{iimg},3));
	inds = inds(1:n_per_pair);
	
	% concatenate existing histograms with randomly selected ones from
	% hists_set
	hists = cat(3,hists,hists_set{iimg}(:,:,inds));

	% keep record of which histograms were used
	train_inds(end+1,:) = inds;
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
save(fullfile(outpath,'hist_samples.mat'),...
	'hists','lp_dist','train_imgs','train_inds');



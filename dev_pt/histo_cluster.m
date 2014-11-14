% clear;

% load data
% outpath = 'U:\projects\mammography\data\histograms\';
load(fullfile(outpath,'hist_samples.mat'));

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

% cluster histograms (L1 distance by default)
K = 50;
[idx,centres] = histo_kmeans(hists,K);

% save data
% outpath = 'U:\projects\mammography\data\histograms\';
filename = sprintf('hist_clusters_%03d.mat',K);
save(fullfile(outpath,filename),'centres','idx');


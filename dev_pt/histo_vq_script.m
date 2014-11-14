clear;

outpath = 'U:\projects\mammography\data\histograms\';
outpath = fullfile(outpath,datestr(now,'yyyy-mm-dd HH.MM.SS'));
mkdir(outpath);

%% define log-polar map
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

%% run miniscripts with precomputed data
fprintf('Sampling precomputed histograms...\n');
histo_collect_precomputed;
fprintf('Clustering...\n');
histo_cluster;
fprintf('Creating gallery...\n');
histo_gallery;
fprintf('Matching...\n');
histo_classify_precomputed;
fprintf('Comparing...\n');
histo_compare_vq;

return

%% run mini scripts
imgview = 'CC';
	histo_collect;
	histo_cluster;
	histo_classify;

imgview = 'ML';
	histo_collect;
	histo_cluster;
	histo_classify;

histo_compare_vq;

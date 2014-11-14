function lp_hist = compute_lp_hist(ori_map,lp_map,nbins_t,nbins_lp)

% debug mode
f_debug = false;
if nargin==0 && nargout==0
	f_debug = true;	clc;
	lp_map = logpolar_map;
	ori_map = randn(size(lp_map))/3 * pi/2;
		ori_map = min(max(ori_map,-pi/2),pi/2);
end

% default parameter values
if nargin<3, nbins_t = 8; end
if nargin<4, nbins_lp = max(lp_map(:)); end

% precompute indices for histogram map
lp_inds = cell(1,nbins_lp);
for ib = 1:nbins_lp, lp_inds{ib} = find(lp_map==ib); end

% compute histogram
t_rng = linspace(-pi/2,pi/2,nbins_t+1);
lp_hist = zeros(nbins_lp,nbins_t+1);
for ib = 1:nbins_lp
	lp_hist(ib,:) = histc(ori_map(lp_inds{ib}),t_rng);
end
% combine last two columns from histc
lp_hist = [lp_hist(:,1:end-2) sum(lp_hist(:,end-1:end),2)];


if f_debug
	% precompute binned orientation map
	t_rng = linspace(-pi/2,pi/2,nbins_t+1); t_rng(end) = [];
	t_map = nan(size(ori_map));
	for it = 1:nbins_t, t_map(ori_map>=t_rng(it)) = it; end

	lp_hist2 = zeros(nbins_lp,nbins_t);
	for it = 1:nbins_t
		ori_mask = (t_map==it);
		for ib = 1:nbins_lp
			lp_hist2(ib,it) = sum(ori_mask(lp_inds{ib}));
		
			rgb(:,:,1) = ori_mask;
			rgb(:,:,2) = (lp_map==ib);
			rgb(:,:,3) = 0;

			figure(1); clf;
			image(rgb); axis('off','image');
		end
	end
	
	d = lp_hist-lp_hist2;
	sum(abs(d(:)))

	clear;
end




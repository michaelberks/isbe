function [lp_map,lp_dist,r_map,t_map] = logpolar_map(r_min,r_max,nbins_r,nbins_t)

% debug mode
f_debug = false;
if nargin==0 && nargout==0
	f_debug = true; clc;
end

% default parameter values
if nargin<1, r_min = 10; end
if nargin<2, r_max = 50; end
if nargin<3, nbins_r = 5; end
if nargin<4, nbins_t = 12; end

% compute radius and angle values
[xx,yy] = meshgrid(-r_max:r_max,-r_max:r_max);
rr = sqrt(xx.*xx+yy.*yy);
tt = atan2(yy,xx);

% radius map (log scale)
r_rng = round(logspace(log10(r_min),log10(r_max),nbins_r+1));
r_map = nan(size(rr));
for ir = 1:nbins_r, r_map(rr>r_rng(ir)) = ir; end
r_map(rr>r_rng(nbins_r+1)) = nan;

% angle map (linear)
t_rng = linspace(-pi,pi,nbins_t+1); t_rng(end) = [];
t_map = nan(size(tt));
for it = 1:nbins_t, t_map(tt>t_rng(it)) = it; end

% complete map of bins
lp_map = (r_map-1)*nbins_t + t_map;

% compute 'distance' between bins if asked for
if nargout>1 || f_debug
% 	t_rng_c = exp(sqrt(-1)*t_rng);
% 	[xx,yy] = meshgrid(t_rng_c,t_rng_c);
% 	d_ang_c = angle(xx.*conj(yy));
% 	d_ang = nbins_t * abs(d_ang_c)/(2*pi);
	
	[xx,yy] = meshgrid(1:nbins_t,1:nbins_t);
	t_dist = min(abs(xx-yy),nbins_t-abs(xx-yy));
	
	[xx,yy] = meshgrid(1:nbins_r,1:nbins_r);
	r_dist = abs(xx-yy);
	
 	lp_dist = kron(ones(nbins_r,nbins_r),t_dist) + ...
			  kron(r_dist,ones(nbins_t,nbins_t));
end

if f_debug
	figure(1); clf; colormap([0 0 0; jet(nbins_t*nbins_r)]);
	subplot(2,2,1); image(uint16(nbins_t*r_map+1)); axis('off','image');
	subplot(2,2,2); image(uint16(nbins_r*t_map+1)); axis('off','image');
	subplot(2,2,3); image(uint16(lp_map+1)); axis('off','image');
	
	img = lp_map;
	for j = 1:size(lp_dist,1), img(lp_map==j) = lp_dist(j,1); end
	subplot(2,2,4); imagesc(img); axis('off','image');
	
	clear;
end

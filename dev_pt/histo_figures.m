clear;

outpath = 'U:\projects\mammography\data\histograms\';
if ~exist(outpath,'dir'), mkdir(outpath); end

%% generate log-polar map image
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

filename = 'A:\data\orientation_maps\rf\2004_screening_processed\normals\001RCC_class.mat';
orimap = load_uint8(filename,true);
orimap = orimap(1:4:end,1:4:end); % subsample
orimap_rgb = ori2rgb(orimap); % convert to RGB image
x = 215; y = 425;

% for whatever reason, pure white lines do not get exported here - hence
% the [0.999,0.999,0.999] colour
figure(1); clf; hold on;
image(orimap_rgb); title(''); axis('off','image','ij');
for r = r_rng
	rectangle(	'position',[x-r,y-r,2*r,2*r],...
				'curvature',[1,1],'edgecolor',0.999*[1,1,1],'linewidth',1);
end
t_x = [r_min*cos(t_rng); r_max*cos(t_rng)];
t_y = [r_min*sin(t_rng); r_max*sin(t_rng)];
plot(x+t_x,y+t_y,'-','color',0.999*[1,1,1],'linewidth',1);
exportfig(fullfile(outpath,'roi_full'));


orimap_roi = orimap_rgb(y-r_max:y+r_max,x-r_max:x+r_max,:);
figure(2); clf; hold on;
image(orimap_roi); title(''); axis('off','image','ij');
for r = r_rng
	rectangle(	'position',[r_max-r+1,r_max-r+1,2*r,2*r],...
				'curvature',[1,1],'edgecolor',0.999*[1,1,1],'linewidth',3);
end
t_x = [r_min*cos(t_rng); r_max*cos(t_rng)];
t_y = [r_min*sin(t_rng); r_max*sin(t_rng)];
plot(r_max+t_x+1,r_max+t_y+1,'-','color',0.999*[1,1,1],'linewidth',3);
exportfig(fullfile(outpath,'roi_zoom'));

%% next


		
		
clc; clear all;

% parameters
s = 1;			% relative size of image with respect to circle
r1 = 0.4;		% inner radius
r2 = 0.9;		% outer radius
N = s*128;		% size of the image
n_lines = 12;	% number of lines

% create grid of radius and theta
img = zeros(N,N);
[xx,yy] = meshgrid(linspace(-s,s,N),linspace(s,-s,N));
radius = sqrt(xx.*xx+yy.*yy);
theta = mod(atan2(yy,xx),pi);

% generate image with n_lines lines per half circle
trng = linspace(0,pi,n_lines+1); 
trng = (trng(1:end-1)+trng(2:end))/2;
for t = trng
	d = theta-t;
	[ignore,indsx] = min(abs(d),[],1);
	inds = sub2ind([N,N],indsx,1:N);
	img(inds) = t;
	[ignore,indsy] = min(abs(d),[],2);
	inds = sub2ind([N,N],1:N,indsy');
	img(inds) = t;
end

% blank out anything outside of the 'doughnut'
img(radius<r1) = 0;
img(radius>r2) = 0;

if (1)
	% use karssemeijer's orientation estimation algorithm
	[ori_map,line_map] = karssemeijer_line_detection(double(img>0),'degrees',0);
	line_map = (line_map>0);
	ori_map = mod(ori_map+pi,pi); ori_map(~line_map) = NaN;
else
	% or use ground truth
	ori_map = img; line_map = (img>0);
end

% show true orientation map and estimated one
figure(1); clf; sbsz = [3,2];
mysubplot(sbsz); imagesc(img); axis('image','off');
mysubplot(sbsz); imagesc(ori_map); axis('image','off');

% baseline: old karssemeijer code (pants)
tic;
[f_i1 f_i2 n_i n_plus N_i p_i K_i] = ...
	karssemeijer_radial_projection(...
		line_map,ori_map,...
		floor(r1*N/2),ceil(r2*N/2),ceil(r1*N/2));
toc;
mysubplot(sbsz); imagesc(f_i1); axis('image','off');
mysubplot(sbsz); imagesc(f_i2); axis('image','off');

profile clear; profile on;

% updated: new karssemeijer code
tic;
[f_i1 f_i2 n_i n_plus N_i p_i K_i] = ...
	karssemeijer_radial_projection_pix(...
		line_map,ori_map,...
		floor(r1*N/2),ceil(r2*N/2)+1,2);
toc;
mysubplot(sbsz); imagesc(f_i1); axis('image','off');
mysubplot(sbsz); imagesc(f_i2); axis('image','off');

% show profile results (if running)
profstat = profile('status');
if strcmp(profstat.ProfilerStatus,'on')
	profile off; profile report;
end

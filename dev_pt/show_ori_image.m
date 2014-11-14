function show_ori_image(ori)
% draw an orientation image using a cyclic colormap to avoid
% discontinuities

if nargin==0
	% generate synthetic radial map
	[xx,yy] = meshgrid(-50:50,-50:50);
	ori = (atan2(-yy,xx)+2*pi)*180/pi; % want rows to go from top to bottom
	ori = mod(ori,180);
end

% number of colours to use
ncols = 32;

% create cyclical colormap
t = linspace(0,2*pi,ncols);
cmap = [(cos(t)+1)/2; (sin(t)+1)/2; zeros(1,ncols)]';
cmap(1,:) = [0 0 0]; % add black for background

% convert to indexed colormap with zeros instead of NaNs
ori = (ncols-1) * (ori/180);

% show image (colour indices start at zero)
imshow(uint8(ori+1),cmap);



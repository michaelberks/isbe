function rgb = ori2rgb(ori,f_norm)

% whether to scale rr to give something closer to hsv2rgb
if nargin<2, f_norm = true; end

f_debug = false;
if nargin==0 && nargout==0
	f_debug = true;
	clc;
	
	% generate synthetic orientation map
	imsz = 400;
	[xx,yy] = meshgrid(-imsz:imsz,-imsz:imsz);
	rr = sqrt(xx.*xx + yy.*yy);
	tt = atan2(yy,xx);
	ori = rr.*exp(sqrt(-1)*tt);

	img(:,:,1) = mod(tt,pi)/pi;
	img(:,:,2) = 1;
	img(:,:,3) = rr/max(rr(:));
	figure(1); clf;
		subplot(2,1,1); image(hsv2rgb(img)); axis('off','image');
		drawnow;
end

tt = mod(angle(ori),pi)/pi;
rr = abs(ori); rr = rr/max(rr(:));
if f_norm, rr = rr.^1.3; end

% convert orientations to colormap indices
cinds = ceil(254*tt)+2; % scale to range 2..256 (1 is reserved for NaNs)
cinds(isnan(cinds)) = 1; % corresponds to black

% copy indices to actual RGB values
cmap = reshape([0 0 0; hsv(255)],[256,1,3]);
rgb = reshape(cmap(cinds,:),[size(cinds),3]);

% scale by strength
rgb = rgb.*repmat(rr,[1,1,3]);

if f_debug
	figure(1);
		subplot(2,1,2); image(rgb); axis('off','image');

	clear;
end
		



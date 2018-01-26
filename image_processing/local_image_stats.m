function [local_mean, local_std, local_x2] = ...
	local_image_stats(image_in, window_size, min_std, mask)
% compute the local mean and standard deviation at each pixel over a square
% window of size window_size. Standard deviation is bounded below by
% min_std. This implementation uses the integral image for quick
% computation

f_debug = (nargin==0 && nargout==0);
if f_debug
	image_in = zeros(5,5); image_in(3,3) = 1;
end

if ~exist('image_in','var'), error('No image supplied'); end
if ~exist('window_size','var'), window_size = 3; end;
if ~exist('min_std','var'), min_std = 0; end;

local_mean = zeros(size(image_in));
local_std = zeros(size(image_in));
local_x2 = zeros(size(image_in));

% check that window_size is odd
if mod(window_size,2)~=1
	error('Window size must be odd');
end

% define size of border needed round the image so we can process every pixel
border = (window_size-1)/2;

if exist('mask','var')
	% apply mask to image
	image_in(~mask) = 0;
	
	% pad mask with zeros
	pad_mask = zeros(size(image_in)+2*border);
	pad_mask(border+1:end-border,border+1:end-border) = mask;

	% compute integral images of maskand pad top/left with zeros
	int_mask = cumsum(cumsum(pad_mask,2));
		int_mask = [zeros(1,size(int_mask,2)); int_mask];
		int_mask = [zeros(size(int_mask,1),1) int_mask];
else
	% precompute number of pixels
	n_pix = window_size*window_size;

	% factor that gives an unbiased estimate for variance
	var_scl	= n_pix/(n_pix-1);
end

% pad image with zeros
padimg = zeros(size(image_in)+2*border);
padimg(border+1:end-border,border+1:end-border) = image_in;

% compute image squared for quick computation of the variance
sqrimg = padimg.*padimg;

% compute integral images and pad top/left with zeros
intimg = cumsum(cumsum(padimg,2));
	intimg = [zeros(1,size(intimg,2)); intimg];
	intimg = [zeros(size(intimg,1),1) intimg];
intsqrimg	= cumsum(cumsum(sqrimg,2));
	intsqrimg = [zeros(1,size(intsqrimg,2)); intsqrimg];
	intsqrimg = [zeros(size(intsqrimg,1),1) intsqrimg];

% for every pixel, p, replace with (p-mean)/std where the mean and std
% are computed over a window of width window_size, centred on p
for y = 1:size(image_in,1)
	for x = 1:size(image_in,2)
		sum_p	= intimg(y+window_size,x+window_size) - ...
				  intimg(y+window_size,x) - ...
				  intimg(y,x+window_size) + ...
				  intimg(y,x);
			  
		sum_p2	= intsqrimg(y+window_size,x+window_size) - ...
				  intsqrimg(y+window_size,x) - ...
				  intsqrimg(y,x+window_size) + ...
				  intsqrimg(y,x);
			  
		if exist('mask','var')
			n_pix	= int_mask(y+window_size,x+window_size) - ...
					  int_mask(y+window_size,x) - ...
					  int_mask(y,x+window_size) + ...
					  int_mask(y,x);
				  
			var_scl = n_pix/(n_pix-1);

			if (n_pix==0)
				local_mean(y,x) = NaN;
				local_std(y,x) = NaN;
				continue;
			end
		end
			
		p_mean	= sum_p/n_pix;
        p2_mean	= sum_p2/n_pix;
		%if p_mean > 255
		%	keyboard
		%end

		local_mean(y,x) = p_mean;
        local_x2(y,x) = var_scl*p2_mean;
		local_std(y,x) = var_scl*(p2_mean - p_mean*p_mean);
	end
end

% turn local variance into unbiased, bounded, local standard deviation
local_std = max(sqrt(local_std),min_std);

if f_debug
	display(image_in);
	display(local_mean);
	display(local_std);
	clear;
end


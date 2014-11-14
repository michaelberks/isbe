function [imgout] = fix_contrast(imgin,method,varargin)
% modify contrast of an image using one of several available method

if nargin<2, method = 'stretch'; end
	
imgout = double(imgin);
switch lower(method)
	case 'stretch',
		imgout = imgout - min(imgout(:));
		imgout = imgout / max(imgout(:));
		
	case 'meanstd',
		imgout = imgout - mean(imgout(:));
		imgout = imgout / std(imgout(:));
		
	case 'localstd',
		wsize = 3;
		if length(varargin)>0, wsize = varargin{1}; end
		minval = 0;
		if length(varargin)>1, minval = varargin{2}; end
				
		% add a border to the image so we can process every pixel
		border = (wsize-1)/2;
		padimg = zeros(size(imgout)+2*border);
		padimg(border+1:end-border,border+1:end-border) = imgout;
		
		% compute image squared for quick computation of the variance
		sqrimg = padimg.*padimg;
		
		% compute integral images and pad top/left with zeros
		intimg = cumsum(cumsum(padimg,2));
			intimg = [zeros(1,size(intimg,2)); intimg];
			intimg = [zeros(size(intimg,1),1) intimg];
		intsqrimg	= cumsum(cumsum(sqrimg,2));
			intsqrimg = [zeros(1,size(intsqrimg,2)); intsqrimg];
			intsqrimg = [zeros(size(intsqrimg,1),1) intsqrimg];

		% precompute number of pixels
		n_pix = wsize*wsize;
		
		% factor that gives an unbiased estimate for variance
		var_scl	= n_pix/(n_pix-1);
		
		% for every pixel, p, replace with (p-mean)/std where the mean and std
		% are computed over a window of width wsize, centred on p
		for y = 1:size(imgout,1)
			for x = 1:size(imgout,2)
				sum_p		= intimg(y+wsize,x+wsize) + intimg(y,x) - ...
									intimg(y+wsize,x) - intimg(y,x+wsize);
				sum_p2	= intsqrimg(y+wsize,x+wsize) + intsqrimg(y,x) - ...
									intsqrimg(y+wsize,x) - intsqrimg(y,x+wsize);
				p_mean	= sum_p/n_pix;
				p_var		= var_scl*(sum_p2/n_pix - p_mean*p_mean);
				p_std		= max(sqrt(p_var),minval);
				
				% note this does not affect subsequent calculations since they are
				% all computed from the integral images which do not change
				imgout(y,x) = (imgout(y,x)-p_mean)/p_std;
			end
		end
		
	otherwise,
		% error?
end

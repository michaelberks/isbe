function [imgout] = normim(imgin, method, varargin)
% NORMIM Modify contrast of an image using one of several available methods
%   'stretch'
%       stretch the image to range [0,1]
%   'stretch_fixed'
%       stretch the image tp range [0,1] where zero input maps to 0.5
%       output
%   'meanstd'
%       normalize with respect to global mean and variance
%   'localstd'
%       normalize with respect to local mean and variance

if nargin<2, method = 'stretch'; end
	
imgout	= double(imgin);
valid	= find(~isnan(imgout));
switch lower(method)
	case 'stretch',
		imgout = imgout - min(imgout(valid));
		if max(imgout(valid))>0
			imgout = imgout / max(imgout(valid));
		end
		
	case 'stretch_fixed',
        maxval = max(abs(imgout(valid)));
		if maxval>0
			imgout = imgout / (2*maxval);
            imgout = imgout + 0.5;
        end

    case 'meanstd',
		imgout = imgout - mean(imgout(valid));
		if std(imgout(valid))>0
			imgout = imgout / std(imgout(valid));
		end
		
	case 'localstd',
		wsize = 3;
		if ~isempty(varargin), wsize = varargin{1}; end
		minval = 0;
		if length(varargin)>1, minval = varargin{2}; end

		% compute local image statistics
		[local_mean, local_std] = local_image_stats(imgout, wsize, minval);
		imgout = (imgout - local_mean) ./ local_std;
		
	otherwise,
		% error?
end

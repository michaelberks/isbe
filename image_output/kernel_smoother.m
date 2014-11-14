function [centres,smoothed_y] = kernel_smoother(x,y,N,kernel,varargin)

if nargin<3, N = 100; end
if nargin<4, kernel = 'gaussian'; end

[sorted_x,inds] = sort(x);
sorted_y = y(inds);

centres = linspace(sorted_x(1),sorted_x(end),N);
smoothed_y = zeros(size(centres));

switch kernel
	case 'gaussian',
		% gaussian kernel density estimate
		
		if length(varargin)<1, w_sigma_scale = 3; end
		
		w_sigma = (centres(2)-centres(1)) * w_sigma_scale;
		for i = 1:N
			normed = (sorted_x - centres(i)) / w_sigma;
			weights = exp(-0.5*normed.*normed);
			smoothed_y(i) = sum(weights.*sorted_y) / sum(weights);
		end

	case 'cummean',
		% cumulative mean
		for i = 1:N
			weights = (sorted_x <= centres(i));
			smoothed_y(i) = sum(weights.*sorted_y) / sum(weights);
		end
		
	otherwise,
		error(['Unknown kernel: ',kernel]);
end

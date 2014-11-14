function [value,prob_im] = objective_function(ref_sample,free_sample,type,nbins,data)
%
% type is passed separately because might want different objective function for affine and NRR
% hence just don't pass params!
%
% ALL objective function written so that small values are GOOD!


% Allowed values of type
type_set = {'sqrdiff' 'absdiff' 'MI' 'normMI' 'MML'};

if(nargin==0)
	value = type_set;
	prob_im = [];
else
	if(length(ref_sample(:))~=length(free_sample(:)))
		value = [];
		prob_im = [];
	else
		free_sample = reshape(free_sample,size(ref_sample));
		%%%%%%%%% START OF SWITCH
		switch type
			case 'sqrdiff'
				prob_im = (ref_sample-free_sample).^2;
				value = mean(prob_im(:));
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			case 'absdiff'
				prob_im = abs(ref_sample-free_sample);
				value = mean(prob_im(:));
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			case 'MI'
				[h1,h2,h12,prob_im] = calc_mutual(ref_sample,free_sample,nbins);
				value = -(h1+h2-h12);
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			case 'normMI'
				[h1,h2,h12,prob_im] = calc_mutual(ref_sample,free_sample,nbins);
				value = -(h1+h2)/h12;
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			case 'MML'
				value = calc_MML(nbins,data);
				prob_im = abs(ref_sample-free_sample);
		%%%%%%%% END OF SWITCH %%%%%%%%%%%%%%
		end
	end
end

%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%
function [h1,h2,h12,im_prob] = calc_mutual(image1,image2,nbins)

% Gray scale histogram
edges = [0:1/nbins:1];
nbins = length(edges) - 1;
dummy = histc(image1(:),edges);
hist1 = dummy(1:nbins);
hist1(nbins) = hist1(nbins) + dummy(nbins+1);
% Normalise the histogram
hist1 = hist1/sum(hist1(:));
where = find(hist1 ~= 0);
h1 = -sum(hist1(where) .* log(hist1(where)));

dummy = histc(image2(:),edges);
hist2 = dummy(1:nbins);
hist2(nbins) = hist2(nbins) + dummy(nbins+1);
% Normlise the histogram
hist2 = hist2/sum(hist2(:));
where = find(hist2 ~= 0);
h2 = -sum(hist2(where) .* log(hist2(where)));

% 2D histogram
hist12 = zeros(nbins);
im_prob = zeros(size(image1));
for i=1:nbins
    index = find((image1 >= edges(i)) & (image1 < edges(i+1)));
    data = image2(index);
    if (~isempty(data))
        dummy = histc(data(:),edges);
	std_im = std(data(:));
	mean_im = mean(data(:));
	if(std_im>0)
		% Im_prob related to inverse of probability!!!
		dummy2 = abs((data-mean_im)/std_im);
		im_prob(index) = 1 - exp(-0.5*(dummy2.^2));
	else
		im_prob(index) = 0;
	end
	if(size(dummy(1:nbins))~=size(hist12(i,:)))
	        hist12(i,:) = (dummy(1:nbins))';    
	else
		hist12(i,:) = dummy(1:nbins);
	end
        hist12(i,nbins) = hist12(i,nbins) + dummy(nbins+1);
    end          
end

hist12 = hist12(:);
hist12 = hist12/sum(hist12(:));
where = find(hist12 ~= 0);
h12 = -sum(hist12(where) .* log(hist12(where)));


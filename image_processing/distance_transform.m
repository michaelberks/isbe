function img = distance_transform(img, cost_x, cost_y)
% Compute distance transform of img

if ~exist('cost_x','var'), cost_x = 1; end;
if ~exist('cost_y','var'), cost_y = cost_x; end;

cost_xy = sqrt(cost_x*cost_x + cost_y*cost_y);

img = padarray(img,[1,1],max(img(:)));

% forward pass
for y = 2:size(img,2)-1
	for x = 2:size(img,1)-1
		img(x,y) = min(min(min(img(x,y),...
                               img(x-1,y-1)+cost_xy),...
                           img(x-1,y)+cost_x),...
                       img(x,y-1)+cost_y);
	end
end

% backward pass
for y = size(img,2)-1:-1:2
	for x = size(img,1)-1:-1:2
		img(x,y) = min(min(min(img(x,y),...
                               img(x+1,y+1)+cost_xy),...
                           img(x+1,y)+cost_x),...
                       img(x,y+1)+cost_y);
	end
end

% trim the edge off
img = img(2:end-1,2:end-1);

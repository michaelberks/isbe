function [pts, widths, inner, outer, img_sz] = ...
    scale_vessel(pts, widths, inner, outer, img_width)
if (nargin==0 && nargout==0), test(); return; end

% Convert data from arbitrary units to pixels

% Add border to limits
xmin = min([inner(:,1); outer(:,1)]);
xmax = max([inner(:,1); outer(:,1)]);
ymin = min([inner(:,2); outer(:,2)]);
ymax = max([inner(:,2); outer(:,2)]);

% Define a border equal to 10% of the height
% (Could be better)
border_sz = 0.1 * (ymax - ymin);

% Compute the scale from units to pixels
scale = (img_width-1) / (xmax-xmin + 2*border_sz);

% Scale everything to pixels
offset = ones(size(pts,1),1) * [-xmin+border_sz -ymin+border_sz];
pts = scale * (pts + offset);
inner = scale * (inner + offset);
outer = scale * (outer + offset);
widths = scale * widths;

% Ensure image sizes are multiples of 8 for pretty videos
img_width = 8*ceil(img_width/8);
img_height = 8*ceil(scale*(ymax-ymin + 2*border_sz)/8);
img_sz = [img_height img_width]; % [rows,columns]


function test()
clc;

pts = generate_path(500);
[widths, dirs, inner, outer] = generate_edges(pts, [], 10);

[pts, widths, inner, outer, imsz] = ...
    scale_vessel(pts, widths, inner, outer, 160);

imsz

figure(1); clf; hold on;
    plot(pts(:,1), pts(:,2), 'b.-');
    plot(pts(1,1), pts(1,2), 'go');
    plot(pts(end,1), pts(end,2), 'ro');
    plot(inner(:,1), inner(:,2), 'r.-');
    plot(outer(:,1), outer(:,2), 'g.-');
    axis('equal','ij',[1,imsz(1),2,imsz(2)]);
    


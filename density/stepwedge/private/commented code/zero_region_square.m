function [ map_zeroed ] = zero_region_square( map,x,y,r )
%ZERO_REGION_SQUARE : Zeros a square region in a 2D array around a point
% 
%
% Inputs:
%			map					map/image to be zeroed [2D array]
%			x 					x position of centre point [#]
%			y					y position of centre point [#]
%			r					half length of square to be zeroed [#]
%
% Outputs:
%			map_zeroed			map with zero region added [2D array]
%
% Example:
%
% Notes:
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com

[Y X] = size(map);
map_zeroed = map;

x1 = x-r; x2 = x+r; y1 = y-r; y2 = y+r;
if x1 < 1; x1=1; end;
if y1 < 1; y1=1; end;
if x2 > X; x2=X; end;
if y2 > Y; y2=Y; end;

map_zeroed(y1:y2,x1:x2)=0;


end


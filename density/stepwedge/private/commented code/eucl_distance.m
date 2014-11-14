function [ distance ] = eucl_distance( x1,y1,x2,y2 )
%calculates euclidean distance between two points

distance = sqrt((x1-x2).^2+(y1-y2).^2);

end


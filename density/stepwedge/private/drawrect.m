function   ok = drawrect(X,Y)
% draws a rectangle
% X is vector containing x coords of opposite corners 
% Y is vector containing y coords of opposite corners

line([X(1) X(1) X(2) X(2) X(1)],[Y(1) Y(2) Y(2) Y(1) Y(1)]);

ok=1;
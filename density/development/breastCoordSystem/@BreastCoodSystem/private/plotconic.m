function plotconic(C)
%PLOTCONIC plot the given conic
%   PLOTCONIV(C) plots the conic C in the figure.

x=1:1000;
y=1:1000;

[X,Y]=meshgrid(x,y);

Z=zeros(size(X));

for idx=1:length(X(:));
  Z(idx)=[X(idx) Y(idx) 1]*C*[X(idx) Y(idx) 1]';
end
contour(X,Y,Z,[0 0],'g-');

function C=parabmin2(x,y,xn,yn)
%PARABMIN2 fits parabola to a given vertex, direction and another point
%   C=PARABMIN2(X,Y,XN,YN) fits the parabola to two points given in vectors
%   X and Y, of which the first elements represent the vectex point.
%   The normal vector pointing to the opening direction of the parabola is
%   (XN,YN). 

%   [C1,C2]=PARABMIN2(...) retuns also the other, usually degenerate, solution.
%
%   Note: sometimes the result may be complex if the point and normal
%   configuration is such. 

n=[xn;yn];
n=n/norm(n);

Rm90=[0 1;-1 0];

xaxis=Rm90*n;
yaxis=n;

R=[xaxis yaxis]';
t=-R*[x(1);y(1)];

D=[R t;0 0 1];

Xp2=D*[x(2);y(2);1];
xp2=Xp2(1);
yp2=Xp2(2);

Cp=[(yp2/xp2^2) 0 0; 0 0 -0.5; 0 -0.5 0];
C=D'*Cp*D;

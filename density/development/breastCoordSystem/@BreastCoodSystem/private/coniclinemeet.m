function [x1,x2]=coniclinemeet(C,l)
%CONICLINEMEET computes the points where a line intersects a conic
%   [X1,X2]=CONICLINEMEET(C,L) returns two points where the conic C 
%   and the line L meet. The solutions are given in the homogeneous
%   form and may be complex. 

%Pick up two points on the line
if l(1)==0 & l(2)==0 %test if l is the l_inf
  x0=[1 0 0];
  xinf=[0 1 0];
else
  x0=[- l(1)*l(3)/(l(1)^2+l(2)^2);-l(2)*l(3)/(l(1)^2+l(2)^2);1];
  xinf=cross([0;0;1],l);
end

a=x0'*C*x0;
b=2*x0'*C*xinf;
c=xinf'*C*xinf;

alpha1=(-b-sqrt(b^2-4*a*c))/(2*a);
alpha2=(-b+sqrt(b^2-4*a*c))/(2*a);

x1=alpha1*x0+xinf; x1=x1/norm(x1);
x2=alpha2*x0+xinf; x2=x2/norm(x2);
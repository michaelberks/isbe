function [x,y,dx,dy]=z2x(z,C,R)

if nargin<3
[U,Lambda]=eig(C(1:2,1:2)); 
if abs(Lambda(2,2))>abs(Lambda(1,1))
    U=U(:,[2 1]);
    Lambda(1,1)=Lambda(2,2);
    Lambda(2,2)=0;
    U(:,2)=-U(:,2);
end
R=U';
else
    U=R';
end
C=[U zeros(2,1); zeros(1,2) 1]'*C*[U zeros(2,1); zeros(1,2) 1];

invC=diag([1,-C(1,1)/C(2,3)/2]);
t=[C(1,3)/C(1,1);(C(1,3)/C(1,1))^2-C(3,3)/C(1,1)];
xv=R'*invC*[z;z^2]-R'*invC*t;
x=xv(1);
y=xv(2);

if nargout>2
    dxv=R'*invC*[1;2*z];
    dx=dxv(1);
    dy=dxv(2);
end
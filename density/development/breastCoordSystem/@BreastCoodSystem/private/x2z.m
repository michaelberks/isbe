function [z,dz]=x2z(x,y,C,xA,yA,R,dR)

if nargin<6
[U,Lambda]=eig(C(1:2,1:2));
if abs(Lambda(2,2))>abs(Lambda(1,1))
    U=U(:,[2 1]);
    Lambda(1,1)=Lambda(2,2);
    Lambda(2,2)=0;
end
if det(U)<0
   U(:,2)=-U(:,2);
end
R=U';
else
    U=R';
end
C=[U zeros(2,1); zeros(1,2) 1]'*C*[U zeros(2,1); zeros(1,2) 1];

z=R(1,1:2)*[x;y]+C(1,3)/C(1,1);

if nargout>1
  dz=[R(1,1:2) dR(1,1:2)*[x-xA;y-yA]];
end


function [s,dsdz,dsdphi]=z2s(z,z0,C,dCdphi,R,dR)

if nargin<5
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
Cnew=[U zeros(2,1); zeros(1,2) 1]'*C*[U zeros(2,1); zeros(1,2) 1];

a=1;
b=(-2*Cnew(2,3)/Cnew(1,1));

A=b/(2*a);

ss=2/b*(intbody(z,A)-intbody(z0,A));
s=abs(ss);

if nargout>1
  dsdz=sign(ss)*2/b*sqrt(A^2+z^2);
  %dR=[R(1,2) -R(2,2); R(1,1) -R(2,1)]; %THE SIGN DOES NOT MATCH WITH THE NUMERICAL DERIVATIVE
  dCnew=[dR' zeros(2,1); zeros(1,2) 0]'*C*[R' zeros(2,1); zeros(1,2) 1]+...
        [R' zeros(2,1); zeros(1,2) 1]'*dCdphi*[R' zeros(2,1); zeros(1,2) 1]+...
        [R' zeros(2,1); zeros(1,2) 1]'*C*[dR' zeros(2,1); zeros(1,2) 0]; %ok
  dbdphi=-2*(dCnew(2,3)/Cnew(1,1)-Cnew(2,3)/Cnew(1,1)^2*dCnew(1,1)); %ok upto the first 2 digits, sensitive to the step size
  [izA,dizA]=intbody(z,A); %ok
  [iz0A,diz0A]=intbody(z0,A); %ok
  dsdphi=2*sign(ss)*(-1/b^2*dbdphi*(izA-iz0A)+1/b*(dizA-diz0A)*(1/2)*dbdphi);
end
  
  
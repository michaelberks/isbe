function [S,PHI,J]=cart2brst(X,Y,xA,yA,xB,yB,xC,yC,c1,c2)
%CART2BRST transform Cartesian to breast coordinates
%   [S,PHI]=CART2BRST(X,Y,XA,YA,XB,YB,XC,YC,C1,C2) transforms the Cartesian coordinate 
%   arrays X,Y to the breast coordinate arrays S and PHI. The breast
%   coordinate frame is specified by the points (XA,YA)--Nipple,
%   (XB,YB)--upper conic, (XC,YC)--lower conic
%   and the upper and lower parabola, C1 and C2, respectively. 
%
%   [S,PHI,J]=CART2BRST(X,Y,XA,YA,XB,YB,XC,YC,C1,C2) additionally returns
%   the Jacobian J of the coordinate transform. 

if xA>=xC %Test, whether the breast is on the left or the right 
  sign_phi=1;
else
  sign_phi=-1;
end

x=X(:);
y=Y(:);
sizeX=size(X);

l=length(x);

s=zeros(l,1);
phi=zeros(l,1);
if nargout>2
  J=zeros(2,2,l);  
end

%We solve the reference line lBC
lBC=cross([xB;yB;1],[xC;yC;1]); 

%Compute the direction vector at (xA,yA) to the direction of (xB,yB);
d_A=-getdirvec(c1,xA,yA,xB,yB); %We however take the opposite direction
phi_0=mod(atan2(d_A(2),d_A(1)),2*pi);

%Find the position where the common axis of the boundary parabolas intersect
%the line BC
l_axis=cross([xA;yA;1],[[0 -1; 1 0]*d_A;0]);
p_axis=cross(lBC,l_axis); p_axis=p_axis/p_axis(3);

for idx=1:l
   X=x(idx); Y=y(idx);
   if isnan(X) | isnan(Y) | ((sign([X Y 1]*c1*[X;Y;1])~=...
           sign(p_axis'*c1*p_axis)...
           & sign(l_axis'*[X;Y;1])==sign(l_axis'*[xB;yB;1]))...
           |...
           (sign([X Y 1]*c2*[X;Y;1])~=...
           sign(p_axis'*c2*p_axis)...
           & sign(l_axis'*[X;Y;1])==sign(l_axis'*[xC;yC;1])))
     phi(idx)=NaN;
     s(idx)=NaN;
     if nargout>2
        J(:,:,idx)=NaN;
     end

   else %Point is inside the boundary curves

   %if idx==5
   %    keyboard
   %end
   %Finding the
   
   %d0=distbc(0,X,Y,xA,yA,xB,yB,xC,yC,phi_0,d_A,c1,c2,l_axis,p_axis,sign_phi);
   %dpi=distbc(pi,X,Y,xA,yA,xB,yB,xC,yC,phi_0,d_A,c1,c2,l_axis,p_axis,sign_phi);   if flag==-6
   %if sign(d0)==dpi
   %    disp('Warning, fzero did not detect sign change of the function')
   %    [tmp,midx]=min([d0 dpi]);
   %    if midx==1
   %        phi(idx)=0;
    %   else
    %       phi(idx)=pi;
    %   end  
   %else
   if sign(distbc(0,X,Y,xA,yA,xB,yB,xC,yC,phi_0,d_A,c1,c2,l_axis,p_axis,sign_phi))...
   ~=sign(distbc(pi,X,Y,xA,yA,xB,yB,xC,yC,phi_0,d_A,c1,c2,l_axis,p_axis,sign_phi))
     [phi(idx),fval,eflag]=fzero(@(PHI) distbc(PHI,X,Y,xA,yA,xB,yB,xC,yC,phi_0,d_A,c1,c2,l_axis,p_axis,sign_phi),[0,pi],optimset('FunValCheck','on'));
   else %The solution is on the end point
     if abs(distbc(0,X,Y,xA,yA,xB,yB,xC,yC,phi_0,d_A,c1,c2,l_axis,p_axis,sign_phi))...
        < abs(distbc(pi,X,Y,xA,yA,xB,yB,xC,yC,phi_0,d_A,c1,c2,l_axis,p_axis,sign_phi))
       phi(idx)=0;
     else
       phi(idx)=pi;
     end
   end
   
   %l_axis were removed form the argument list as they were not used
   %end

   %end
   
   %Then we only need to convert phi to s
   if nargout<3
     cAE=phi2c(phi(idx),xA,yA,xB,yB,xC,yC,phi_0,sign_phi);
     z=x2z(X,Y,cAE);
     
     %To compute the relative distance we need the distance from the
     %pectoral muscle
     xphi=xC+phi(idx)/pi*(xB-xC); yphi=yC+phi(idx)/pi*(yB-yC);
     zphi=x2z(xphi,yphi,cAE);
     sphi=z2s(zphi,0,cAE);
     
     s(idx)=z2s(z,0,cAE)/sphi;
     
   else %Compute the Jacobian too
     [cAE,dcAE,R,dR]=phi2c(phi(idx),xA,yA,xB,yB,xC,yC,phi_0,sign_phi); % dcAE verified to be OK
     [z,dz]=x2z(X,Y,cAE,xA,yA,R,dR); %dz should be ok - the sign of dz/dphi was coorected with the hack by giving R and dR as input

     xphi=xC+phi(idx)/pi*(xB-xC); yphi=yC+phi(idx)/pi*(yB-yC);
     [zphi,dzphi]=x2z(xphi,yphi,cAE,xA,yA,R,dR);
     [sphi,D1sphi,D2sphi]=z2s(zphi,0,cAE,dcAE,R,dR);
     
     [sabs,dsdz,dsdphi]=z2s(z,0,cAE,dcAE,R,dR); %should be ok - the hack is here too to take the correct the sign of the derivative of the rotation matrix   
     s(idx)=sabs/sphi;
     
     f1=[X Y 1]*cAE*[X;Y;1];
     D1f1=2*[X Y 1]*cAE(1:2,:)';
     D2f1=[0 [X Y 1]*dcAE*[X;Y;1]];
     D1f2=[dsdz*dz(1:2)];
     D2f2=[-1 dsdz*dz(3)+dsdphi];
     D1f=[D1f1; D1f2];
     D2f=[D2f1; D2f2];
     Df=[D1f D2f];
     Jabs=-inv(D2f)*D1f;
     
     %Normalisation according to the relative distance
     dsphidphi=D1sphi*(dzphi(1)*1/pi*(xB-xC)+...
                       dzphi(2)*1/pi*(yB-yC)+...
                       dzphi(3)) + D2sphi;
     Jn=[1/sphi -sabs*sphi^(-2)*dsphidphi; 0 1];
     J(:,:,idx)=Jn*Jabs; %/sphi;
   end
   end
end

S=reshape(s,sizeX);
PHI=reshape(phi,sizeX);


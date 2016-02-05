function [X,Y,J]=brst2cart(S,PHI,xA,yA,xB,yB,xC,yC,c1,c2)
%BRST2CART transfors the breast coordinates into the Cartesian coordinates
%   [X,Y]=BRST2CART(S,PHI,XA,YA,XB,YB,XC,YC,C1,C2) transforms the 
%   breast coordinates (PHI,S)into the cartesian coordinate frame. 
%   The breast coordinate frame is specified by the points (XA,YA), (XB,YB), (XC,YC)
%   and the upper and lower parabola, C1 and C2, respectively.
%
%   [S,PHI,J]=BRST2CART(S,PHI,XA,YA,XB,YB,XC,YC,C1,C2) additionally returns
%   the Jacobian J of the coordinate transform. 

if xA>=xC %Test, whether the breast is on the left or the right 
  sign_phi=1;
else
  sign_phi=-1;
end

sizePHI=size(PHI);

PHI=PHI(:);
S=S(:);

X=zeros(size(PHI));
Y=zeros(size(PHI));
l=length(PHI);
if nargout>2
  J=zeros(2,2,l);  
end

%Compute the direction vector at (xA,yA) to the direction of (xB,yB);
d_A=-getdirvec(c1,xA,yA,xB,yB); %We however take the opposite direction
phi_0=mod(atan2(d_A(2),d_A(1)),2*pi);

for idx=1:l
   phi=PHI(idx); s=S(idx);
   if isnan(phi) | isnan(s) | phi<0 | phi>pi | s<0 %Test if outside the range
       X(idx)=NaN;
       Y(idx)=NaN;
       if nargout>2
         J(:,:,idx)=NaN;
       end
   else % Inside the range
       %Solve the parabola corresponding to phi
       if nargout<3
         cAE=phi2c(phi,xA,yA,xB,yB,xC,yC,phi_0,sign_phi);
       else %Compute the Jacobian too
         [cAE,dcAE,R,dR]=phi2c(phi,xA,yA,xB,yB,xC,yC,phi_0,sign_phi);  
       end
       
       xphi=xC+phi/pi*(xB-xC);
       yphi=yC+phi/pi*(yB-yC);
       if nargout<3
         zphi=x2z(xphi,yphi,cAE); %,xA,yA);
         sphi=z2s(zphi,0,cAE);
       else
         [zphi,dzphi]=x2z(xphi,yphi,cAE,xA,yA,R,dR);  
         [sphi,D1sphi,D2sphi]=z2s(zphi,0,cAE,dcAE,R,dR);
       end
       
       if rank(cAE)==1 %The degenerate parabola (line)
         d=[xphi-xA;yphi-yA]; d=d/norm(d);
         X(idx)=xA+d(1)*s;
         Y(idx)=yA+d(2)*s;
       else           
         z=fzero(@(Z) implicitfunction(Z,s*sphi,cAE),0);
         
         if nargout<3
           [x1,y1]=z2x(z,cAE); %Two solutions
           [x2,y2]=z2x(-z,cAE); 
         else
           [x1,y1]=z2x(z,cAE,R); %otherwise the z-result is not comparable 
           [x2,y2]=z2x(-z,cAE,R);  
         end
         %We pick the one that is on the same branch as the pectoral muscle
         %reference point

         if sign(zphi)==sign(z)
           X(idx)=x1;
           Y(idx)=y1;           
         else
           X(idx)=x2;
           Y(idx)=y2;
         end
       end
       if nargout>2 %The Jacobian computations
         [z,dz]=x2z(X(idx),Y(idx),cAE,xA,yA,R,dR);
         [tmp,dsdz,dsdphi]=z2s(z,0,cAE,dcAE,R,dR);
         
         f1=[X(idx) Y(idx) 1]*cAE*[X(idx);Y(idx);1];
         D1f1=2*[X(idx) Y(idx) 1]*cAE(1:2,:)';
         D2f1=[0 [X(idx) Y(idx) 1]*dcAE*[X(idx);Y(idx);1]];
         D1f2=[dsdz*dz(1:2)];
         D2f2=[-1 dsdz*dz(3)+dsdphi];
         D1f=[D1f1; D1f2];
         D2f=[D2f1; D2f2];
         Df=[D1f D2f];
         Jabs=-inv(D1f)*D2f;

         dsphidphi=D1sphi*(dzphi(1)*1/pi*(xB-xC)+...
                       dzphi(2)*1/pi*(yB-yC)+...
                       dzphi(3)) + D2sphi;
         sabs=sphi*s;
         Jn=inv([1/sphi -sabs*sphi^(-2)*dsphidphi; 0 1]);
         J(:,:,idx)=Jabs*Jn;
       end
   end
end

X=reshape(X,sizePHI);
Y=reshape(Y,sizePHI);

function result=implicitfunction(z,s,C)

if s==0
    result=0;
else
    result=z2s(z,0,C)-s;
end
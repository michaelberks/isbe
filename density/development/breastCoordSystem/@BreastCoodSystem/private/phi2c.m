function [C,dC,R,dR]=phi2c(phi,xA,yA,xB,yB,xC,yC,phi_0,sign_phi)

xphi=xC+phi/pi*(xB-xC);
yphi=yC+phi/pi*(yB-yC);

R=[cos(sign_phi*phi+phi_0) -sin(sign_phi*phi+phi_0); ...
   sin(sign_phi*phi+phi_0)  cos(sign_phi*phi+phi_0)]';
D=[R -R*[xA;yA];...
   0 0 1];

Cp=[(R(2,:)*[xphi-xA;yphi-yA])/(R(1,:)*[xphi-xA;yphi-yA])^2 0 0;...
    0 0 -0.5;
    0 -0.5 0];

C=D'*Cp*D;

if nargout>1
  dR=[-sign_phi*sin(sign_phi*phi+phi_0) -sign_phi*cos(sign_phi*phi+phi_0); ...
      sign_phi*cos(sign_phi*phi+phi_0)  -sign_phi*sin(sign_phi*phi+phi_0)]';  
  dD=[dR -dR*[xA;yA];0 0 0];
  dCp=[...%((dR(2,:)*[xphi-xA;yphi-yA])*(R(1,:)*[xphi-xA;yphi-yA])^2-(R(2,:)*[xphi-xA;yphi-yA])*2*(R(1,:)*[xphi-xA;yphi-yA])*(dR(1,:)*[xphi-xA;yphi-yA]))/((R(1,:)*[xphi-xA;yphi-yA])^4) 0 0;
    (dR(2,:)*[xphi-xA;yphi-yA]+R(2,:)*[1/pi*(xB-xC);1/pi*(yB-yC)])/(R(1,:)*[xphi-xA;yphi-yA])^2-2*(R(2,:)*[xphi-xA;yphi-yA])*(R(1,:)*[xphi-xA;yphi-yA])^(-3)*(dR(1,:)*[xphi-xA;yphi-yA]+R(1,:)*[1/pi*(xB-xC);1/pi*(yB-yC)]) 0 0;...
    0 0 0;
    0 0 0];
  dC=dD'*Cp*D+D'*dCp*D+D'*Cp*dD;    
end
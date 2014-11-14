function [E_S] = CalculateEquivalentPhaseFieldEnergy(gx, gy, d)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% function CalculateEquivalentPhaseFieldEnergy
% input is the whole image
%  syntax: [EPFE] = CalculateEquivalentPhaseFieldEnergy(input)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

if nargin < 3
    d = 3;
end

% win_size = 4*d;
% pad = floor(win_size/2);
pad = ceil(2*d);
beta = 0.001;

[r c] =  size(gx);

gx_pad = padarray(gx, [pad pad]);
gy_pad = padarray(gy, [pad pad]);

%Pre-allocate energy map
E_S_phi_x = zeros(r, c);

%Go through offsets for local window;
for xo = -pad:pad
    for yo = -pad:pad
        
        %Compute psi term
        dxy = sqrt(xo^2 + yo^2) / d;
        if dxy < 2
            psi = 0.5*(2 - dxy + sin(pi*dxy)/pi);
            
            %Get shifted phi_x
            sr = pad + 1 + yo;
            er = sr + r - 1;
            sc = pad + 1 + xo;
            ec = sc + c - 1;
            
            gxx = gx_pad(sr:er, sc:ec);
            gyy = gy_pad(sr:er, sc:ec);
            
            E_S_phi_x = E_S_phi_x + (gxx.*gx + gyy.*gy)*psi; 
        end
    end
end
E_S = -beta*sum(E_S_phi_x(:));

% % % [FX FY]=gradient(input);
% % 
% % 
% %     [Psi_Term]=CalculateInteractionFunction(Z);
% % %     [Phi_TermX]=CalculateGradientDotProductFunction(FX);
% %     [Phi_Term]=CalculateGradientDotProductFunction(Z);
% % %     [Phi_TermY]=CalculateGradientDotProductFunction(FY);
% % 
% % % get the composite Phi term by summing the X and Y components(?)
% % % Phi_Term=Phi_TermX+Phi_TermY;
% % EPFE=(-000.1/2).*Psi_Term.*Phi_Term;
% % 
% % end
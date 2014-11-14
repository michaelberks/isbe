function [Xilp,Xtop]=dtcwt2ilp(Yh)
% [Xtop Xilp]=dtcwt2ilp(Yh)
% Converts the DTCWT wavelet coefficients from dtwavexfm2 to a set of ILP
% coefficients Xilp; also bounces back Xtop to complete the transform.
%
% Ryan Anderson, October 2004
% Nick Kingsbury, October 2004
%
levels=size(Yh,1)-1;
Xtop=Yh{levels};

for lev=levels-1:-1:1

    % Interpolate up Yh{4} using complex interpolator.
    Yi = zeros(size(Yh{lev}));
    
    % Set up the expected phase shifts for each subband:
    % w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.5; 
    w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; % Nominally j * pi/2, but reduced a bit due to asymmetry of subband freq responses.
    p = [1 3]/4;  % Interpolation points
    
    for k = 1:6,
    % clipping extra space, for the moment
        % temp = cpxinterp(cpxinterp(Yh{lev+1}(:,:,k), p, w(k,1)).', p, w(k,2)).';
        temp = cpxinterp2(Yh{lev+1}(:,:,k), p-0.5, w(k,:),'spline');
        Yi_size=size(Yi);
        Yi(:,:,k) = temp(1:Yi_size(1),1:Yi_size(2));
    end
    
    % Double the angle of Yi.
    Yid = Yi .* Yi ./ (abs(Yi) + 1e-6);
    
    % Generate vectors at difference angles of Yh{3} and Yid
    Xilp{lev,1} = Yh{lev} .* conj(Yid);

end

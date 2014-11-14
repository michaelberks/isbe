function [Xilp, dt_top]=dtcwt2ilp3(dual_tree)
% [Xtop Xilp]=dtcwt2ilp(dual_tree)
% Converts the DTCWT wavelet coefficients from dtwavexfm2 to a set of ILP
% coefficients Xilp; also bounces back Xtop to complete the transform.
%
% Ryan Anderson, October 2004
% Nick Kingsbury, October 2004
%
levels = size(dual_tree,1)-1;
dt_top = dual_tree{levels};

for lev = 1:levels
    for ori = [2 5]
        dual_tree{lev}(:,:,ori) = -dual_tree{lev}(:,:,ori);
    end
end

Xilp = cell(levels-1,1);

for lev=levels-1:-1:1
    
    [r_level c_level] = size(dual_tree{lev}(:,:,1));
    Xilp{lev} = zeros(r_level, c_level, 6);
    
    % Interpolate coarser level using complex interpolator.
    
    % Set up the expected phase shifts for each subband:
    % w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; 
    w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; % Nominally j * pi/2, but reduced a bit due to asymmetry of subband freq responses.
    p = [1 3]/4;  % Interpolation points
    
    for ori = 1:6,
        % clipping extra space, for the moment
        temp = cpxinterp2(dual_tree{lev+1}(:,:,ori), p-0.5, w(ori,:),'spline');        
        dt_int = temp(1:r_level, 1:c_level);
        clear temp;
        
        % Double the angle of interpolated coefficients
        dt_int2 = dt_int .* dt_int ./ (abs(dt_int) + 1e-6);
        
        % Generate vectors at difference angles of current level coefficients and
        % phase-doubled, interpolated coarser level coefficients
        angle_diff = dual_tree{lev}(:,:,ori) .* conj(dt_int2);
        
        %for sub-bands 4,5,6 derotate vector by pi/2
        if ori >= 4
            angle_diff = complex(imag(angle_diff), -real(angle_diff));
        end
        
        %Now map the 2nd and 3rd complex quadrants onto the 1st and 4th by
        %flipping coefficients with negative real part
        idx = real(angle_diff) < 0; %locations to swap
        angle_diff(idx) =...
            complex(-real(angle_diff(idx)), imag(angle_diff(idx)));
        
        %Copy into output ILP structure
        Xilp{lev}(:,:,ori) = angle_diff;
    end
    


end

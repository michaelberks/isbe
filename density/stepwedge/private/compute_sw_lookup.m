function sw_lookup = compute_sw_lookup(step_g, step_h, n_bits)
%COMPUTE_SW_LOOKUP interpolate step heights and associated grey levels
%to assign a height to every grey-level on a N bit scale (i.e. from
%0-4095 for N = 12)
%
% Arguments:
%
% 'step_g'
%   - Grey-level associated with each step in the wedge
%
% 'step_h'
%   - height of each step
%
% 'n_bits (12)'
% - Bit level of grey-scale (12 = 0:4095)
%
% Outputs: none
%
% 'sw_lookup'
%   - step height associated with each grey-level
%
% See also: STEPWEDGE
%
% Created: 09-May-2006
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5125 
% Copyright: (C) University of Manchester

%Set default bit level and workout number of grey levels
if nargin < 3
    n_bits = 12;
end
g_max = (2^n_bits) - 1;

%Preallocate output
sw_lookup = zeros(g_max+1,1);
n_steps = length(step_h);

%Precompute step differentials to act as multiplier in interpolation
step_d = diff(step_h) ./ diff(step_g);

for g = 0:g_max
    % if i_pix lies outside the range of pixel values in the wedge 
    % then set x_sw to zero or max height 
    if g > step_g(1)
        sw_lookup(g +1 ) = 0;
        
    elseif g == step_g(1)
        sw_lookup(g + 1) = step_h(1);
        
    elseif g <= step_g(n_steps)
        sw_lookup(g + 1) = step_h(n_steps);    

    else
        %Find first step with grey value higher than g
        ii = 1;
        while step_g(ii) > g
            ii = ii + 1;
        end
       
        %Interpolate heights between this and previous step
        sw_lookup(g + 1) = step_h(ii-1) + step_d(ii-1) * (g - step_g(ii-1));
     
     end
end

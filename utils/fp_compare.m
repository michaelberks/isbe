function [equality] = fp_compare(a, b, range)
%
% FP_COMPARE Compare two floating point numbers
%
% % E = FP_COMPARE(A, B) If A >= B-(100*EPS) and A <= B+(100*EPS)
% then E is true, otherwise it is false, where EPS is the value returned
% by the EPS function, which is the distance from 1.0 to the next largest
% floating point number for this architecture. 
%
% E = FP_COMPARE(A, B, R) If A >= B-(0.5 * R) or A <= B+(0.5 * R)
% then E is true, otherwise it is false. R should be reasonably small,
% say 10*EPS, to ensure that the comparison is reasonable.
% If R > 0.000002, then an error is reported.
%
% See Also:
%
% EPS

if nargin == 2
	range = 200 * eps;
elseif nargin == 3 && range > 0.000002
	error('The range specified is too large')
end

equality = (a >= b-(0.5 * range) & a <= b+(0.5 * range));

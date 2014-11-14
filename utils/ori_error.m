function [orierr,err_stats] = ori_error(ori1, ori2, method, angleformat)
%ORI_ERROR Compute error between two orientations
%   [orierr,err_stats] = ori_error(ori1, ori2)
%
% Inputs:
%      ori1 - First set of orientations (e.g. ground truth)
%
%      ori2 - Second set of orientations (e.g. estimated)
%
%
% Outputs:
%      orierr - Absolute error, point for point
%
%      err_stats - summary statistics of the errors *in degrees*
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-Feb-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

if nargin<3, method = 'correct'; end
if nargin<4, angleformat = 'rad'; end

switch method
	case 'correct',
		% orientations should be in the form cos(2T) + i.sin(2T)
		% the error is then abs(angle(ori2.*conj(ori1))/2)
		
		% if in angle form then convert to complex form: cos(2T) + i.sin(2T)
		if isreal(ori1)

			% convert to radians if asked for
			if any(strcmp(angleformat,{'deg','degrees'}))
				ori1 = complex(cosd(2*ori1),sind(2*ori1));
			else
				ori1 = complex(cos(2*ori1),sin(2*ori1));
			end
		end
			
		% if in angle form then convert to complex form: cos(2T) + i.sin(2T)
		if isreal(ori2)

			% convert to radians if asked for
			if any(strcmp(angleformat,{'deg','degrees'}))
				ori2 = complex(cosd(2*ori2),sind(2*ori2));
			else
				ori2 = complex(cos(2*ori2),sin(2*ori2));
			end
		end
		
		% compute error (in radians)
		orierr = angle(ori2.*conj(ori1))/2;
	
	otherwise,
		% orientations are treated as angles, we compute the difference as it
		% is and convert to the range [-90..90] degrees
		
		% convert to angles in degrees
		if ~all(isreal(ori1)), ori1 = angle(ori1)*(180/pi); end
		if ~all(isreal(ori2)), ori2 = angle(ori2)*(180/pi); end

		% compute error (in radians)
		orierr = mb_mod(ori1-ori2,180) * pi/180;
end

if nargout>1
    
    err_stats = ori_error_stats(orierr);
end
function new_im = non_maximal_supp(im, orientations_in, orientation_type)
% NON_MAXIMAL_SUPP - Suppress non-maximal values
%   @im = non_maximal_supp(im, orientations, orientation_type)@
%
%   Suppress (set to zero) values in the 2D image @im@ that are not
%   greater than the two values on each side.  The points to compare
%   are the two that are on the normal to the @orientations@ at each
%   point.  The optional parameter @orientation_type@ says whether the
%   @orientations@ are given in @'degrees'@ (the default) or
%   @'radians'@.
  
  if nargin < 3
    orientation_type = 'degrees';
  end

  orientations = orientations_in;
    
  switch orientation_type
   case 'degrees'
    %angles = (45/2):45:180+(45/2);
    angles = 45:45:180;
    half_revolution = 180;
   case 'radians'
    %angles = pi/8:pi/4:pi+(pi/8);
    angles = pi/4:pi/4:pi;
    half_revolution = pi;
   otherwise
    error(['Unknown orientation type ''' orientation_type ''''])
  end

  % Bring all the orientations within %$[0, pi]$% (or the equivalent
  % for other orientation types).
  min_o = min(orientations(:));
  if min_o < 0
    orientations = orientations + half_revolution * -floor(min_o/half_revolution);
  end
  
  if max(orientations(:) > half_revolution)
    orientations = mod(orientations, half_revolution);
  end

  % Prepare the normal directions.
  x_inc = size(im, 1);
  y_inc = 1;

  inc = [ -y_inc, -x_inc - y_inc, -x_inc, -x_inc + y_inc, y_inc, ...
	  x_inc + y_inc];
  % Since the normal `vectors' are added and subtracted, the sign
  % does not matter, and having them all positive makes it easier
  % to test if the point is within the matrix.
  %inc = abs(inc);
  
  new_im = zeros(size(im));
  max_i = numel(im);

  for i=1:max_i,
    if im(i)
      for j=1:length(angles),
	if orientations(i) <= angles(j)
	  % Non-interpolating version.
% $$$ 	  if ~((i + inc(j) < max_i & im(i + inc(j)) >= im(i)) ...
% $$$ 		| (i - inc(j) > 0 & im(i - inc(j)) >= im(i)))
% $$$ 	    new_im(i) = im(i);
% $$$ 	  end
% $$$ 	  break
	  rad = pi/180;
	  % This should not be orientation.
	  % Interpolating version.
	  if orientations(i) < 45
	    d1 = tan(rad*orientations(i));
	    d2 = 1 - d1;
	  elseif orientations(i) < 90
	    d2 = tan(rad*(90 - orientations(i)));
	    d1 = 1 - d2;
	  elseif orientations(i) < 135
	    d1 = tan(rad*(orientations(i) - 90));
	    d2 = 1 - d1;
	  elseif orientations(i) < 180
	    d2 = 1/tan(rad*(orientations(i) - 90));
	    d1 = 1 - d2;
	  elseif orientations(i) == 180
	    d1 = 0;
	    d2 = 1;
	  else
	    error('orientation out of range.')
	  end
% $$$ 	  inc(j), inc(j+1)
% $$$ 	  d1,d2,(orientations(i)+90)
% $$$ 	  keyboard

	  if (i + inc(j)) <= max_i && (i + inc(j + 1)) <= max_i ...
		&& (i + inc(j)) > 0 && (i + inc(j + 1)) > 0
	    p1 = ((1 - d1) * im(i + inc(j)) ...
		  + (1 - d2) * im(i + inc(j + 1)));
	  else
	    p1 = -inf;
	  end
	  
	  if (i - inc(j)) <= max_i && (i - inc(j + 1)) <= max_i ...
		&& (i - inc(j)) > 0 && (i - inc(j + 1)) > 0
	    p2 = ((1 - d1) * im(i - inc(j)) ...
		  + (1 - d2) * im(i - inc(j + 1)));
	  else
	    p2 = -inf;
	  end

	  if im(i) > p1 && im(i) > p2;
	    new_im(i) = im(i);
	  end
	  
	  break
	end
      end
    end
  end
  

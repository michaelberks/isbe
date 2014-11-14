function inds_str = compute_haar_inds(feature_type,w,h,scl)

% default parameters
if nargin<2, w = 3; end
if nargin<3, h = w; end
if nargin<4, scl = 1; end

inds_str = struct(	'iimg',[],...
					'iimg45_horiz',[],...
					'iimg45_vert',[] );
				
switch lower(feature_type)
	case 'diag',
		inds_str.iimg = ...
			[	h w scl;		0 w -scl; 
				h 0 -scl;		0 0 scl;
				-1 -1 scl;		-h-1 -1 -scl;
				-1 -w-1 -scl;	-h-1 -w-1 scl;
				h -1 -scl;		0 -1 scl;
				h -w-1 scl;		0 -w-1 -scl;
				-1 w -scl;		-h-1 w scl; 
				-1 0 scl;		-h-1 0 -scl ];
			
	case 'diag45',
		inds_str.iimg45_horiz = ...
			[	0 -1 scl;		-h -h-1 -scl; 
				w -w-1 -scl;	w-h -w-h-1 scl;
				0 0 scl;		-w w -scl; 
				h h -scl;		h-w w+h scl ];
		inds_str.iimg45_vert = ...
			[	-1 0 -scl;		-h-1 -h scl; 
				-w-1 w scl;		-w-h-1 w-h -scl;
				0 0 -scl;		w -w scl; 
				h h scl;		w+h h-w -scl ];
end


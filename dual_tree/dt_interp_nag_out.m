function [full_tree] = dt_interp_nag_out(knot_mag, knot_im, knot_re, dt_dims, rows, cols, levels, whole_image)
%DT_INTERP_NAG_OUT interpolate each level of a DT-CWT to the full pixel
%grid at a given set of image locations given pre-commuted knot points (NAG
%algorithm)
%   [] = dt_interp_nag_out()
%
% Inputs:
%      dual_tree- DT-CWT cell array
%
%      rows/cols - row/column indices of points to interpolate
%
%      levels - 1D array specifying the levels of the dual-tree to include
%      in the full tree (to save memory it is useful not include the 1st
%      level for example)
%
% Outputs:
%      full_tree - (m, n, 6, num_levels)-array such that (:,:,b,l) is the
%      array of interpolated dual-tree coefficients for b-th oriented
%      sub-band in the levels(l)-th level of the original dual tree
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 11-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%If whole_image flag is on, interpret rows/cols as the dimensions of an
%image to interpolate - default is off
if nargin < 8
    whole_image = 0;
end

if nargin < 7
    levels = 1:(size(dual_tree, 1) - 1);
end

%Get number of levels to put in full tree
num_levels = length(levels);

%Get size of image/pixel subset
if whole_image
    m = rows(1);
    n = cols(1);
else
    [m n] = size(rows);
end

%Pre-allocate space for full-tree - can be an array since each level has
%same number of coeffs
full_tree = zeros(m, n, 6, num_levels);

%Set central frequencies of bands for interp method
w = [-3 -1; -sqrt(5) -sqrt(5); -1 -3; 1 -3; sqrt(5) -sqrt(5); 3 -1]*pi/2.15;

%for each sub-band, interpolate up the coefficients to the full size
for lev = 1:num_levels
    
    l = levels(lev);
    
    %Set up interpolated coordinates for full grid
    st = .5 - 2^-(l + 1);
    
    if whole_image
        xi = repmat(st + (2^-l)*(1:n), m, 1);
        yi = repmat(st + (2^-l)*(1:m)', 1, n);
    else
        xi = st + (2^-l)*cols;
        yi = st + (2^-l)*rows;
    end
    
    dy = dt_dims(l,1);
    dx = dt_dims(l,2);
    
    %Set up additional arguments needed for NAG interpolation
    mu =    [0 0 0 0 2:dy-1 dy+1 dy+1 dy+1 dy+1]';
    lamda = [0 0 0 0 2:dx-1 dx+1 dx+1 dx+1 dx+1]';
    
    %For each band interpolate DT coefficients
    for band = 1:6
        
        %make central frequencies complex
        jw = sqrt(-1) * w(band,:);
        
        if whole_image %use NAG algorithm e02df
            %Interpolate magnitudes
            z_mag_int = reshape(e02df(xi(1,:), yi(:,1), lamda, mu, knot_mag{l}(:,band)), size(xi));

            %Interpolate phases
            z_phase_unwrap = reshape(atan2(...
                e02df(xi(1,:), yi(:,1), lamda, mu, knot_im{l}(:,band)),... interpolated imaginary coeffs
                e02df(xi(1,:), yi(:,1), lamda, mu, knot_re{l}(:,band))),... interpolated real coeffs
                size(xi));
        else %use NAG algorithm e02de
            %Interpolate magnitudes
            z_mag_int = reshape(e02de(xi(:), yi(:), lamda, mu, knot_mag{l}(:,band)), size(xi));

            %Interpolate phases
            z_phase_unwrap = reshape(atan2(...
                e02de(xi(:), yi(:), lamda, mu, knot_im{l}(:,band)),... interpolated imaginary coeffs
                e02de(xi(:), yi(:), lamda, mu, knot_re{l}(:,band))),... interpolated real coeffs
                size(xi));
        end

        % Rewrap the phases and save in output structure        
        full_tree(:,:,band,lev) = z_mag_int .* exp(xi*jw(2)+yi*jw(1) + i*z_phase_unwrap);
    end
end

%finally, apply the phase correction to each of the sub-bands - multiply
%bands by - note this still needs testing (it affects the ILP calculations)
% NK suggests {1,-1,1,j,-j,j}, but I think it should be {1,-1,1,-j,j,-j} -
% after checking I agree with NK!!!!
full_tree(:,:,[1 3],:) =  i*full_tree(:,:,[1 3],:);
full_tree(:,:,[4 6],:) =   -full_tree(:,:,[4 6],:);
full_tree(:,:,2,:)     = -i*full_tree(:,:,2,:);

full_tree(:,:,[1 3 4 6],1) =   -full_tree(:,:,[1 3 4 6],1);

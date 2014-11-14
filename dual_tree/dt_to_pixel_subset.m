function [full_tree] = dt_to_pixel_subset(dual_tree, rows, cols, varargin)

%DT_TO_FULL_IMAGE interpolate each level of a DT-CWT to the full pixel grid
%   [full_tree] = dt_to_full_image(dual_tree)
%
% Inputs:
%      dual_tree- DT-CWT cell array
%
%      args.levels - 1D array specifying the args.levels of the dual-tree to include
%      in the full tree (to save memory it is useful not include the 1st
%      level for example)
%
%      args.interpmethod - string description of interpolation method passed to
%      Matlab's interp2 function
%
% Outputs:
%      full_tree - (m, n, 6, num_levels)-array such that (:,:,b,l) is the
%      array of interpolated dual-tree coefficients for b-th oriented
%      sub-band in the args.levels(l)-th level of the original dual tree
%
%
% Example:
%
% Notes:
%
% See also: COMPLEX_INTERP_FULL_IMAGE, CPXINTERP2
%
% Created: 05-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

args = u_packargs(varargin, '0',...
    'levels', [],...
    'interpmethod', 'cubic',...
    'band_fequencies', 1i*[-3 -1; -sqrt(5) -sqrt(5); -1 -3; 1 -3; sqrt(5) -sqrt(5); 3 -1]*pi/2.15,...
    'unwrap_phase', 1,...
    'interp_mag_phase', 1,...
    'correct_phase', 1);

if isempty(args.levels)
    args.levels = 1:(size(dual_tree, 1) - 1);
end

%Get number of args.levels to put in full tree
num_levels = length(args.levels);

%Get size of pixel subset
[m n] = size(rows);

%Pre-allocate space for full-tree - can be an array since each level has
%same number of coeffs
full_tree = zeros(m, n, 6, num_levels);

%Get unique set of pixels from input
[unique_pix, dummy, original_idx] = unique([rows(:) cols(:)], 'rows'); clear dummy;
original_idx = reshape(original_idx, m, n);

%for each sub-band, interpolate up the coefficients to the full size
for i_level = 1:num_levels
    for i_band = 1:6        
        if args.unwrap_phase
            band_frequency = args.band_fequencies(i_band,:);
        else
            band_frequency = [];
        end
        unique_coeffs =...
            complex_interp_pixel(dual_tree{args.levels(i_level)}(:,:,i_band),...
                unique_pix(:,1), unique_pix(:,2), args.levels(i_level), ...
                'band_frequency', band_frequency,...
                'interpmethod', args.interpmethod,...
                'interp_mag_phase', args.interp_mag_phase);
            
        full_tree(:,:,i_band,i_level) = unique_coeffs(original_idx);
        
%         full_tree(:,:,i_band,i_level) =...
%             complex_interp_pixel(dual_tree{args.levels(i_level)}(:,:,i_band),...
%                 rows, cols, args.levels(i_level), ...
%                 'band_frequency', band_frequency,...
%                 'interpmethod', args.interpmethod,...
%                 'interp_mag_phase', args.interp_mag_phase);
    end
end

%finally, apply the phase correction to each of the sub-bands - multiply
%bands by - note this still needs testing (it affects the ILP calculations)
% NK suggests {j,-j,j,-1,1,-1}, but I think it should be {1,-1,1,-j,j,-j} -
% after checking I agree with NK!!!!
if args.correct_phase
    full_tree(:,:,[1 3],:) =  1i*full_tree(:,:,[1 3],:);
    full_tree(:,:,2,:)     = -1i*full_tree(:,:,2,:);
    full_tree(:,:,[4 6],:) =   -full_tree(:,:,[4 6],:);
end


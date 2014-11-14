function [modified_region target_centre target_dt] =...
    mb_transfer_mass_dual_tree(varargin)
%MB_TRANSFER_MASS_DUAL_TREE *Insert a one line summary here*
%   [] = mb_transfer_mass_dual_tree(varargin)
%
% MB_TRANSFER_MASS_DUAL_TREE uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 05-Mar-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ... % non strict mode
		  {...  % The mandatory arguments
          'Mass'...
          },... % The optional arguments
          'TargetImage', [],...
          'TargetDualTree', [],...
          'Location', [0 0], ...
          'CoarsestLevel', 5, ...
          'MassDilate', 0, ...
          'Rotate', 0,...
          'Plot', 0 ...
          );

%Get the target dual-tree - either supplied in the arguments or generated
%from the supplied target image
if ~isempty(args.TargetImage)
    % convert TargetImage to doubles, maintain the pixel value ranges (0-255)
    args.TargetImage = double(args.TargetImage); % TargetImage is now a double
    
    % Calculate the dual_tree for the image to be synthesised
    [target_dt] = dtwavexfm2(args.TargetImage, args.CoarsestLevel+1,'near_sym_b','qshift_b');

elseif ~isempty(args.TargetDualTree)
    target_dt = args.TargetDualTree;
    args = rmfield(args, 'TargetDualTree');
else
    error('Either a target image or target dual_tree must be supplied');
end

% Generate the mass dual-tree from the mass structure;

% Get mask of mass - we should expand this cover a larger region I think...
mass_bw = poly2mask(args.Mass.mass_outline(:,1),...
    args.Mass.mass_outline(:,2),...
    size(args.Mass.mass_ROI,1), size(args.Mass.mass_ROI,2)); 

mass_centre = mean(args.Mass.mass_outline);
%If we're rotating, do this now
if args.Rotate
    
    %Pad args.Mass.background_ROI
    [r c] = size(args.Mass.background_ROI);
    pad_vec = ceil(0.5*(sqrt(r^2 + c^2) - [r c]));
    mass_image_pad = padarray(args.Mass.background_ROI, pad_vec, 'symmetric');
    
    %Rotate
    mass_image_pad_rot = imrotate(mass_image_pad, 180*(args.Rotate/pi), 'bilinear', 'crop');
    
    %Take central region
    args.Mass.background_ROI = mass_image_pad_rot(pad_vec(1)+(1:r), pad_vec(2)+(1:c));
    
    %No need to pad mass_bw as zeros on edges anyway
    mass_bw = imrotate(mass_bw, 180*(args.Rotate/pi), 'nearest', 'crop');
    
    %correct the mass centroid
    r_theta = [cos(args.Rotate) sin(args.Rotate); -sin(args.Rotate) cos(args.Rotate)];
    mass_centre = (r_theta*(mass_centre - [c r]/2)')' + [c r]/2;
end

%Take the DT-CWT of the mass image
mass_dt = dtwavexfm2(args.Mass.background_ROI,...
    args.CoarsestLevel+1,'near_sym_b','qshift_b');

% Get row/columns of the mass region
%mass_bw = imresize(mass_bw, 0.5);
mass_bw = imdilate(mass_bw, strel('disk', args.MassDilate));

[mass_r mass_c] = find(mass_bw);

%effectively the upper-left point of the bounding rectangle that will house
%the mass region we are inserting - it may make sense to add
%the mass centroid to define the centre of the mass in the target
%region and return this
% make the offset a multiple of 32
target_offset = args.Location - mass_centre;
target_offset = target_offset - rem(target_offset, 32);
target_centre = target_offset + mass_centre;

%Now transfer the dual_tree data across from the mass to the target in each
%level
for level = 1:args.CoarsestLevel
    
    %Downsample the row/column coefficients of the region w.r.t the mass
    mass_r_level = ceil(mass_r ./ 2^level);
    mass_c_level = ceil(mass_c ./ 2^level);
    
    %Downsample target offset (shouldn't need to ceiling this as it should
    %be divisible by 16)
    level_offset = target_offset ./ 2^level;
    
    %Add the offset to get row/column w.r.t to the target
    target_r_level = mass_r_level + level_offset(2);
    target_c_level = mass_c_level + level_offset(1);
    
    if args.Plot
        figure; plot(target_c_level, target_r_level, 'rx'); axis image;
    end
    
    %Convert row/column subs to indices
    mass_idx = sub2ind(size(mass_dt{level}(:,:,1)), mass_r_level, mass_c_level);
    target_idx = sub2ind(size(target_dt{level}(:,:,1)), target_r_level, target_c_level);
    
    %Now transfer across - because of the way this is indexed we have to
    %this subband at a time using temprorary arrays to hold each sub-band
    for band = 1:6
        mass_band = mass_dt{level}(:,:,band);
        target_band = target_dt{level}(:,:,band);
        target_band(target_idx) = mass_band(mass_idx);
        target_dt{level}(:,:,band) = target_band;
    end
end

[modified_region] = dtwaveifm2(target_dt, 1, 'near_sym_b', 'qshift_b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%My plan goes a little something like this...

% 1) Import 4th (5th?) level of scaling coefficients from mass region into
% target region

% 2) Filter by the dual_tree method to obtain new 4th dual_tree
% coefficients

% 3) Now substitute the ILP/ICP transform coefficients from 3rd level
% upwards

% 4) Reconstruct the full dual-tree and then reconstruct the region

% That's bullshit. The whole point is we want the coarse level maintained
% not swapped... think again.
function [frames, cell_positions] = make_synthetic_flow_frames(vessel_name, mean_flow, varargin)
%MAKE_SYNTHETIC_FLOW_COMPARISON_VIDEO generates synthetic flow videos using real
%capillary shapes
%   [] = make_synthetic_flow_comparison_video('normalapex0008', 3, 1.5, 'adjacent_videos', 0)
%
% MAKE_SYNTHETIC_FLOW_COMPARISON_VIDEO uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names.
%
% Mandatory Arguments:
% - vessel_name: name of the real vessel contour used as the shape (defines
% the vessel's centre line, inner and outer edges and locations of any
% apices in variables: vessel_centre, inner_edge, outer_edge (all n x 2)
% and apex_idx (m x 1, where vessel_centre(apex_idx(i),:) are the x,y
% coordinates of the i'th apex). For normal 'single loop' capillaries m = 1
%
% - mean_flow: the maximum flow (effectively the number of pixels a red
% blood moves in a signle frame) in the synthesised flow field. Flow in the
% rest of the vessel will be scaled to the inverse square of vessel width
% so that a constant theoretical volumetric flow rate is maintained
% throught the vessel. One set of coloured tags in the final will move at the
% exact rate defined by this flow field
%
%
% See also:
%
% Created: 15-Feb-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'contour_dir',          'C:\isbe\nailfold\data\rsa_study\set12g\vessel_contours\',... %Where the contour data files live
    'overwrite',            true,... %If true, saves the new flow field over the top of any existing data for the contour in the video directory
    'flow_map_dir',         [],...
    'cell_density',         8,... %Density of cells in the synthesised video - try adjusting to see how it affects the ability to perceive flow
    'num_frames',           120,... %Together with frame_rate defines video length
    'vessel_contrast',      0.5,...
    'speckle_noise',        0.01,... %Changes how noisy the videos are, increasing speckle noise decreases video quality, may be worth experimenting with to see if effects how easy the task is
    'max_jitter',           [4 2],...
    'scale_cell_width',     0,... %Can choose to scale red blood cell size to the vessel width, physiologically this doesn't make much sense, so best leave at 0
    'fixed_cell_width',     10,... %If above is 0, this defines blood cell width
    'debug',                false,...
    'plot', 0);
clear varargin;

if isempty(args.flow_map_dir)

    %Load in the vessel contour - will load variables vessel_centre,
    %inner_edge, outer_edge and apex_idx
    load([args.contour_dir vessel_name '_vc.mat'],...
        'vessel_centre', 'inner_edge', 'outer_edge', 'apex_idx');
    apex_idx = sort(apex_idx); %#ok
    
    %Make a flowmap for the vessel shape
    [flowmap, mask, vessel_centre, widths] = ...
        create_flowmap_profile(vessel_centre, inner_edge, outer_edge, apex_idx, []); %#ok
else
    %Load pre-computed flow-map
    load([args.flow_map_dir vessel_name '_fm.mat'],...
        'flowmap', 'mask', 'vessel_centre', 'widths');
end

%Flow mmap is scaled so that its max flow = 1. Scale instead so the mean
%flow = 1, then scale by the mean_flow input parameter
flat_map = max(flowmap, [], 3);
flat_mask = any(mask,3);
flowmap = mean_flow*flowmap / nanmean(abs(flat_map(flat_mask)));

%Set the cell size and compute number of cells to use
if args.scale_cell_width
    cell_sz = min(widths(:));
else
    cell_sz = args.fixed_cell_width;    
end
%Doing this by area only seems to work well for normal size vessels,
%because their volume is commensurate with their width, for giants we
%need to consider cell and vessel volume
%n_pixels_per_cell = cell_sz^2;
%n_mask_pixels = sum(mask(:)>0);
cell_vol = cell_sz^3;
vessel_vol = sum(mask(:)>0) * mean(widths);
num_cells = ceil(args.cell_density * vessel_vol / cell_vol);

%Generate cell positions for each frame
cell_positions = generate_cell_positions_layers(flowmap, mask, ...
    vessel_centre, widths, ...
    num_cells, args.num_frames+1);

%Get size of flow map
[flow_h, flow_w, ~] = size(flowmap);

%Make background texture
cloud_bg = noiseonf(max(flow_h, flow_w), 3);
background = normim(cloud_bg(1:flow_h,1:flow_w), 'stretch_fixed');

if all(args.max_jitter)
    frame_jitter = rand(args.num_frames, 2);
    g_j = gaussian_filters_1d(4);
    g_j = g_j(:) / sum(g_j);
    frame_jitter = conv2(frame_jitter, g_j, 'same');
end

%--------------------------------------------------------------------------
%Now we can make our video
%--------------------------------------------------------------------------
frames = zeros(flow_h, flow_w, args.num_frames);
for i_fr = 1:args.num_frames;

    %Make the synthetic frame
    cloud_vessel = noiseonf(max(flow_h, flow_w), 1.0);
    contrast = args.vessel_contrast *(0.5 + 0.5*normim(cloud_vessel(1:flow_h,1:flow_w), 'stretch_fixed'));

    cell_pos_i = cell_positions(:,:,i_fr:i_fr+1);
    if all(args.max_jitter)
        cell_pos_i(:,1) = cell_pos_i(:,1) + frame_jitter(i_fr,1)*args.max_jitter(1);
        cell_pos_i(:,2) = cell_pos_i(:,2) + frame_jitter(i_fr,2)*args.max_jitter(2);
    end
    frame_i = make_frame(cell_pos_i, ...
                     [flow_h flow_w], cell_sz, ...
                     [], background, contrast, inf, 0); %
    frame_i = frame_i - min(frame_i(:));
    frame_i = frame_i / max(frame_i(:));
    if args.speckle_noise
        frames(:,:,i_fr) = imnoise(frame_i, 'gaussian', 0, args.speckle_noise);
    else
        frames(:,:,i_fr) = frame_i;
    end

end
g_lims = [min(frames(:)) max(frames(:))];
g_range = g_lims(2) - g_lims(1);
frames = (frames-g_lims(1))/g_range;
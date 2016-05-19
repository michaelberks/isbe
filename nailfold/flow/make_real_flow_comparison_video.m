function [] = make_real_flow_comparison_video(vessel_name, speed_factor, varargin)
%MAKE_REAL_FLOW_COMPARISON_VIDEO adds cell markers to real flow videos to
%compare observed flow rate
%capillary shapes
%   [] = make_real_flow_comparison_video('001wellcome_2015_02_27_L1_09_38_22_s10_v01', 1.5, 'adjacent_videos', 1)
%
% MAKE_REAL_FLOW_COMPARISON_VIDEO uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names.
%
% Mandatory Arguments:
% - vessel_name: name of the real vessel
%
%
% - speed_factor: the factor by which the 'wrong' display tags are scaled.
% This scaling is applied equally throughout the flow field, so for example
% if speed_factor = 1.5, the 'wrong' tags will be 50% faster than the true
% flow throughout the vessel
%
% Optional Arguments: see unpack arguments below
%
% Outputs:
%
% Example: [] = make_real_flow_comparison_video('001wellcome_2015_02_27_L1_09_38_22_s10_v01', 1.5, 'adjacent_videos', 1)
%
% Notes: 
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
    'data_dir',             '\\nasr.man.ac.uk\epsrss$\snapped\replicated\MSI\flow_project\flow_data\',... %Where the contour data files live
    'results_dir',          '\\nasr.man.ac.uk\epsrss$\snapped\replicated\MSI\flow_project\flow_results\',... %Where the contour data files live
    'centrelines_dir',      '\\nasr.man.ac.uk\epsrss$\snapped\replicated\MSI\flow_project\flow_centrelines\',... %Where the contour data files live
    'video_dir',            '\\nasr.man.ac.uk\epsrss$\snapped\replicated\MSI\flow_project\flow_videos\real\',... %Where the videos (and associated flow data) are saved to
    'speckle_noise',        0.01,... %Changes how noisy the videos are, increasing speckle noise decreases video quality, may be worth experimenting with to see if effects how easy the task is
    'snap_to_centre',       true,... %Leave true, not snapping to centre doesn't work properly!
    'tag_spacing',          250,... %Approx number of pixels along the vessel path between starting positions of the display tags
    'scale_tag_width',      0,... %Can scale tag size to vessel width, can experiment to see which you think makes the task easier
    'fixed_tag_width',      2, ... %If above is zero, defines the radius of the tags
    'solid_tags',           true,... %If false, creates rings rather than solids circles
    'adjacent_videos',      false,... %If false produces an overlay video with one vessel and red green tags, if true creates adjacent videos with red tags.
    'temp_frames_dir',      [],... if empty, uses the system temporary dir
    'frame_rate',           30,...
    'debug',                false,...
    'plot', 0);
clear varargin;

flow_data_dir = 'O:\flow_project\flow_data\';
flow_results_dir = 'O:\flow_project\flow_results\';
flow_centrelines_dir = 'O:\flow_project\flow_centrelines\';

%We need this created (if not already) to save data to
if ~exist(args.video_dir, 'dir')
	mkdir(args.video_dir);
end

load([flow_data_dir vessel_name '.mat'], 'cropped_frames');
load([flow_results_dir vessel_name '.mat'], 'flow_results');
load([flow_centrelines_dir vessel_name '.mat'], 'vessel_centre');

%Get size of flow map
flowmap = flow_results.flowPyramidEst{1};

%Run morph dilation over the flow field by an approx radius of cell size,
%this gives a fairer reflection of the predicted flow
flowmap = exp(1i*angle(flowmap)).*imdilate(abs(flowmap), strel('disk', 4));

%Smooth the flow
g_prob = gaussian_filters_1d(2);
g_prob = g_prob / sum(g_prob);
flowmap = conv2(g_prob', g_prob, flowmap, 'same');

[flow_h flow_w] = size(flowmap);

%Prepare the frames
if flow_h > size(cropped_frames,1) %#ok
    cropped_frames(flow_h,:,:) = min(cropped_frames(:));
end
if flow_w > size(cropped_frames,2)
    cropped_frames(:,flow_w,:) = min(cropped_frames(:));
end
g_range = flow_results.g_lims(2) - flow_results.g_lims(1);
cropped_frames = 255*(cropped_frames-flow_results.g_lims(1))/g_range;    

%Set the display tag radius - fixed or based on the frame width?
if args.scale_tag_width
    tag_radius = widths(apex_idx(1)) * args.scale_tag_width * 0.25;
else
    tag_radius = args.fixed_tag_width;    
end
tag_radius2 = tag_radius^2;

%Set the first point, the select equally spaced further points on the
%vessel
first_pt = ceil(tag_radius) + ceil(tag_radius*rand);

vessel_centre_hi = spline_contour(vessel_centre, [], 0.1);
num_pts = size(vessel_centre_hi,1); 

num_tags = ceil(num_pts / args.tag_spacing);
idx = round(linspace(first_pt, num_pts, num_tags+1));
idx(end) = [];

%Tag points define x,y coordinate and flow layer
tag_pts = zeros(num_tags, 3, 2);
tag_pts(:,1:2,1) = vessel_centre_hi(idx,:);
tag_pts(:,1:2,2) = vessel_centre_hi(idx+2,:);

tag_pts(:,3,1) = idx;
tag_pts(:,3,2) = idx+2;

if args.debug
    tag_history = zeros(num_tags, 3, 2, args.num_frames);
    tag_history(:,:,1) = tag_pts;
end

%Set up points matrices to compute tag masks
xx = repmat(1:flow_w, flow_h, 1);
yy = repmat((1:flow_h)', 1, flow_w);

%Make a temporary directory to store frame in
if isempty(args.temp_frames_dir)
    temp_frames_dir = tempname;
else
    temp_frames_dir = args.temp_frames_dir;
end
mkdir(temp_frames_dir);

%Choose which colour tags are the 'true' tags
if rand > 0.5
    scaling = [1.0 speed_factor];
else
    scaling = [speed_factor 1.0];
end

%--------------------------------------------------------------------------
%Now we can make our video
%--------------------------------------------------------------------------
num_frames = size(cropped_frames,3);
for i_fr = 1:num_frames;

    if args.adjacent_videos
        vid_frame = zeros(flow_h, 2*flow_w + 16, 3, 'uint8'); 
    else
        vid_frame = zeros(flow_h, flow_w, 3, 'uint8'); 
    end
    
    %Take real frame
    vid_frame_i = cropped_frames(:,:,i_fr);
                 
    if i_fr > 1
        tag_history(:,:,:,i_fr) = tag_pts;
    end

    for i_dim = 1:2
        %Add the display tags to this frame
        frame_i = vid_frame_i;

        %For each tag
        for i_tag = 1:num_tags

            %Make cell mask at current position
            ix = tag_pts(i_tag,1,i_dim);
            iy = tag_pts(i_tag,2,i_dim);
            ii = tag_pts(i_tag,3,i_dim);
            
            xy2 = (xx-ix).^2 + (yy-iy).^2;
            cell_mask = xy2 < tag_radius2;
            if ~args.solid_tags
                cell_mask = cell_mask & xy2 > (tag_radius2/4);
            end
            frame_i(cell_mask) = 255;    

            %If we're still in a valid position, move to next position
            %based on flow estimate
            if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h;                 
                
                fxyi = flowmap(round(iy), round(ix));
                
                if ~isnan(fxyi)
                    jx = ix + scaling(i_dim)*real(fxyi);
                    jy = iy + scaling(i_dim)*imag(fxyi);
                    
                    if args.snap_to_centre
                        %Allow to move only to vessel points in this layer
                        %or the first 5 points of the next layer to avoid
                        %crossing on to the next layer (while still alowing
                        %it to transition layers)!
                        available_pts = ii+(1:12);
                        available_pts(available_pts > num_pts) = [];
                        
                        d = ...
                            (vessel_centre_hi(available_pts,1)-jx).^2 +...
                            (vessel_centre_hi(available_pts,2)-jy).^2;
                        [~,idx] = min(d);
                        
                        ji = idx + ii;
                        
                        if ji == num_pts
                            jx = vessel_centre_hi(first_pt,1);
                            jy = vessel_centre_hi(first_pt,2);
                            ji = 1;
                        else                            
                            jx = vessel_centre_hi(ji,1);
                            jy = vessel_centre_hi(ji,2);
                        end                        
                        
                    else
                        ji = ii;
                    end
                    tag_pts(i_tag,1,i_dim) = jx;
                    tag_pts(i_tag,2,i_dim) = jy;
                    tag_pts(i_tag,3,i_dim) = ji;
                else
                    %Otherwise start at the beginning
                    tag_pts(i_tag,:,i_dim) = [vessel_centre_hi(first_pt,1:2) first_pt];
                end
            else
                %Otherwise start at the beginning
                tag_pts(i_tag,:,i_dim) = [vessel_centre_hi(first_pt,1:2) 1 first_pt];
            end
        end

        %Put frames either on top or by one another
        if args.adjacent_videos
            vid_frame(:, (1:flow_w) + (i_dim-1)*(flow_w + 16), 1) = uint8(frame_i);
            vid_frame(:, (1:flow_w) + (i_dim-1)*(flow_w + 16), 2) = uint8(vid_frame_i);
            vid_frame(:, (1:flow_w) + (i_dim-1)*(flow_w + 16), 3) = uint8(vid_frame_i);
        else
            vid_frame(:, :, i_dim) = uint8(frame_i);
        end

    end
    if ~args.adjacent_videos
        vid_frame(:, :, 3) = uint8(vid_frame_i);
    end

    %Write the frame to the temporary directory
    imwrite(vid_frame,...
        [temp_frames_dir '\frame' zerostr(i_fr,4) '.bmp']);
end

%Make path to save video
if args.adjacent_videos
    video_path = [args.video_dir vessel_name ...
        '_f1_' num2str(scaling(1),3) '_f2_' num2str(scaling(2),3) '_adjacent.mp4'];
else
    video_path = [args.video_dir vessel_name ...
        '_f1_' num2str(scaling(1),3) '_f2_' num2str(scaling(2),3) '_overlay.mp4'];
end
cmd = ['ffmpeg -y -r ' num2str(args.frame_rate) ' -i "' temp_frames_dir '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' video_path '"'];
system(cmd);

delete([temp_frames_dir '\*']);
rmdir(temp_frames_dir);
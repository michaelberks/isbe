function [seq_dirs] = phantom_flow(root_dir, seq_dirs, do_make_videos, do_compute_flow)
%PHANTOM_FLOW make flow videos and compute flow estimates for phantom data
%   [] = select_vessels_from_candidates(varargin)
%
%
% Mandatory Arguments:
%
% Optional Arguments:
%   root_dir - path to folder where sequences are stored
%
%   seq_dirs - cell list of sequences to process. Can be left empty in
%   which case all sequences containing frames will be detected and
%   processed
%
%   do_make_videos - flag to choose whether or not to make videos (0 = no,
%   1 = yes, default = yes)
%
%   do_compute_flow - flag to choose whether or not to compute flow estimates
%   (0 = no, 1 = yes, default = yes)
%
% Outputs:
%
%   seq_dirs - list of sequences folders containing frames, can be saved
%   for future use to save time searching for new folders each time
%
% Example: seq_dirs = phantom_flow('O:\Summer Placement 2016\Microbeads\Polystyrene\', [], 0, 1);
%
% Notes:
%
% See also:
%
% Created: 20-Jul-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 


%% Get sequence directories where frames are if not supplied as input
if ~exist('root_dir', 'var') || isempty(root_dir)
    root_dir = 'O:\Summer Placement 2016\Microbeads\';
end
if ~exist('seq_dirs', 'var') || isempty(seq_dirs)
    seq_dirs = make_seq_dirs(root_dir, cell(0,1));
    display('Found frames in the following sequence directories:');
    display(seq_dirs);
end    
if ~exist('do_make_videos', 'var')
    do_make_videos = true;
end
if ~exist('do_compute_flow', 'var')
    do_compute_flow = true;
end 

for i_s = 1:length(seq_dirs)
    seq_dir = seq_dirs{i_s};
    frames_list = dir([seq_dir '*.bmp']);

    num_repeats = length(frames_list) / 180;
    frames = zeros(480, 160, 180);
    
    display(['Checking frames in ' seq_dirs{i_s} ', found ' num2str(num_repeats) ' repeat(s)']);
    for i_rpt = 1:num_repeats
        mkdir([seq_dir '\v' zerostr(i_rpt,2)]);
        
        if do_make_videos
            display(['Making video for ' seq_dirs{i_s} ' repeat ' num2str(i_rpt)]);
            for i_f = (i_rpt-1)*180 + (1:180)
                copyfile([seq_dir frames_list(i_f).name],...
                    [seq_dir '\v' zerostr(i_rpt,2) '\frame' zerostr(i_f - 180*(i_rpt-1),4) '.bmp']);
            end
            cmd = ['ffmpeg -y -r 60 -i "' seq_dir '\v' zerostr(i_rpt,2) ...
                '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' seq_dir 'v' zerostr(i_rpt,2) '\flow.mp4"'];
            system(cmd);
            delete([seq_dir 'v' zerostr(i_rpt,2) '\*.bmp']);
        end
        if do_compute_flow
            display(['Loading frames for ' seq_dirs{i_s} ' repeat ' num2str(i_rpt)]);
                 
            for i_f = (i_rpt-1)*180 + (1:180)
                frame = imread([seq_dir frames_list(i_f).name]);
                frames(:,:,i_f - (i_rpt-1)*180) = double(frame(:,240 + (1:160)));
            end

            %Contrast normalise the frames
            g_lims = prctile(frames(:), [1 99]);
            g_range = g_lims(2) - g_lims(1);

            display(['Computing flow for ' seq_dirs{i_s} ' repeat ' num2str(i_rpt)]);
            %Compute flow results
            flow_results = [];
            [flow_results.flowPyramidEst, flow_results.flowConfidence] = ...
                estimate_flow_multilevel(255*(frames-g_lims(1))/g_range, [],[],1:4);   
            flow_results.g_lims = g_lims;   
            save([seq_dir '\v' zerostr(i_rpt,2) '\flow_results.mat'], 'flow_results');
        end       
    end
end

%% ------------------------------------------------------------------------
function seq_dirs = make_seq_dirs(this_dir, seq_dirs)
[contains_frames subdir_list] = check_for_frames(this_dir);

if contains_frames
    seq_dirs(end+1,1) = {this_dir};
end
for i_dir = 1:length(subdir_list)
    seq_dirs = make_seq_dirs(subdir_list{i_dir}, seq_dirs);
end

%% ------------------------------------------------------------------------
function [contains_frames subdir_list] = check_for_frames(dirname)


contains_frames = false;
subdir_list = cell(0,1);
bmplist = dir([dirname '*.bmp']);

if ~isempty(bmplist)
    contains_frames = true;
    return;
end

filelist = dir(dirname);
for i_f = 3:length(filelist)
    if filelist(i_f).isdir
    	subdir_list(end+1,1) = {[dirname filelist(i_f).name '\']}; %#ok
    end
end

%% ------------------------------------------------------------------------
% root_dir = 'O:\flow_project\camera_frames\500um_videos\';
% %%
% seq_dirs = {...
%     [root_dir '\4ml\h_500_1\2016_05_03\'];
%     [root_dir '\5ml\h_500_1\2016_05_03\'];
%     [root_dir '\5.5ml\h_500_1\2016_05_03\'];
%     [root_dir '\6ml\h_500_1\2016_05_03\'];
%     [root_dir '\6.5ml\h_500_1\2016_05_03\'];
%     [root_dir '\7ml\h_500_1\2016_05_03\'];
%     [root_dir '\7.5ml\h_500_1\2016_05_03\'];
%     [root_dir '\8ml\h_500_1\2016_05_03\'];
%     [root_dir '\8.5ml\h_500_1\2016_05_03\'];
%     [root_dir '\9ml\h_500_1\2016_05_03\'];
%     [root_dir '\9.5ml\h_500_1\2016_05_03\'];
%     [root_dir '\10ml\h_500_1\2016_05_03\'];
%     [root_dir '\70ml\h_500_1\2016_05_03\']};
% 
% %%
% for i_s = 9:length(seq_dirs)
%     seq_dir = seq_dirs{i_s};
%     frames = dir([seq_dir '*.bmp']);
% 
%     num_repeats = length(frames) / 180;
%     for i_rpt = 1:num_repeats
%         mkdir([seq_dir 'v' zerostr(i_rpt,2)]);
%         for i_f = (i_rpt-1)*180 + (1:180)
%             copyfile([seq_dir frames(i_f).name],...
%                 [seq_dir 'v' zerostr(i_rpt,2) '\frame' zerostr(i_f - 180*(i_rpt-1),4) '.bmp']);
%         end
%         cmd = ['ffmpeg -y -r 60 -i "' seq_dir 'v' zerostr(i_rpt,2) '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' seq_dir 'v' zerostr(i_rpt,2) '\flow.mp4"'];
%         system(cmd);
%         delete([seq_dir 'v' zerostr(i_rpt,2) '\*.bmp']);
%     end
% end
% %%
% for i_s = 8%1:length(seq_dirs)
%     seq_dir = seq_dirs{i_s};
%     frames_list = dir([seq_dir '*.bmp']);
% 
%     num_repeats = length(frames_list) / 180;
%     for i_rpt = 2%1:num_repeats
%         frames = zeros(480, 160, 180);
%         
%         for i_f = (i_rpt-1)*180 + (1:180)
%             frame = imread([seq_dir frames_list(i_f).name]);
%             frames(:,:,i_f - (i_rpt-1)*180) = double(frame(:,240 + (1:160)));
%         end
%         
%         %Contrast normalise the frames
%         g_lims = prctile(frames(:), [1 99]);
%         g_range = g_lims(2) - g_lims(1);
% 
%         %Compute flow results
%         flow_results = [];
%         [flow_results.flowPyramidEst, flow_results.flowConfidence] = ...
%             estimate_flow_multilevel(255*(frames-g_lims(1))/g_range, [],[],1:4);   
%         flow_results.g_lims = g_lims;   
%         save([seq_dir 'v' zerostr(i_rpt,2) '\flow_results.mat'], 'flow_results');
%         
%     end
% end







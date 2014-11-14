function [image_list, args] = create_image_lists(varargin)
% Create a structure that holds the filenames for all images of interest to
% a predictor

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[image_list, args] = func(varargin);


%% The function
function [image_list, args] = func(varargin)

args = u_packargs(varargin, ...
    '0', ... % non-strict mode
    ... % the mandatory arguments
   {'image_root', ...
    'image_dir'}, ... 
    ... % What images do we want
    'image_set', {'image','fov_mask','fg_mask','ori','width'}, ...
    ... % the optional arguments
	'fov_mask_dir', [], ...
	'fg_mask_dir', [], ...
	'ori_dir', [], ...
    'width_dir', [], ...
    ... % estimated measurements (added to image list if defined)
    'probability_dir', [], ...
    'ori_est_dir', [], ...
    'width_est_dir', [],...
    'shift_images', 0);

% add estimated quantities to list of images if defined
if ~isempty(args.probability_dir),
    args.image_set{end+1} = 'probability';
end
if ~isempty(args.ori_est_dir),
    args.image_set{end+1} = 'ori_est';
end
if ~isempty(args.width_est_dir),
    args.image_set{end+1} = 'width_est';
end    
    
% remove any duplicated image flags
args.image_set = unique(args.image_set);

% Create empty image list
image_list = cell2struct(cell(length(args.image_set),1), args.image_set);

% Create structure with filenames in it
for fields = {'image', 'fg_mask', 'fov_mask', 'ori', 'width', ...
              'probability', 'ori_est', 'width_est'}
    fieldname = fields{1};

    % If fieldname is not in the desired image_set then skip to next field
    if ~any(strcmp(args.image_set, fieldname)), continue; end

    % If image folder has not been given then skip to next field
    subdir = args.([fieldname,'_dir']);
    if isempty(subdir), continue; end

    % Search the path for .mat files and add filenames to the structure
    datapath = prettypath(subdir);
    d = dir([datapath,'/*.mat']);
    if ~isempty(d)
        % deal() does something like this but I couldn't get it 
        % quite right
        for i = 1:numel(d)
            image_list(i).(fieldname) = [datapath '/' d(i).name];
        end            
    else
        display(['Warning, no ' fieldname ' files found matching the path ' datapath,'/*.mat']);
    end    
end
if args.shift_images
    image_list = circshift(image_list, [0 args.shift_images]);
end


%% Test script
function test_script()
clc;
args = default_args();
sampling_args = get_sampling_args_from(args);
image_list = create_image_lists(sampling_args);
disp(length(image_list));
image_list(1)
    
function [bg, bg_index] = get_random_background(args, varargin)
% Randomly select background - will load a matrix 'bg'

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[bg, bg_index] = func(args, varargin{:});


%% The function
function [bg, bg_index] = func(args, bg_list)

if exist('bg_list', 'var')
    % Select a random background from the given list
    bg_index = ceil(length(bg_list)*rand);
    bg_filename = bg_list(bg_index).name;
else
    % Select a random background from all possible backgrounds
    if isempty(args.bg_zeros)
        bg_zeros = floor(log10(args.num_bgs)) + 1;
    else
        bg_zeros = args.bg_zeros;
    end
    bg_index = ceil(args.num_bgs*rand);
    bg_filename = [args.bg_stem zerostr(bg_index, bg_zeros) '.' args.bg_fmt];
end

% prepend directory
bg_filename = [args.bg_dir '/' bg_filename];
bg_filename = prettypath(bg_filename);

% if background not found then skip and try another
if ~exist(bg_filename,'file')
    warning([bg_filename,' not found']);
    bg = [];
    return;
end

% load background according to specified format
[pathname,filebase,extension] = fileparts(bg_filename);
switch extension
    case '.mat', bg = u_load(bg_filename);
    case '.png', bg = double(imread(bg_filename));
end									


%% Test script
function test_script()
clc;

rootpath = 'U:\projects\mammography\data\synthetic_backgrounds\';
args.bg_dir = [rootpath,'real512\train/'];
args.bg_fmt = 'mat';
args.bg_stem = 'bg';
args.bg_zeros = 5;
args.num_bgs = 10;

[bg, bg_index] = get_random_background(args);
figure(1); clf; hold on; colormap(gray(256));
imagesc(bg); axis('image','ij');
disp(bg_index);

args.bg_dir = [rootpath,'smooth512\train/'];
bg_list = dir([args.bg_dir,'*.mat']);
[bg, bg_index] = get_random_background(args, bg_list);
figure(2); clf; hold on; colormap(gray(256));
imagesc(bg); axis('image','ij');
disp(bg_index);


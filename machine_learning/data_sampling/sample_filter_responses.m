function [samples] = sample_filter_responses(responses, rows, cols, varargin)
% SAMPLE_FILTER_RESPONSES Given a 
%   [dt_samples] = sample_g2d_data(dt, rows, cols, rotate)
%
% Inputs:
%      responses - [n_rows x n_cols x n_levels x n_filters] matrix of
%      filter responses at multiple scales
%
%      rows - row subscripts of sample locations
%
%      cols - column subscripts of sample locations
%
%      win_size - dimensions of window to sample
%
%      pca - parameters of Principal Component Analysis model
%
%      subtract_mean - normalize responses with respect to their mean value
%
%
% Outputs:
%
%      samples - the sampled feature responses, concatenated over the local
%      neighbourhood, over scales and over filters
%
% Example:
%
%
% See also:
%
% Created: 23-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[samples] = func(responses, rows, cols, varargin{:});


%% The function
function [samples] = func(responses, rows, cols, varargin)

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    'do_max', 0,...
    'rotate',0,...
    'win_size', 3,...
    'subtract_mean', false);

win_size = args.win_size; 
[n_rows, n_cols, n_levels, n_filters] = size(responses);
n_pixels = win_size * win_size;


%First check if we need to rotate the filters or take the maximum response
if args.do_max || args.rotate
    [responses n_filters] = reorder_bands(responses, args.do_max);
end
    
% rows and cols must be column vectors
rows = rows(:);
cols = cols(:);

% Make copies of sample rows and cols at positions of local window patch
[rr, cc, num_samples] = get_indices_from_winsize(rows, cols, win_size);

% Note: patch_idx is a matrix
patch_idx = sub2ind([n_rows, n_cols], rr, cc);

% allocate space
samples = zeros(num_samples, n_filters * n_levels * n_pixels);

cols = 1:n_pixels;
for i_filter = 1:n_filters
    for i_level = 1:n_levels
        filter_response = responses(:,:,i_level,i_filter);
		samples(:, cols) = filter_response(patch_idx);
        cols = cols + n_pixels;
    end
end

% Subtract mean if requested
if args.subtract_mean
    samples = bsxfun(@minus, samples, mean(samples, 2));
end


%% Old version for comparison
function [samples] = old_func(responses, rows, cols, varargin)

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    'win_size', 3,...
    'pca', [], ...
	'sigma_range', []); % not used but prevents warnings

win_size = args.win_size; 
[n_rows, n_cols, n_levels, n_filters] = size(responses);
n_pixels = win_size * win_size;

% rows and cols must be column vectors
rows = rows(:);
cols = cols(:);

% Make copies of sample rows and cols at positions of local window patch
[rr, cc, num_samples] = get_indices_from_winsize(rows, cols, win_size);

% Note: patch_idx is a matrix
patch_idx = sub2ind([n_rows, n_cols], rr, cc);

% allocate space
samples = zeros(num_samples, n_filters * n_levels * n_pixels);

cols = 1:n_pixels;
for i_filter = 1:n_filters
    for level = 1:n_levels
        responses_d = reshape(responses(:,:,:,i_filter), [], n_levels);
		samples(:, cols) = reshape(responses_d(patch_idx, level), ...
								   [num_samples, n_pixels]);
        cols = cols + n_pixels;
	end
end


%% Test script
function test_script()
clc;

image_in = ceil(rand(256,256)*255);

inds = ceil(rand(100,1)*numel(image_in));
[rows,cols] = ind2sub(size(image_in), inds);
rows = min(max(rows, 2), size(image_in,1)-1);
cols = min(max(cols, 2), size(image_in,2)-1);

fprintf('Testing...');
old_samples = old_func(image_in, rows, cols, {'win_size', 3});
new_samples = func(image_in, rows, cols, {'win_size', 3});

success = all(old_samples(:) == new_samples(:));
if success, fprintf('success\n');
else        error('fail');
end

fprintf('Testing sample_pixel_data...');
old_samples = sample_pixel_data(image_in, rows, cols, {'win_size', 3});
new_samples = func(image_in, rows, cols, {'win_size', 3, ...
                                          'subtract_mean', true});

success = all(old_samples(:) == new_samples(:));
if success, fprintf('success\n');
else        error('fail');
end

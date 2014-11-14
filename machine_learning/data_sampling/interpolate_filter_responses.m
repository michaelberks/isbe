function [samples] = interpolate_filter_responses(responses, rows, cols, varargin)
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
    'subtract_mean', false,...
    'interp_method', 'cubic');

%Responses is an n_level cell, the top level of which is the full size
%image
win_size = args.win_size; 
[n_rows, n_cols, n_filters] = size(responses{1});
n_levels = length(responses);
n_pixels = win_size * win_size;


%First check if we need to rotate the filters or take the maximum response
if args.do_max || args.rotate
    for i_level = 1:n_levels
        %Reshape responses so the bands form columns of a matrix
        responses{i_level} = reshape(responses{i_level}, [], n_filters);

        %Work which band has maximum absolute response for each sample
        [d, max_band] = max(abs(responses{i_level}), [], 2); clear d; %#ok

        %Loop through each band rotating/selecting the max as appropriate
        for i_filter = 1:n_filters

            %Get which samples were maximal in this band
            band_idx = max_band==i_filter;

            if args.do_max
                %Note taking max in this way select maximum absolute response
                %whilst preserving sign
                responses{i_level}(band_idx,1) = responses{i_level}(band_idx,i_filter);
            else
                %We must be rotating (we'll assume the user hasn't being silly
                %enough to select do_max and rotate!) so shift such that
                %i_filter occupies the 1st column
                responses{i_level}(band_idx,:) = circshift(responses{i_level}(band_idx,:), [0 1-i_filter]);
            end
        end
        if args.do_max
            %Throw away the remaining bands
            responses{i_level}(:,2:n_filters) = [];
            n_filters = 1;
        end
        %Reshape responses
        responses{i_level} = reshape(responses, n_rows, n_cols, n_levels, n_filters);
    end
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
    
    %For level 1 we don't need to interpolate, just copy the responses
    filter_response = responses{1}(:,:,i_filter);
    samples(:, cols) = filter_response(patch_idx);
    cols = cols + n_pixels;
    
    for i_level = 2:n_levels
        %Interpolate responses
        
        scaling = 2^(i_level-1);
        offset = 1 + (i_level-1) / (2^i_level);
        filter_response = padarray(responses{i_level}(:,:,i_filter), [1 1], 'replicate');
		samples(:, cols) = interp2(filter_response,...
            offset+cc/scaling, offset+rr/scaling, args.interp_method);
        cols = cols + n_pixels;
    end
end

% Subtract mean if requested
if args.subtract_mean
    samples = bsxfun(@minus, samples, mean(samples, 2));
end

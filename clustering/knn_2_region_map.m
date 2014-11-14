function [dist_map] = knn_2_region_map(varargin)
%KNN_2_REGION_MAP *Insert a one line summary here*
%   [] = knn_2_region_map(varargin)
%
% KNN_2_REGION_MAP uses the U_PACKARGS interface function
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
% Created: 25-Jun-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'region1',... % the mandatory arguments
    'region2',...
    'num_samples2'}, ...
    'k', 10,...
    'win_size', 3, ...
    'num_levels', 6,...
    'do_max', 0,...
    'feature_type', 'all',...
    'max_size', 128,...
    'mask', []);

%Sample from region 2
[row2 col2] = size(args.region2);
sample_idx = randperm(row2*col2)';
sample_idx(args.num_samples2+1:end) = [];
[rows2 cols2] = ind2sub([row2 col2], sample_idx);
sample2 = sample_data(args.region2, rows2, cols2, args.win_size, args.num_levels, args.feature_type, args.do_max);

%compute number of parts we need to break image into
[row1 col1] = size(args.region1);
r_parts = ceil(row1 / args.max_size);
c_parts = ceil(col1 / args.max_size);

%preallocate distance map
dist_map = zeros(size(args.region1));

%Go through each segment
for rp = 1:r_parts
    for cp = 1:c_parts
        
        sr = 1 + (rp-1)*args.max_size;
        er = min(rp*args.max_size, row1);
        sc = 1 + (cp-1)*args.max_size;
        ec = min(cp*args.max_size, col1);
        
        %Get rows/cols subscripts for this part
        [part_cols part_rows] = meshgrid(sc:ec,sr:er);
        part_idx = sub2ind([row1 col1], part_rows, part_cols);
        
        %Check whether we've been given a mask to select specific pixels
        if ~isempty(args.mask)
            
            %Throw away pixels not belonging to the mask
            part_rows(~args.mask(part_idx)) = [];
            part_cols(~args.mask(part_idx)) = [];
            part_idx(~args.mask(part_idx)) = [];
            
            %check whether there's any pixels left to process
            if isempty(part_rows)
                continue;
            end
        end
        
        %Get sample points from this part
        sample1 = sample_data(args.region1, part_rows(:), part_cols(:), args.win_size, args.num_levels, args.feature_type, args.do_max);
        
        %Now get inter/intra sample distances and save in distance map
        dist_map(part_idx) = knn_2_sample_distance(sample1, sample2, args.k);
        
    end
end

function [training_data] =...
        sample_data(image_in, rows, cols, win_size, num_levels, feature_type, do_max)

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size
win_idx = -pad_w:pad_w;
num_pts = length(rows);
    
% Create DT-CWT of image
dt = dtwavexfm2b(image_in, num_levels); clear image_in;

%Make copies of sample rows and cols at positions of local window patch
rr = repmat(rows*ones(1,win_size) + ones(num_pts,1)*win_idx, 1, win_size);
cc = kron(cols*ones(1,win_size) + ones(num_pts,1)*win_idx, ones(1,win_size));
    
%Get interpolated dual-tree coefficients
dt_samples = dt_to_pixel_subset(dt, rr, cc); clear dt rr cc;

if do_max
    %get the maximum response across orientations
    dt_samples = squeeze(max(dt_samples, [], 3)); clear dt;
end
    
temp_samples=reshape(dt_samples, num_pts, []);
    
switch feature_type
    case 'all'
        training_data = [abs(temp_samples) angle(temp_samples)];
        
    case 'real'
        training_data = real(temp_samples);

    case 'mag'
        training_data = abs(temp_samples);

    case 'phase'
        training_data = angle(temp_samples);

    otherwise
        warning(['Feature type: ', args.feature_type, ' not recognised']); %#ok
        training_data = [abs(temp_samples) angle(temp_samples)];
end
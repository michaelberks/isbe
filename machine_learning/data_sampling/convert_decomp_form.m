function [data] = convert_decomp_form(data, old_decomp_args, new_decomp_args)
%CONVERT_DECOMP_FORM *Insert a one line summary here*
%   [data_out] = convert_decomp_form(data_in, decomp_args, new_args)
%
% Inputs:
%      data - *Insert description of input variable here*
%
%      decomp_args - *Insert description of input variable here*
%
%      new_args - *Insert description of input variable here*
%
%
% Outputs:
%      data - *Insert description of input variable here*
%
%
% Example:
%
% Notes: 
%
% See also:
%
% Created: 01-Oct-2012
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
new_decomp_args = u_packargs(new_decomp_args, '0',...
    'levels', [],...
    'bands', [], ...
    'win_size', 3,...
    'do_max', 0,...
    'rotate', 0,...
    'feature_type', 'all');

%Compute the number of samples in each data block for the old decomp type
[n_bands n_pixels n_levels] = get_sample_block_sizes(old_decomp_args);
n_samples = size(data, 1);

if (n_bands*n_pixels*n_levels) ~= size(data, 2)
    %Something has gone wrong
    error('Data size does not match old decomposition args');
end

%Reshape the data
switch old_decomp_args.decomp_type{1}
    case 'dt'
        
        %Check data is in full complex form 
        [data band_reduction] = convert_complex_representation_back(data, ...
            old_decomp_args.feature_type);
        n_bands = n_bands / band_reduction;
        
        %Reshape the data into a 4D vector
        data = reshape(data, n_samples, n_pixels, n_bands, n_levels);
        
        %Swap the last 2 dims (for DT only)
        data = permute(data, [1 2 4 3]);
        
    case {'gabor', 'gabori'}
        %Check data is in full complex form 
        [data band_reduction] = convert_complex_representation_back(data, ...
            old_decomp_args.feature_type);
        n_bands = n_bands / band_reduction;
        
        %Reshape the data into a 4D vector
        data = reshape(data, n_samples, n_pixels, n_levels, n_bands);
        
    otherwise %Don't have to mess around with complex forms
        %Reshape the data into a 4D vector
        data = reshape(data, n_samples, n_pixels, n_levels, n_bands);
        
end

%Get the indices for the pixel block given a new window size and extract
%these data
pixel_idx = get_win_idx_subset(old_decomp_args.win_size, new_decomp_args.win_size);
data = data(:,pixel_idx,:,:);
n_pixels = new_decomp_args.win_size^2;

%Now select the levels we want
if ~isempty(new_decomp_args.levels)
    if any(~ismember(new_decomp_args.levels, 1:n_levels))
        error('Levels selected for new decomposition are not valid');
    else
        data = data(:,:,new_decomp_args.levels,:);
        n_levels = length(new_decomp_args.levels);
    end
end

%Now select the bands we want
if ~isempty(new_decomp_args.bands)
    if any(~ismember(new_decomp_args.levels, 1:n_bands))
        error('Bands selected for new decomposition are not valid');
    else
        data = data(:,:,:,new_decomp_args.bands);
        n_bands = length(new_decomp_args.bands);
    end
end

%Now see if we need to rotate or do max
if new_decomp_args.do_max || new_decomp_args.rotate
    [data n_bands] = reorder_bands(data, new_decomp_args.do_max);
end

%Reshape the data so each feature vector is a single row
if strcmpi(old_decomp_args.decomp_type, 'dt')
        %Swap the last 2 dims (for DT only)
        data = permute(data, [1 2 4 3]);
end
data = reshape(data, n_samples, n_pixels*n_levels*n_bands);

%Finally complex types convert the complex form
if isfield(old_decomp_args, 'feature_type') && ...
    ~strcmpi(old_decomp_args.feature_type, new_decomp_args.feature_type)

    data = convert_complex_representation(data,...
            new_decomp_args.feature_type, new_decomp_args.win_size);
        
end
        

        

function [] = correct_mammo_dropoff(mammo_dir, mask_dir, save_dir, s_width)
%CORRECT_MAMMO_DROPOFF *Insert a one line summary here*
%   [] = correct_mammo_dropoff(mammo_dir, mask_dir, s_width)
%
% Inputs:
%      mammo_dir - *Insert description of input variable here*
%
%      mask_dir - *Insert description of input variable here*
%
%      s_width - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 03-Dec-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if nargin < 4
    s_width = 15;
end

mammo_list = dir([mammo_dir '/*.mat']);
mammo_names = get_mammo_info(mammo_list);

[mask_names missing_idx] =...
    match_mammo_names(mask_dir, mammo_names);
    
if ~isempty(missing_idx)
    display('warning: missing masks');
    mammo_list(missing_idx) = [];
    mask_names(missing_idx) = [];
end

num_mammos = length(mammo_list);
for ii = 1:num_mammos
    display(['processing ' num2str(ii) ' of ' num2str(num_mammos)]);
    
    mammo = double(u_load([mammo_dir '/'  mammo_list(ii).name]));
    mask = u_load([mask_dir '/' mask_names{ii}]);
    
    new_mask = mask;
    orig_area = sum(mask(:));
    new_area = orig_area;
    while new_area > 0.9*orig_area;
        new_mask = imerode(new_mask, strel('disk', 1));
        new_area = sum(new_mask(:));
    end
    
    Te = mean(mammo(new_mask));
    smooth_mammo = imfilter(mammo, fspecial('average', s_width));
    
    edge_idx = smooth_mammo < Te;
    mammo_processed = mammo;
    mammo_processed(edge_idx) = uint8(mammo(edge_idx) - smooth_mammo(edge_idx) + Te);
    mammo_processed(~mask) = uint8(Te); %#ok
    
    save([save_dir '/' mammo_list(ii).name], 'mammo_processed');
    clear new_mask mask mammo mammo_processed
end

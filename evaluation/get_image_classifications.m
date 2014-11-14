function [classifications image_n_pts] =...
    get_image_classifications(image_dir, class_dir, label_type)
%GET_IMAGE_CLASSIFICATIONS Compute ROC curve for a set of test images that have
%been classified
%   [classifications] = compute_image_orientation_errors(label_dir,class_dir,label_name)
%
% Inputs:
%      image_dir- directory containing test images with associated labels
%
%      class_dir- directory containing classification maps of the test images
%
%
% Outputs:
%
%
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if nargin < 3
    label_type = 'centre_line';
end

%Assume labels are sto
label_dir = [image_dir '/labels/'];

%get directory listing of input images
if ~isempty(label_dir) && ~strcmp(label_dir(end), filesep)
    label_dir = [label_dir filesep];
end
label_list = dir([label_dir, '*.mat']);

%get directory of classification images
if ~isempty(class_dir) && ~strcmp(class_dir(end), filesep)
    class_dir = [class_dir filesep];
end
class_list = dir([class_dir, '*.mat']);

%Check lists are same length
if length(class_list) ~= length(label_list)
    error('Number of labels and probability images differ');
end

%pre-allocate outputs
classifications = [];
image_n_pts = zeros(length(class_list),1);

%For each image
for ii = 1:length(class_list)
    
    %Load probability image
    class_map = load_uint8([class_dir, class_list(ii).name]);
    
    %Load image label
    s = load([label_dir label_list(ii).name]); %#ok
    
    switch label_type
        case 'centre_line'
            label = s.label_centre & s.label < 2;
            
        case 'all_line'
            label = s.label == 1;
            
        case 'not_centre'
            label = s.label == 1 & ~s.label_centre;

        otherwise
            warning('Label type not recognised, using centre line'); %#ok
            label = s.label_centre(:);
            
    end
    
    %Get image classifications at labels and add to main store
    classifications = [classifications; class_map(label)]; %#ok
    image_n_pts(ii) = sum(label(:));
    
end

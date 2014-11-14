function [] = classify_image_set2(job_id, level, image_dir, output_dir)
%
%
%
% Inputs:
%      fold- *Insert description of input variable here*
%
%      n_fold- *Insert description of input variable here*
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
% Created: 03-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Load in random forest and start and end indices
load([mberksroot, 'classification/data/spicules/rf_all_data_3_', zerostr(level,1), '_all_', zerostr(job_id,2), '.mat']);

%Get list of images
image_list = dir([mberksroot image_dir '*.mat']);

%Make output directory
mkdir([mberksroot output_dir]);

%3. Classify each image using forest
for ii = start_n:end_n
    image_in = u_load([mberksroot image_dir, image_list(ii).name]);
    [prob_spic] = classify_image(...
        'image_in', image_in, ...
        'forest', random_forest,...
        'forest_type', 'isbe_boot'); %#ok
    
    %Save the image in the output directory
    save([mberksroot output_dir 'prob_spic', zerostr(ii,3), '.mat'], 'prob_spic');
end


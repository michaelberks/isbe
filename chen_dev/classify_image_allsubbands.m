function [probability_image, labels] = classify_image_allsubbands(image_in, forest, forest_type, num_levels, win_size)
%CLASSIFY_IMAGE *Insert a one line summary here*
%   [probability_image] = classify_image(image_in,forest)
%
% Inputs:
%      image_in - Image to classify
%
%      forest - Random forest classifier (model)
%
%
% Outputs:
%      probability_image - Pixel-wise probability of belonging to a bar
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 27-Jan-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

if nargin < 3
    forest_type = 'isbe';
end
if nargin < 5
    win_size = 3;
end
if nargin < 4
    num_levels = 4;
end
if nargin < 5
    win_size = 3;
end

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size
size_sample_vector = 2*win_size*win_size*num_levels*6;
win_idx = -pad_w:pad_w;

[ROW COL] = size(image_in);

%Create storage for test data
test_data = zeros(ROW*COL, size_sample_vector);

% Create DT-CWT of image
dt = dtwavexfm2b(image_in,num_levels);

%interpolation DT-CWT to full pixel grid
dt_inter = dt_to_full_image(dt); clear dt;

%pad the edges of subbands and label
pad_dt = padarray(dt_inter, [pad_w pad_w], 'replicate'); clear dt_inter;

curr_sample = 1;
%loop through image pixels selecting points
for cc = (1:COL) + pad_w
    for rr = (1:ROW) + pad_w

        %Extract and sample vector
        sample_vector = pad_dt(rr+win_idx, cc+win_idx, :);

        %Convert sample vector to mag/phase form and save in test data
        test_data(curr_sample,:)= [abs(sample_vector(:))' angle(sample_vector(:))'];
        curr_sample = curr_sample + 1;
    end
end
clear pad_dt;

%Now classify the data using either ours or Breiman's code
if strcmpi(forest_type, 'breiman')
    [labels votes] = classRF_predict(test_data, forest);
    probability_image = reshape(votes(:,1) / forest.ntree, size(image_in));
    labels = reshape(labels, size(image_in));
elseif strcmpi(forest_type, 'isbe')
    [labels votes] = mb_random_forest_class_predict(forest, test_data);
    probability_image = reshape(votes(:,1) / length(forest.trees), size(image_in));
    labels = reshape(labels, size(image_in));
elseif strcmpi(forest_type, 'isbe_boot')
    [labels votes] = mb_random_forest_class_predict_boot(forest, test_data);
    probability_image = reshape(votes(:,1) / length(forest.trees), size(image_in));
    labels = reshape(labels, size(image_in));
else
    error('Type of random forest not recognised');
end


%probability_image = test_data;

function compute_rf_pair_classification(idx, data_type)
%COMPUTE_RF_PAIR_CLASSIFICATION function to batch process RF pair
%classifcation on hydra
%   compute_rf_pair_classification(idx, data_type)
%
% Inputs:
%      idx - specify pair number to process (from directory listing)
%
%      data_type - string to data type we're processing (e.g
%      'abnormals')
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
% Created: 02-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

display(['compute_rf_pair_classification running: ' datestr(now)]);
%--------------------------------------------------------------------------

idx = idx + 20;

rand('twister', sum(100*clock));

%Set mammogram/segmentation directories
mammo_dir = ['/data/zchen/mammograms/2004_screening/' data_type '/mat/']; %
seg_dir = ['/data/zchen/segmentations/2004_screening/' data_type '/']; %

%Get list of CC mammogrmas (note we're assuming we have strict pairs of
%mammograms in specified dirs)
r_list = dir([mammo_dir '*RCC*.mat']);
l_list = dir([mammo_dir '*LCC*.mat']);

%Load in left/right mammograms
forest_args.image1 = u_load([mammo_dir r_list(idx).name]);
forest_args.image2 = fliplr(u_load([mammo_dir l_list(idx).name]));

%Load in segmentations
seg_r = u_load([seg_dir r_list(idx).name(1:end-4) '_segmentation.mat']);
seg_l = u_load([seg_dir l_list(idx).name(1:end-4) '_segmentation.mat']);

%Resize segmentations
seg_r.breast_border = segment_breast_resize(size(forest_args.image1), seg_r);    
seg_l.breast_border = segment_breast_resize(size(forest_args.image2), seg_l);

%create masks of breast region for each mammograms
forest_args.mask1 = roipoly(forest_args.image1, seg_r.breast_border(:,1), seg_r.breast_border(:,2));
sub_mask1 = false(size(forest_args.mask1));
sub_mask1(2:2:end, 2:2:end) = true;
forest_args.mask1 = forest_args.mask1 & sub_mask1; clear sub_mask1;

forest_args.mask2 = fliplr(roipoly(forest_args.image2, seg_l.breast_border(:,1), seg_l.breast_border(:,2)));
sub_mask2 = false(size(forest_args.mask2));
sub_mask2(2:2:end, 2:2:end) = true;
forest_args.mask2 = forest_args.mask2 & sub_mask2; clear sub_mask2;

forest_args.num_levels = 4; %Note that the full mammograms are 1/2 resolution compared to previously

forest_args.win_size = 3;

forest_args.do_max = 0;

forest_args.num_train = 5e4;

forest_args.feature_type = 'all';

forest_args.d = 20;

forest_args.split_min = 200;

forest_args.tree_dir = ['/data/zchen/contralateral_images_binary/', data_type, '/'];

forest_args.save_path = ['/data/zchen/contralateral_images_binary/', data_type, '/random_forest', zerostr(idx, 3), '.mat'];

forest_args.n_trees = 100;

forest_args.do_test1 =1; 

forest_args.do_test2 =1; 

[random_forest] = mb_random_forest_two_images(forest_args);

%--------------------------------------------------------------------------
display(['compute_rf_pair_classification completed: ' datestr(now)]);
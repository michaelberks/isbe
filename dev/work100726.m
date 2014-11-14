image_num = 1;

rand('twister', sum(100*clock));


% sample_args.pair_name = 'M:\asymmetry_project\data\contralateral\002LCC';

image_dir = 'C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\';

image_list = dir([image_dir, '*.mat']);

sample_args.pair_name = [image_dir, image_list(image_num).name];

sample_args.win_size = 1;

sample_args.num_levels = 5;

sample_args.do_max = 0;

sample_args.num_train = 5e3;

sample_args.feature_type = 'ilp';

% 
forest_args.sampling_method = 'sample_contralateral_train';

forest_args.sampling_method_args = sample_args;

forest_args.d = 20;

forest_args.split_min = 200;

forest_args.tree_dir = ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\trees\', zerostr(image_num, 3), '/'];

forest_args.save_path = ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\image_pair', zerostr(image_num, 3), '/random_forest.mat'];

forest_args.do_err = 1;

forest_args.n_trees = 2;

%%
%build forest
[random_forest abnormal_roi_abnormalvotes normal_roi_abnormalvotes abnormal_roi_totalvotes normal_roi_totalvotes] = mb_random_forest_contra_train(forest_args);

save([forest_args.tree_dir, 'votes_image.mat'], 'abnormal_roi_abnormalvotes', 'normal_roi_abnormalvotes', 'abnormal_roi_totalvotes', 'normal_roi_totalvotes');




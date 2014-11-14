
for ii = 1:160
    [image_out, label, pars] = create_bar_image;
    
    parameters(ii) = pars; %#ok
    
    save(['C:\isbe\dev\classification\data\bg+bar_128\bar', zerostr(ii,3)],...
        'image_out', 'label');
    display(['Saved image ', num2str(ii)]);
end
%%
generate_bar_training_data(...
    'num_images', 100,...
    'image_dir', 'C:\isbe\dev\classification\data\bg+bar_128_test\',...
    'num_levels', 4,...
    'compute_dt', 1);

%%
profile on;
[training_data training_labels] = sample_image_training_data(...
    'num_samples', 4e5,...
    'image_dir', 'C:\isbe\dev\classification\data\bg+bar_128\',...
    'total_samples', 128*128*160);
profile viewer
%%
profile on;
sampling_method_args.num_samples = 4e5;
sampling_method_args.image_dir = 'C:\isbe\dev\classification\data\bg+bar_128\';
sampling_method_args.total_samples = 128*128*160;

forest_args.sampling_method = 'sample_image_training_data';
forest_args.sampling_method_args = sampling_method_args;
forest_args.d = 10;
forest_args.n_trees = 1;
forest_args.tree_dir = 'C:\isbe\dev\classification\rf\rf_test';
[random_forest] = mb_random_forest_class_train(forest_args);
profile viewer;
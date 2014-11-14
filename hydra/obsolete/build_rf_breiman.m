function [] = build_rf_breiman(image_dir, n_samples, n_trees, d)
%BUILD_RF *Insert a one line summary here*
%   [] = build_rf(job_idx,image_dir)
%
% Inputs:
%      job_idx- *Insert description of input variable here*
%
%      image_dir- *Insert description of input variable here*
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

sampling_method_args.num_samples = n_samples;
sampling_method_args.image_dir = [mberksroot, 'classification/data/' image_dir];
sampling_method_args.total_samples = 128*128*160;

[training_data training_labels] = sample_image_training_data(sampling_method_args);
tic; rf_brei_bar = classRF_train(training_data, training_labels, n_trees, d); toc;
save([mberksroot, 'classification/rf/rf_breiman.mat'], 'rf_brei_bar');

for ii = 1:100
    load([mberksroot, 'classification/data/bg+bar_128_test/bar', zerostr(ii,3), '.mat']);
    [probability_image] = classify_image(image_out, rf_brei_bar, 'breiman'); %#ok
    save([mberksroot, 'classification/data/bg+bar_128_test_results/bar', zerostr(ii,3), '_results.mat'], 'probability_image');
end



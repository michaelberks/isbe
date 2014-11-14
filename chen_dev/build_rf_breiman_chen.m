% function [] = build_rf_breiman(model_dir, num_samples, n_trees, d)
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

if strcmpi(computer,'PCWIN') |strcmpi(computer,'PCWIN64')
    %% windows path
    bg_dir='E:\DTCWTmini\data\normal_smooth128\';
    save_dir='E:\DTCWTmini\data\cpiculeimage\';
   image512_dir='E:\DTCWTmini\data\normal_512\';
else
    %% Linux path
    bg_dir='/home/zchen/data/normal_smooth128/';
    save_dir='/home/zchen/data/cpiculeimage/';
    image512_dir='/home/zchen/data/normal_512/';
end
   
num_samples = 400000;
n_trees = 200;
d = 10;
sampling_method_args.num_samples = num_samples;
sampling_method_args.bg_dir = bg_dir; 
sampling_method_args.save_path = save_dir; 
[training_data training_labels parameters] = sample_spicule_training_data_chen(sampling_method_args);
tic; model = classRF_train(training_data, training_labels, n_trees, d); toc;
% save([save_dir, 'model_spicule.mat'], 'model');
% model_h1 = struct;
% model_h1.nrnodes=model.nrnodes;
% model_h1.ntree=model.ntree;
% model_h1.xbestsplit=model.xbestsplit;
% model_h1.classwt=model.classwt;
% model_h1.cutoff=model.cutoff;
% model_h1.treemap=model.treemap;
% model_h1.nodestatus=model.nodestatus;
% 
% model_h2 = struct;
% model_h2.nodeclass=model.nodeclass;
% model_h2.bestvar=model.bestvar;
% model_h2.ndbigtree=model.ndbigtree;
% model_h2.mtry=model.mtry;
% model_h2.orig_labels=model.orig_labels;
% model_h2.new_labels=model.new_labels;
% model_h2.nclass=model.nclass;
% model_h2.outcl=model.outcl;
% model_h2.counttr=model.counttr;
% model_h2.proximity=model.proximity;
% model_h2.localImp=model.localImp;
% model_h2.importance=model.importance;
% model_h2.importanceSD=model.importanceSD;
% model_h2.errtr=model.errtr;
% model_h2.inbag=model.inbag;
% model_h2.votes=model.votes;
% model_h2.oob_times=model.oob_times;
% save([save_dir, 'model_spicule_breiman_h1.mat'], 'model_h1');
% save([save_dir, 'model_spicule_breiman_h2.mat'], 'model_h2');
% 

for ii = 1:100
    test_image=u_load([save_dir, 'test_image', int2str(ii), '.mat']);
     [probability_image] = classify_image(test_image, forest_model, 'breiman'); %#ok
    predict(ii).probability_image = probability_image; 
    disp(ii);
end
    save([save_dir, 'predict_spicule_breiman'], 'predict');

bglist=[3, 6, 26, 49, 64];  % normal_512
%predict the test data
for ii = 1:length(bglist)
    test_image=u_load([image512_dir, 'bg', zerostr(bglist(ii),3), '.mat']);
    [probability_image] = classify_image(test_image, forest_model, 'isbe'); %#ok
    predict(ii).probability_image = probability_image;
    disp(ii);
end
save([save_path, 'predict_spicule_breiman_512'], 'predict');

%function analyse_aam_results
%% 1) First get the apex candidates and model scores from the output files of the AAM fit
rsa_dir = 'rsa_study/';
apex_gt_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_gt/'];
aam_dir = [nailfoldroot 'data/' rsa_dir 'test/aam/'];
create_folder([aam_dir 'candidates']);

test_list = dir([apex_gt_dir '*_gt.mat']);
num_images = length(test_list);

for i_im = 1:num_images;
    
    im_num = test_list(i_im).name(1:6);
    display(['Processing image ' num2str(i_im) ' of ' num2str(num_images)]);
    
    
    %load([apex_gt_dir im_num '_gt.mat'], 'distal_xy', 'distal_width', 'non_distal_xy', 'undefined_xy');
    %load([centre_dir im_num '_vc.mat']);
    
    %---------------------------------------------
    apex_list = dir([aam_dir im_num '/*aam.mat']);
        
    num_candidates = length(apex_list);
    candidate_xy = zeros(num_candidates,2);
    candidate_scores = zeros(num_candidates,1);

    for i_ap = 1:num_candidates
        apex_candidate = u_load([aam_dir im_num '/' apex_list(i_ap).name]);
        candidate_xy(i_ap,:) = apex_candidate.vessel_xy(16,:) + ...
            [apex_candidate.sc apex_candidate.sr] - 1;
        candidate_scores(i_ap,:) = apex_candidate.model_score;
    end
    save([aam_dir 'candidates/' im_num '_candidates.mat'], 'candidate_xy', 'candidate_scores');

end
%% 2) Now we can run through all the images, work which apices were detected
% We also analyse:
%   - What measure detected each apex
%   - What % of detections were true
%   - For apices we missed, why were they discarded
make_detection_results_struc(...
    'results_name', 'aam_results_', ...
    'candidates_dir',    'C:\isbe\nailfold\data\rsa_study\test\aam\candidates\',... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test/apex_gt/'],...
    'results_dir', [nailfoldroot 'data/rsa_study/test/results/'],...
    'selected_gt', [],...
    'selected_candidates', []);

%Results codes:
% Model fit following shape prior/cross-correlation matching: aam_results_20131031T161927
% Model fit following RF offset predictions: aam_results_20131106T110627

%% 3) Start doing some analysis
analyse_detection_results(...
    'results_name', 'aam_results_20131119T094236',...
    'candidates_dir',    'C:\isbe\nailfold\data\rsa_study\test\aam\candidates\',... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test/apex_gt/'],...
    'results_dir', [nailfoldroot 'data/rsa_study/test/results/'],...
    'selected_images', [],...
    'use_only_gradeable', 1,...
    'min_num_markers', 2,...
    'max_missing_markers', inf,...
    'overall_summary', 1,...
    'summary_by_grade', 1,...
    'summary_by_shape', 1,...
    'summary_by_size', 1,...
    'compute_rocs', 0,...
    'analysis_by_width', 0);
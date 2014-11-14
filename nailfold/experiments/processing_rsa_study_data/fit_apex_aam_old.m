function [] = fit_apex_aam_old(start_idx, end_idx)
%ANALYSE_APEX_IMAGE_PROPERTIES *Insert a one line summary here*
%   [] = analyse_apex_image_properties()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if isdeployed    
    task_idx = unixenv('SGE_TASK_ID', 1);
    start_idx = 20*(task_idx-1) + 1;
    end_idx = 20*task_idx;
else
    if ~exist('start_idx', 'var'); start_idx = 1; end
    if ~exist('end_idx', 'var'); end_idx = 10; end
end

%Set up image lists
if ispc
    rsa_dir = 'rsa_study/';
else
    rsa_dir = [];
end

%% 1) Fit an AAM to the candidates

corr_scores_dir = [nailfoldroot 'data/' rsa_dir 'test/corr_scores/'];
shape_scores_dir = [nailfoldroot 'data/' rsa_dir 'test/shape_scores/'];
image_dir = [nailfoldroot 'data/' rsa_dir 'test/images/'];
aam_dir = [nailfoldroot 'data/' rsa_dir 'test/aam/'];
model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];

do_initialisation = 1 < 0;
do_aam_fit = 1 > 0;

test_list = dir([corr_scores_dir '*_cs.mat']);
num_images = length(test_list);

if do_initialisation
    %Load in model templates
    load([model_dir 'apex_templates.mat']);
    load([model_dir 'mean_shape.mat'], 'mean_shape');
    g2_template = apex_template.g2;
    clear apex_template;    
end

%
if ~exist('start_idx', 'var'); start_idx = 1; end;
if ~exist('end_idx', 'var'); end_idx = num_images; end;

%Loop though each image
for i_im = start_idx:end_idx
    
    im_num = test_list(i_im).name(1:6);
    display(['Processing image ' num2str(i_im) ' of ' num2str(num_images)]); 
    
    if do_initialisation
        load([shape_scores_dir im_num '_ss.mat'], 'maxima_shape_xy');
        load([corr_scores_dir im_num '_cs.mat'], 'maxima_corr_xy');
        maxima_shape_xy(101:end,:) = []; 
        candidates_xy = [maxima_shape_xy; maxima_corr_xy];

        [candidates_xy] = ...
            apply_local_exclusion(candidates_xy,...
            [zeros(size(maxima_shape_xy,1),1); ones(size(maxima_corr_xy,1),1)], 50);

        initialise_aam_candidates([image_dir im_num '.mat'], candidates_xy, ...
            'g2d_template', g2_template,...
            'max_num_candidates', inf,...
            'thetas', -15:3:15,...
            'scales', [0.8:0.2:2 3 4],...
            'mean_shape', mean_shape,...
            'aam_dir', [aam_dir im_num '/']);
    end
    if do_aam_fit
    
        fit_aam_to_candidates(...
            'aam_dir', [aam_dir im_num '/'],...
            'aam_exe', 'ncm_sandpit_mb',...
            'aam_path', [model_dir 'aam/orig/2/vessel_apex_orig.smd'],...
            'g2d_sigma', 4, ...
            'delete_candidate_patches', 0);
    end
       
end
        
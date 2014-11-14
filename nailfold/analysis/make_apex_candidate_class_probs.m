function [] = make_apex_candidate_class_probs(varargin)
%EXTRACT_VESSEL_CENTRES *Insert a one line summary here*
%   [] = extract_vessel_centres()
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
args = u_packargs(varargin, 0, ... % the user's input
    {'image_names'},        ...
    'model_id',             [],...
    'model_root',           [nailfoldroot,'models/apex'], ...
    'model_name',           'class_map',...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'displacement_dir',     'bob',...
    'candidates_dir',       'bob',...
    'label_dir',            'bob',...
    'plot',                 0,...
    'fig_dir',              []);

%Set up model name
if isempty(args.model_id)
    model_id = datestr(now, 30);
else
    model_id = args.model_id;
end
model_dir = [args.model_root '/final_MAP/' model_id '/'];
create_folder(model_dir);

%Set up dir paths
displacement_dir = [args.data_dir '/' args.displacement_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];
label_dir = [args.data_dir '/' args.label_dir '/'];

num_images = length(args.image_names); 

%Loop through images loading in the relevant data
all_candidate_scores = [];
all_candidate_disp = [];
all_candidate_labels = [];
all_candidate_class = false(0,1);
for i_im = 1:num_images
    im_name = args.image_names{i_im} ;
    
    load([candidates_dir im_name '_candidates.mat'], 'candidate_rescores');
    load([label_dir im_name '_label.mat'], 'candidates_class', 'candidate_labels');
    load([displacement_dir im_name '_dis.mat'], 'candidate_displacements');
    
    all_candidate_scores = [all_candidate_scores; candidate_rescores]; %#ok
    all_candidate_class = [all_candidate_class; candidates_class]; %#ok
    all_candidate_labels = [all_candidate_labels; candidate_labels]; %#ok
    all_candidate_disp = [all_candidate_disp; candidate_displacements]; %#ok
end

%Now buildmodel of class probs
combined_dist = cell(1,3);
joint_conditional_probs = zeros(100,100,3);
posterior_class_probs = zeros(100,100,3);
[xx yy] = meshgrid(linspace(0, 1, 100), linspace(-200, 200, 100));
grid_xy = [xx(:) yy(:)]; clear xx yy;    
conditional_class_probs = zeros(100,100,3);
conditional_total_prob = zeros(100,100) + eps;
prior_class_probs = zeros(3,1);

for i_cl = 1:3
    cl_idx = all_candidate_labels==i_cl & all_candidate_disp < inf;   
    
    %Estimate distribution of each class
    combined_dist{i_cl} = build_2d_kernel_distribution(...
        [all_candidate_scores(cl_idx), all_candidate_disp(cl_idx)], grid_xy, []);
    joint_conditional_probs(:,:,i_cl) = reshape(combined_dist{i_cl}.D_f,100,100);
    %Compute class prior
    prior_class_probs(i_cl) =...
        sum(cl_idx) / size(cl_idx,1);

    %Compute class conditional probs
    conditional_class_probs(:,:,i_cl) = ...
        joint_conditional_probs(:,:,i_cl) * prior_class_probs(i_cl);

    %Add to total probs to get normalising factor
    conditional_total_prob = conditional_total_prob +...
        conditional_class_probs(:,:,i_cl);
    
end    

for i_cl = 1:3
    %Now we can compute posterior probs for each clas
    posterior_class_probs(:,:,i_cl) = ...
        conditional_class_probs(:,:,i_cl) ./ conditional_total_prob;
end
   
%Finally we can compute the MAP
[max_posterior_probs posterior_class] = max(posterior_class_probs, [], 3); 

%Save the ouput model structure
class_map.x = reshape(combined_dist{1}.x,100,100);
class_map.y = reshape(combined_dist{1}.y,100,100);
class_map.joint_conditional_probs = joint_conditional_probs;
class_map.conditional_class_probs = conditional_class_probs;
class_map.posterior_class_probs = posterior_class_probs;
class_map.post_class = posterior_class - 1;
class_map.post_probs = max_posterior_probs;

save([model_dir args.model_name '.mat'], 'class_map');

%Plot stuff if we have to...
if args.plot
    colors = 'rgb';

    figure; hold on;
    for i_cl = 1:3
        mesh(reshape(combined_dist{i_cl}.x,100,100),  reshape(combined_dist{i_cl}.y,100,100),...
            reshape(combined_dist{i_cl}.D_f,100,100), 'edgecolor', colors(i_cl));
    end
    xlabel('X (candidate scores)');
    ylabel('Y (candidate displacements)');
    zlabel('p(X&Y|C_i)');
    title('Joint class pdfs over candidate displacement and score');
    if ~isempty(args.fig_dir)
        saveas(gcf, [args.fig_dir 'joint_pdfs_.fig']);
    end

    figure; hold on;
    for i_cl = 1:3

        mesh(reshape(combined_dist{i_cl}.x,100,100),  reshape(combined_dist{i_cl}.y,100,100),...
            posterior_class_probs(:,:,i_cl), 'edgecolor', colors(i_cl));
    end

    xlabel('X (candidate scores)');
    ylabel('Y (candidate displacements)');
    zlabel('p(Ci|X&Y)');
    title('Posterior class probabilities over candidate displacement and score');
    if ~isempty(args.fig_dir)
        saveas(gcf, [fig_dir 'posterior_pdfs_.fig']);
    end

    figure; 
    mesh(...
        reshape(combined_dist{i_cl}.x,100,100),...
        reshape(combined_dist{i_cl}.y,100,100), ...
        conditional_total_prob);
    xlabel('X (candidate scores)');
    ylabel('Y (candidate displacements)');
    zlabel('p(X&Y)');
    title('Joint PDFs over candidate displacement and score');
    if ~isempty(args.fig_dir)
        saveas(gcf, [fig_dir 'data_pdfs.fig']);
    end

    figure; imgray(posterior_class); colormap(eye(3));
    xlabel('X (candidate scores)');
    ylabel('Y (candidate displacements)');
    title('Class map of MAP');
    
    %Do the marginal stuff...
    %figure;
    %subplot(1,2,1); hold all;
    %for i_cl = 1:3
    %end
    %for i_cl = 1:3
    %end
end
    
    

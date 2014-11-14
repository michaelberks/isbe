function [orientation_errors mean_image_errors line_contrasts] =...
    compute_image_orientation_errors(image_dir, ori_dir, label_type, resize, do_save)
%COMPUTE_IMAGE_ORIENTATION_ERRORS Compute ROC curve for a set of test images that have
%been classified
%   [orientation_errors] = compute_image_orientation_errors(label_dir,ori_dir,label_name)
%
% Inputs:
%      image_dir- directory containing test images with known orientations
%
%      ori_dir- directory containing orientation maps of the test images
%
%
% Outputs:
%      mean_errors
%
%      std_errors
%
%      total_samples
%
%
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if nargin < 5
    do_save = 1;
end
if nargin < 4
    resize = 0;
end
if nargin < 3
    label_type = 'centre_line';
end

%Assume labels are sto
label_dir = [image_dir '/labels/'];

%get directory listing of input images
if ~isempty(label_dir) && ~strcmp(label_dir(end), filesep)
    label_dir = [label_dir filesep];
end
label_list = dir([label_dir, '*.mat']);

%get directory of probability images
if ~isempty(ori_dir) && ~strcmp(ori_dir(end), filesep)
    ori_dir = [ori_dir filesep];
end
prob_list = dir([ori_dir, '*.mat']);

%Check lists are same length
if length(prob_list) ~= length(label_list)
    error('Number of labels and probability images differ');
end

%pre-allocate outputs
orientation_errors = [];
if nargout > 1
    mean_image_errors = zeros(length(prob_list),1);
end
do_contrasts = 0;
if nargout > 2
    line_contrasts = [];
    
    %Try and load the parameters file to workout the line intensity in each
    %image
    if exist([image_dir '/parameters/parameters.mat'], 'file')
        parameters = u_load([image_dir '/parameters/parameters.mat']);
        do_contrasts = 1;
    else
        display('No parameter file found, image contrasts will not be returned');
    end
end  

%For each image
for ii = 1:length(prob_list)
    
    %Load probability image
    orientation_map = load_uint8([ori_dir, prob_list(ii).name]);
    
    if resize
        orientation_map = imresize(orientation_map, resize, 'bilinear');
    end
    
    %Check if complex - if not assume map of degrees on range [0 180] and
    %do nothing
    do_spread = false;
    if ~isreal(orientation_map)
        %Complex - need to get arg of map then convert to degrees. Also get
        %the magnitude of the map as this measures the spread of
        %predictions of the orientations
        orientation_spread = abs(orientation_map);
        orientation_map = mod(180*angle(orientation_map)/pi,180);
        do_spread = 1;
    else
        orientation_map = mod(180*orientation_map/pi,180);
    end
    
    %Load image label
    s = load([label_dir label_list(ii).name]); %#ok
    
    switch label_type
        case 'centre_line'
            label = s.label_centre & s.label < 2;
            
        case 'all_line'
            label = s.label == 1;
            
        case 'not_centre'
            label = s.label == 1 & ~s.label_centre;

        otherwise
            warning('Label type not recognised, using centre line'); %#ok
            label = s.label_centre(:);
            
    end
    
    %Get difference (modulo 180) between orientations for labelled pixels
    image_errors = ...
        mb_mod(orientation_map(label) - s.label_orientation(label), 180);
    
    %Add the image errors to the overall list - for complex input also add
    %the measure of orientation spread in the 2nd column
    if do_spread
        orientation_errors = [orientation_errors;...
            image_errors orientation_spread(label)]; %#ok
    else
        orientation_errors = [orientation_errors; image_errors]; %#ok
    end
    
    if nargout > 1
        mean_image_errors(ii) = mean(abs(image_errors));
    end
    
    %If we're computing contrasts we need to construct the contrast image
    if do_contrasts
        contrast_im = zeros(size(label));
        for jj = 1:length(parameters(ii).curr_para)
            
            %Make each bar from saved parameters
            [bar_im bar_label] = create_sin_curve(...
                parameters(ii).curr_para(jj).halfwidth,...
                parameters(ii).curr_para(jj).contrast,...
                parameters(ii).curr_para(jj).radius,...
                parameters(ii).curr_para(jj).orientation,...
                parameters(ii).curr_para(jj).squash,...
                parameters(ii).curr_para(jj).row,...
                parameters(ii).curr_para(jj).col,...
                parameters(ii).curr_para(jj).centre_x,...
                parameters(ii).curr_para(jj).centre_y);
            
            %Add to the contrast image (don't need to worry about
            %intensities at the junctions as these won't be sampled from
            %contrast_im = contrast_im + bar_im;
            contrast_im(bar_label) = parameters(ii).curr_para(jj).contrast;
        end
        %Now sample the contrasts
        line_contrasts = [line_contrasts; contrast_im(label)]; %#ok
    end       
end


if do_save
    mkdir([ori_dir label_type 'errors/']);
    save([ori_dir label_type 'errors/ori_errors.mat'], 'orientation_errors');
end

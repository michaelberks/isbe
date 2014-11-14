function [training_data training_labels] = sample_image_training_data(varargin)
%SAMPLE_IMAGE_TRAINING_DATA *Insert a one line summary here*
%   [] = sample_image_training_data(varargin)
%
% SAMPLE_IMAGE_TRAINING_DATA uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 28-Jan-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
			 {'num_samples',... % the mandatory arguments
             'image_dir',...
             'total_samples'}, ...
             'win_size', 3,...
             'num_levels', 4,...
             'image_list', [], ...
             'image_type', '.mat',...
             'save_path', []);

win_size = args.win_size; % the size of scanning window
num_levels = args.num_levels; %number of levels of DT-CWT from which to extract coefficients
num_samples = args.num_samples; %total number of feature vectors to sample
total_samples = args.total_samples; %total number of feature vectors in global population
image_dir = args.image_dir; %directory in which images are stored
image_list = args.image_list;

if ~isempty(image_dir) && ~strcmp(image_dir(end), filesep)
    image_dir = [image_dir filesep];
end
if isempty(image_list)
    image_list = dir([image_dir, '*', args.image_type]);
end

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size
size_sample_vector = 2*win_size*win_size*num_levels;
win_idx = -pad_w:pad_w;

%Get sample indices for each image and save
%sample_idx = sort(randsample(total_samples, num_samples));
sample_idx = randperm(total_samples);
sample_idx(num_samples+1:end) = [];
sample_idx = sort(sample_idx);


if ~isempty(args.save_path)
    mkdir(args.save_path);
    save([args.save_path, '_sample_idx'], 'sample_idx');
end

%Get list of images in image dir
%Create data
training_data = zeros(num_samples, size_sample_vector);
training_labels = zeros(num_samples, 1);

curr_sample = 1;
%loop through each image sampling data
for kk = 1:length(image_list)
    
    %load in image
    S = load([image_dir image_list(kk).name]);
    image_in = S.image_out;
    label = S.label;
    max_dt = S.max_dt;
    clear S;
    
    %get size of image
    [ROW COL] = size(image_in);
    
    %See whether the loaded data contains a pre-computed DT-CWT, if not
    %we'll need to compute it now
    if isempty(max_dt)
        % Create DT-CWT of image
        dt = dtwavexfm2(image_in,num_levels,'near_sym_b','qshift_b');

        %interpolate DT-CWT to full pixel grid
        dt_inter = dt_to_full_image(dt); clear dt;

        %get the maximum response across orientations
        max_dt = squeeze(max(dt_inter, [], 3)); clear dt;
    end

    %pad the edges of subbands and label
    pad_dt = padarray(max_dt, [pad_w pad_w], 'replicate'); clear max_dt;
    pad_label = padarray(label, [pad_w pad_w], 'replicate');

    %Get indices for this image and convert to row,col
    image_idx = sample_idx(sample_idx <= ROW*COL);
    num_samples_image = length(image_idx);

    [rows cols] = ind2sub([ROW COL], image_idx);

    %Remove samples from main sample indices and subtract image size
    sample_idx(1:num_samples_image) = [];
    sample_idx = sample_idx - ROW*COL;

    %loop through image pixels selecting points
    for ii = 1:num_samples_image

        rr = rows(ii) + pad_w;
        cc = cols(ii) + pad_w;

        %Extract label and sample vector
        sample_label = pad_label(rr, cc);
        sample_vector = pad_dt(rr+win_idx, cc+win_idx, :);

        %Convert sample vector to mag/phase form and save in main data
        training_data(curr_sample,:) = [abs(sample_vector(:))' angle(sample_vector(:))'];
        %training_data(curr_sample,1:end/2) = sample_vector(:).';
        training_labels(curr_sample,:)= sample_label;

        %increment sample count
        curr_sample = curr_sample + 1;
    end
    clear pad_dt;
end
%Convert into magnitude and phase
% training_data(:,1+end/2:end) = angle(training_data(:,1:end/2));
% training_data(:,1:end/2) = abs(training_data(:,1:end/2));
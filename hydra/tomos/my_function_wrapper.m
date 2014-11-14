function my_function_wrapper(sge_idx)
%Wrapper function that takes the $SGE_IDX batch job variable and uses this
%to construct arguments to pass to my_function

%e.g. assume id_idx = 1,...,5, and use this to call different datasets
%labelled dataset_1, dataset_2, etc.
data_file_path = ['/san/staff/student/mberks/data/dataset_', num2str(sge_idx)];

% Possibly set up some other arguments for my functions
other_args.a = 1;

% Can make this conditional on sge_idx - so for example set a particular parameter to
% 0 for the odd indexed datasets, or 1 for the even

if rem(sge_idx, 2)
    other_args.b = 0;
else
    other_args.b = 1;
end

% Now call my_function
model = my_function(data_file_path, other_args); %#ok

%Can use id_idx to save outputs to different location
model_file_path = ['/san/staff/student/mberks/models/model_', num2str(sge_idx)]; %#ok

%save(model_file_path, 'model');

%tidy up and quit matlab if we're running on hydra (as opposed to testing
%on PC)
clear;
if ~ispc, exit; end



function args_out = remove_prefix(prefix, args_in)
% Parse a structure of arguments, removing the specified prefix from any of
% the field names. e.g. 'boost_n_levels' becomes 'n_levels' if
% prefix='boost_'.

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

args_out = func(prefix, args_in);


%% The function
function args_out = func(prefix, args_in)

fnames = fieldnames(args_in);
prefix_length = length(prefix);

for i = 1:length(fnames)
    if strcmp(fnames{i}(1:prefix_length), prefix)
        args_out.(fnames{i}(prefix_length+1:end)) = args_in.(fnames{i});
    else
        args_out.(fnames{i}) = args_in.(fnames{i});
    end
end


%% Test script
function test_script()
clc;

s.boost_n_levels = 10;
s.boost_shrinkage = true;
s.my_name = 'phil';

s2 = remove_prefix('boost_',s);
disp(s2);
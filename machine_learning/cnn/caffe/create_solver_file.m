function [] = create_solver_file(solver_filename, varargin)
%CREATE_SOLVER_FILE *Insert a one line summary here*
%   [] = create_solver_file(varargin)
%
% CREATE_SOLVER_FILE uses the U_PACKARGS interface function
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
% Created: 12-Sep-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, 0, ...
    'train_net', [],... 
    'test_net', [],... 
    'test_initialization', 'false',...
    'test_iter', '0',...
    'test_interval', '0',...
    'base_lr', '0.01',...
    'momentum', '0.9',...
    'weight_decay', '0.0005',...
    'lr_policy', '"inv"',...
    'gamma', '0.0001',...
    'power', '0.75',...
    'display', '100',...
    'max_iter', '10000',...
    'snapshot', '5000',...
    'snapshot_prefix', [],...
    'solver_mode', 'CPU');
clear varargin;

fid = fopen(solver_filename, 'wt');

param_names = fieldnames(args);

for i_p = 1:length(param_names)
    p_name = param_names{i_p};
    
    if ( (strcmpi(p_name, 'test_iter') || strcmpi(p_name, 'test_interval'))...
            && strcmpi(args.(p_name), '0') )
        fprintf(fid, '# %s: %s\n', p_name, args.(p_name));
    elseif ~isempty(args.(p_name))
        fprintf(fid, '%s: %s\n', p_name, args.(p_name));
    end
end
fclose(fid);

    
function [file_list] = dir_to_file_list(dir_in, pathstr)
%DIR_TO_FILE_LIST convert struc returned from DIR to a cell array of full
%filepath strings
%   [file_list] = dir_to_file_list(dir_in)
%
% Inputs:
%      dir_in - Either a string template, on which DIR will be called or a
%       file list struc as returned froma previous call to DIR
%
%      pathstr - optional, folder to which filenames will be appended (i.e.
%       to be used when a struc is given as input or to override the path
%       in the input string
%
%
% Outputs:
%      file_list - a cell array containing full filepaths
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 14-Aug-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if ischar(dir_in)
    if ~exist('pathstr', 'var');
        pathstr = [fileparts(dir_in) filesep];
    end
    dir_in = dir(dir_in);
else
    if ~isfield(dir_in, 'name')
        error('Input directory structure must contain name as a field');
    end
    if ~exist('pathstr', 'var');
        pathstr = [];
    end
end

nfiles = length(dir_in);
file_list = cell(nfiles, 1);

for i_file = 1:nfiles
    file_list{i_file} = [pathstr dir_in(i_file).name];
end

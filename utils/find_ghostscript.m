function [cmd] = find_ghostscript(default_path, default_exe)
%FIND_GHOSTSCRIPT Finds the latest version of ghostscript installed -
%Windows only
%   [cmd] = find_ghostscript(default_path)
%
% Inputs:
%      default_path - default program files path, if empty, defaults to
%      C:/Program Files (x86)/gs/
%
%
% Outputs:
%      cmd - fullpath to ghostscript executable, if empty, GhostScript
%      wasn't found
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 10-Feb-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('default_path', 'var') || isempty(default_path)
    default_path = 'C:\Program Files (x86)\gs\';
end
if ~exist('default_exe', 'var') || isempty(default_exe)
    default_exe = '\bin\gswin32c.exe';
end

% Find the full path to a ghostscript executable
cmd = '';

% Find all subdirs in the default location
dir_list = dir(default_path);
ver_num = 0;

for a = 1:numel(dir_list)
    % If there are multiple versions, use the newest
    ver_num2 = sscanf(dir_list(a).name, 'gs%g');
    if ~isempty(ver_num2) && ver_num2 > ver_num
        cmd2 = [default_path dir_list(a).name default_exe];
        if exist(cmd2, 'file') == 2
            cmd = ['"' cmd2 '"'];
            ver_num = ver_num2;
        end
    end
end

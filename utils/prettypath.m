function path_out = prettypath(path_in)
% Function to take a path that may contain backslashes and repeated file
% separators and return a pretty path with none of that nonsense and a
% terminating forward slash
%
% Should not be used with format specifiers (e.g. '\n', '\t')

% First, switch backslashes to forward slashes
path_out = strrep(path_in,'\','/');

% Repeatedly replace pairs of slashes with a single one
% Note: not good for URLs (e.g. http://www.manchester.ac.uk)
changed = true;
while changed
    old_path = path_out;
    path_out = strrep(path_out,'//','/');
    changed = ~strcmp(path_out, old_path);
end

% % Add terminating forward slash
% if (path_out(end) ~= '/')
%     path_out(end+1) = '/';
% end

function subdirs2path(current_dir) %SUBDIRS2PATH - check subdirs and add directory contents to path recursively
%
% SUBDIRS2PATH(current_dir) - check subdirs and add directory contents to path recursively
%
%       Inputs : current_dir - the directory to check
%
%       Outputs : none

% add this directory to the path
path(path,current_dir)

%Now check if any child directories to add
d2 = dir(current_dir);
d = {};
for i = 1:length(d2)
    if ~d2(i).isdir || strcmp(d2(i).name,'.') || strcmp(d2(i).name,'..') || strcmp(d2(i).name,'.svn') ...
            || strcmp(d2(i).name,'private') || strcmp(d2(i).name,'common') || strcmp(d2(i).name,'src') ...
            || strcmp(d2(i).name,'distrib')
        % do nothing
    else
        d = [d; {[current_dir '/' d2(i).name]}];
    end
end

for i = 1:length(d)
    subdirs2path(d{i})
end


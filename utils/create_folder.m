function path_name = create_folder(path_name)
% Create folder for path_name and return a tidied up version of the path
% name.

if ~isempty(path_name)
    % Tidy up the path first
    path_name = prettypath(path_name, true);
    
    % Create the folder if it doesn't exist already
    if ~exist(path_name, 'dir')
        mkdir(path_name);
		if ~ispc
			fileattrib(path_name,'+w','g'); % make folder writable...
			fileattrib([path_name,'..'],'+w','g'); % ...and its parent, too
		end
    end
end

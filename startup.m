%This is a script to be used by the server version of matlab
% to set my paths so I can use my own mfiles (adapted from Kola's).

%Turn off change notification warnings
system_dependent DirChangeHandleWarn Never;

%Set the path to the default path
%path(pathdef);

% this is get_username() except it isn't in the path yet
if ispc
	username_str = getenv('username');
elseif isunix
	[retval,username_str] = unix('echo $USER');
	username_str(end) = [];
end
	
switch username_str
	case {'mberks', 'Michael Berks'}
		%determine platform
		if ispc
			%my pc is configured as I want for now so just set window docking here    
			set(0,'DefaultFigureWindowStyle','docked')
			set(0,'DefaultFigureColormap', jet(256));
			
			home = 'C:\isbe\';
			path(path,strcat(home,'matlab_code\trunk'));
			subdirs2path(strcat(home,'matlab_code\trunk'));

			% add paths from isbe repository
			path(path,strcat(home,'matlab_code\isbe_repository\amb'));
			path(path,strcat(home,'matlab_code\isbe_repository\Geometry'));
			path(path,strcat(home,'matlab_code\isbe_repository\Models'));
			path(path,strcat(home,'matlab_code\isbe_repository\Stats'));
			path(path,strcat(home,'matlab_code\isbe_repository\Utils'));
		else
			%we're on the server so need to set paths (might change the
			%configuration of this later
			home = '/home/mberks/';
			path(path,strcat(home,'matlab_code/trunk'));
			subdirs2path(strcat(home,'matlab_code/trunk'));

			% add paths from isbe repository
			path(path,strcat(home,'matlab_code/amb'));
			path(path,strcat(home,'matlab_code/Geometry'));
			path(path,strcat(home,'matlab_code/Models'));
			path(path,strcat(home,'matlab_code/Stats'));
			path(path,strcat(home,'matlab_code/Utils'));
        end
        
    case 'momeemb2'
        if ispc
			%This should never trigger but hey...
		else
			%we're on the central uni cluster so need to set paths (might change the
			%configuration of this later
            if isdeployed
                display('Don''t do anything for deployed applications!');
            else
                path(path,strcat('matlab_code/trunk'));
                subdirs2path(strcat('matlab_code/trunk'));
            end
            
        end

	case 'ptresadern',
		%determine platform
		if ispc
			%my pc is configured as I want for now so just set window docking here    
			set(0,'DefaultFigureWindowStyle','docked')
			set(0,'DefaultFigureColormap', jet(256));
			
			home = 's:\projects\mammography\';
			
			% I have a startup.m in the c:\matlab folder to add paths as I need
			% them
		else
			%we're on the server so need to set paths (might change the
			%configuration of this later
			home = '/home/ptresadern/';
			path(path,strcat(home,'matlab/trunk'));
			subdirs2path(strcat(home,'matlab/trunk'));
			disp(path);
		end
end

clear

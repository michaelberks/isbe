function myroot = asymmetryroot(shared)
% get correct root variable for pc/cluster    

if ~exist('shared','var'), shared = false; end

if ischar(shared)
	switch lower(shared)
		case 'local', shared = false;
		case 'shared', shared = true;
	end
end

if ispc %work/home desktop
	switch get_username
		case {'mberks', 'Michael Berks'}
			if shared,	myroot = 'Z:\asym\';
			else		myroot = 'C:\isbe\asymmetry_project\';
			end
		case 'ptresadern'
			if shared,	myroot = 'A:\asym\';
			else		myroot = 'U:\projects\mammography\';
			end
		case 'Phil Tresadern' % phil's home computer
			if shared,	myroot = 'D:\Projects\mammography\';
			else		myroot = 'D:\Projects\mammography\';
			end
		otherwise
			myroot = 'E:\asymmetry_project\';
	end
else % hydra cluster
	switch get_username
		case {'mberks', 'ptresadern'}
			myroot = '/san/images/asym/asym/';
		case 'momeemb2'
			myroot = 'scratch/asym/';
		otherwise
			myroot = '/san/images/asym/asym/';
	end

	%myroot = '/home/mberks/asymmetry_project/';
end
	

function vessels = read_vessels_from(filename)

vessels = {};

f_debug = (nargin == 0 && nargout == 0);
if (f_debug)
	filename = 'U:\projects\nailfold\images\tmp.txt';
end

fid = fopen(filename);
while ~feof(fid)
	vessel_no = textscan(fid,'vessel %d');
	vessel_no = vessel_no{1};
	
	% check that file is in expected format
	if isempty(vessel_no)
		return
	end

	if ~isempty(vessel_no)
		ignore = fgetl(fid); % opening brace
		vessel_str = textscan(fid,'%[^}]');
		if ~isempty(vessel_str{1})
			vessel_pts = str2num(vessel_str{1}{1});
			vessels{vessel_no} = vessel_pts;
		end
		ignore = fgetl(fid); % closing brace
	end
end
fclose(fid);

if f_debug
	clear;
end

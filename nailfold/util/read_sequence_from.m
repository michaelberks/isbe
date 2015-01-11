function [sequence] = read_sequence_from(filename)
% Parse the *_markup.txt file output by the ncm_qmarkup.exe application
% into a structure whose fields have corresponding names.

if (nargin == 0 && nargout == 0), test_script(); return; end

sequence = [];

% Get cell structure of strings, removing any comments
fid = fopen(filename,'r');
    frewind(fid);
    s = textscan(fid, '%s', 'commentstyle', '//'); s = s{1};
fclose(fid);

% Turn it into a single string with no line breaks
s = sprintf('%s ', s{:});

% Find a token followed by a colon
[tok,pos] = textscan(s, '%s:'); tok = tok{1};
s = s(pos+1:end);

tok = tok{1}(1:end-1); % get just the field name
switch lower(strtrim(tok))
    case {'ncm_image_sequence'}
        [iOpen,iClosed] = find_braces(s);
        substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
        s = s(iClosed(1)+1:end);
        %try
            sequence = read_ncm_image_sequence(substr);
        %catch
            display(['File: ', filename]);
        %end
        
    otherwise
        error(['Unexpected token: ',tok]);
end

function [ncm_sequence] = read_ncm_image_sequence(s)
% Break down the ncm_annotation structure

ncm_sequence = struct(...
    'version',[],...
    'sequence_id',[],...
    'sequence_name', [],...
    'time_started',[],...
    'time_finished',[],...
    'session_id',[],...
    'hand',[],...
    'digit',[],...
    'num_frames',[],...
    'acceptable_quality',[],...
    'image_dir',[],...
    'preview_mosaic',[],...
    'full_mosaic',[],...
    'preview_count_map',[],...
    'full_count_map',[],...
    'ncm_camera_properties', [],...
    'ncm_motor_properties',[],...
    'frames', []);                  

while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'version','sequence_id','session_id', 'num_frames', 'digit', 'acceptable_quality'}
            % numerical value
            [value,pos] = textscan(s, '%d'); 
            s = s(pos+1:end);
            ncm_sequence.(tok) = value{1};

        case {'sequence_name', 'hand', 'preview_mosaic', 'full_mosaic', 'preview_count_map' 'full_count_map'}
            % string value
            % Observer name sometimes has a space in it. Grrrr.
            [str2,tok2] = strtok(s, ':');
            [toks,ignore] = textscan(str2, '%s');
            value = '';
            for i = 1:length(toks{1})-1
                value = [value ' ' toks{1}{i}];
            end
            s = [toks{1}{end}, tok2];
            ncm_sequence.(tok) = strtrim(value);
            
        case {'image_dir'}
            % string value, might have a ':\' or ':/' in it (e.g. drive),
            % might have spaces. Grrrr.
            [str2,tok2] = strtok(s, ':');
            
            if (tok2(2)=='\' || tok2(2) == '/')
                [str3,tok2] = strtok(tok2(2:end), ':');
                str2 = [str2 ':' str3];
            end
            
            [toks,ignore] = textscan(str2, '%s');
            
            value = '';
            for i = 1:length(toks{1})-1
                value = [value ' ' toks{1}{i}];
            end
            s = [toks{1}{end}, tok2];
            ncm_sequence.(tok) = strtrim(value);

        case {'time_started', 'time_finished'}
            % multiple strings
            [value,pos] = textscan(s, '%s', 3); 
            s = s(pos+1:end);
            ncm_sequence.(tok) = strtrim(sprintf('%s ',value{:}{:}));
         
        case {'ncm_camera_properties'}
            % camera properties structure
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_sequence.ncm_camera_properties = read_camera_properties(substr);
            
        case {'ncm_motor_properties'}
            % motor properties structure
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_sequence.ncm_motor_properties = read_motor_properties(substr);
            
        case {'frames'}
            % frames structure
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_sequence.frames = read_frames(substr);

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

function [frames] = read_frames(s)
% Break down the ncm_annotation structure

% Initialize with an empty ncm_vessel
frames = read_ncm_frame('');

n_frames = 0;
while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'ncm_frame'}
            % vessels structure
            n_frames = n_frames + 1;
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            frames(n_frames) = read_ncm_frame(substr);

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

if (n_frames == 0)
    frames = [];
end

function [ncm_camera_properties] = read_camera_properties(s)
% Break down the ncm_vessel structure

ncm_camera_properties = struct(...
    'camera_type',[],...
    'camera_id',[],...
    'frame_rate',[],...
    'exposure_time',[],...
    'gain',[],...
    'frame_height',[],...
    'frame_width',[]); 

while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'camera_type', 'camera_id'}
            % string values with spaces
            [str2,tok2] = strtok(s, ':');
            [toks,ignore] = textscan(str2, '%s');
            value = '';
            for i = 1:length(toks{1})-1
                value = [value ' ' toks{1}{i}];
            end
            s = [toks{1}{end}, tok2];
            ncm_camera_properties.(tok) = strtrim(value);
            
        case {'frame_rate', 'exposure_time', 'gain', 'frame_height', 'frame_width'}
            % numerical value
            [value,pos] = textscan(s, '%d'); 
            s = s(pos+1:end);
            ncm_camera_properties.(tok) = value{1};

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

function [ncm_motor_properties] = read_motor_properties(s)
% Break down the ncm_vessel structure
                    
ncm_motor_properties = struct(...
    'apt_x_id',[],...
    'apt_y_id',[],...
    'apt_z_id',[],...
    'x_reversed',[],...
    'y_reversed',[],...
    'z_reversed',[]);

while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'apt_x_id', 'apt_y_id', 'apt_z_id'}
            % string value
            [value,pos] = textscan(s, '%s', 1); 
            s = s(pos+1:end);
            ncm_motor_properties.(tok) = strtrim(value{1}{1});
        case {'x_reversed', 'y_reversed', 'z_reversed'}
            % numerical value
            [value,pos] = textscan(s, '%d'); 
            s = s(pos+1:end);
            ncm_motor_properties.(tok) = value{1};

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

function [ncm_frame] = read_ncm_frame(s)
% Break down the ncm_vessel structure

ncm_frame = struct(...
    'frame_name',[],...
    'frame_number',[],...
    'time',[],...
    'motor_x',[],...
    'motor_y',[],...
    'motor_z',[],...
    'sharpness',[],...
    'align_x',[],...
    'align_y',[],...
    'aligned_to',[],...
    'align_status',[]);

while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'frame_name'}
            % string value
            [value,pos] = textscan(s, '%s', 1); 
            s = s(pos+1:end);
            ncm_frame.(tok) = strtrim(value{1}{1});
        case {'frame_number', 'time', 'motor_x', 'motor_y', 'motor_z'...
                'sharpness', 'align_x', 'align_y','aligned_to', 'align_status'}
            % numerical value
            [value,pos] = textscan(s, '%f'); 
            s = s(pos+1:end);
            ncm_frame.(tok) = value{1};

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

function [ncm_vessel_properties] = read_ncm_vessel_properties(s)
% Break down the ncm_vessel_properties structure

ncm_vessel_properties = struct('version',[],...
                               'is_distal',[],...
                               'size',[],...
                               'shape',[]);

while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'version'}
            % numerical value
            [value,pos] = textscan(s, '%d'); 
            s = s(pos+1:end);
            ncm_vessel_properties.(tok) = value{1};

        case {'is_distal'}
            % binary string value ('yes'/'no')
            [value,pos] = textscan(s, '%s', 1); 
            s = s(pos+1:end);
            is_yes = strcmp(strtrim(value{1}{1}),'yes');
            ncm_vessel_properties.(tok) = is_yes;

        case {'size','shape'}
            % string value
            [value,pos] = textscan(s, '%s', 1); 
            s = s(pos+1:end);
            ncm_vessel_properties.(tok) = strtrim(value{1}{1});

        otherwise
            error(['Unexpected token: ',tok]);
    end
end



function [iOpen, iClosed] = find_braces(s)
% Find positions of opening and closing braces
iOpen = findstr(s,'{');
iClosed = findstr(s,'}');
if (numel(iOpen) ~= numel(iClosed))
    error('Unequal number of opening and closing braces');
end

% Match opening braces to closing ones
allInds = sort([iOpen iClosed]);
newOpen = zeros(size(iOpen));
newClosed = zeros(size(iClosed));
n = 1;
for i = length(allInds):-1:1
    % look for first opening where the next index is a close
    if ~isempty(find(iOpen == allInds(i),1)) && ...
       any(iClosed == allInds(i+1))
        % Found match
        newOpen(n) = allInds(i);
        newClosed(n) = allInds(i+1);
        allInds(i:i+1) = [];
        n = n + 1;
    end
end
[iOpen,inds] = sort(newOpen);
iClosed = newClosed(inds);


%% Test script
function test_script()
clc;

pathname = '';
filename = ['example_markup.txt'];

if ~exist(filename,'file')
    [filename,pathname] = uigetfile('*.txt');
end

if ~isempty(filename)
    markup = read_markup_from([pathname,filename]);
    disp(markup);
end

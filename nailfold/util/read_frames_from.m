function [markup] = read_frames_from(filename)
% Parse the *_markup.txt file output by the ncm_qmarkup.exe application
% into a structure whose fields have corresponding names.

if (nargin == 0 && nargout == 0), test_script(); return; end

markup = [];

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
    case {'ncm_annotation'}
        [iOpen,iClosed] = find_braces(s);
        substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
        s = s(iClosed(1)+1:end);
        try
            markup = read_ncm_annotation(substr);
        catch
            display(['File: ', filename]);
        end
        
    otherwise
        error(['Unexpected token: ',tok]);
end

function [ncm_annotation] = read_ncm_annotation(s)
% Break down the ncm_annotation structure

ncm_annotation = struct('version',[],...
                        'observer',[],...
                        'timestamp',[],...
                        'time_taken',[],...
                        'image_grade',[],...
                        'vessels',[],...
                        'haemorrhages',[]);

while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'version','time_taken'}
            % numerical value
            [value,pos] = textscan(s, '%d'); 
            s = s(pos+1:end);
            ncm_annotation.(tok) = value{1};

        case {'observer','image_grade'}
            % string value
            % Observer name sometimes has a space in it. Grrrr.
            [str2,tok2] = strtok(s, ':');
            [toks,ignore] = textscan(str2, '%s');
            value = '';
            for i = 1:length(toks{1})-1
                value = [value ' ' toks{1}{i}];
            end
            s = [toks{1}{end}, tok2];
            ncm_annotation.(tok) = strtrim(value);

        case {'timestamp'}
            % multiple strings
            [value,pos] = textscan(s, '%s', 5); 
            s = s(pos+1:end);
            ncm_annotation.timestamp = strtrim(sprintf('%s ',value{:}{:}));
            
        case {'vessels'}
            % vessels structure
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_annotation.vessels = read_vessels(substr);

        case {'haemorrhages'}
            % haemorrhages structure
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_annotation.haemorrhages = read_haemorrhages(substr);

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

function [vessels] = read_vessels(s)
% Break down the ncm_annotation structure

% Initialize with an empty ncm_vessel
vessels = read_ncm_vessel('');

n_vessels = 0;
while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'ncm_vessel'}
            % vessels structure
            n_vessels = n_vessels + 1;
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            vessels(n_vessels) = read_ncm_vessel(substr);

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

if (n_vessels == 0)
    vessels = [];
end


function [ncm_vessel] = read_ncm_vessel(s)
% Break down the ncm_vessel structure

ncm_vessel = struct('version',[],...
                    'ncm_vessel_properties',[],...
                    'anchor',[],...
                    'points',[],...
                    'apices',[]);

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
            ncm_vessel.(tok) = value{1};

        case {'ncm_vessel_properties'}
            % structure
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_vessel.(tok) = read_ncm_vessel_properties(substr);

        case {'anchor'}
            % 2D point
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_vessel.(tok) = read_2d_points(substr);

        case {'points'}
            % 3D points (x, y, width)
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_vessel.(tok) = read_3d_points(substr);

        case {'apices'}
            % apices structures
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_vessel.(tok) = read_apices(substr);

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

function [apices] = read_apices(s)
% Read a vector of apex structures

% Initialize with an empty ncm_apex
apices = read_ncm_apex('');
            
n_apices = 0;
while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'ncm_apex'}
            % apex structure
            n_apices = n_apices + 1;
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            apices(n_apices) = read_ncm_apex(substr);

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

function [ncm_apex] = read_ncm_apex(s)
% Read ncm_apex structure from string

ncm_apex = struct('version',[],...
                  'inner_point',[],...
                  'outer_point',[]);

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
            ncm_apex.(tok) = value{1};

        case {'inner_point','outer_point'}
            % 2D point
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_apex.(tok) = read_2d_points(substr);

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

function [haemorrhages] = read_haemorrhages(s)

% Initialize with an empty ncm_haemorrhage
haemorrhages = read_ncm_haemorrhage('');

n_haemorrhages = 0;
while ~isempty(s)
    % Find a token followed by a colon
    [tok,pos] = textscan(s, '%s:'); tok = tok{1};
    s = s(pos+1:end);

    tok = tok{1}(1:end-1); % get just the field name
    switch lower(strtrim(tok))
        case {'ncm_haemorrhage'}
            % vessels structure
            n_haemorrhages = n_haemorrhages + 1;
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            haemorrhages(n_haemorrhages) = read_ncm_haemorrhage(substr);

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

if (n_haemorrhages == 0)
    haemorrhages = [];
end


function [ncm_haemorrage] = read_ncm_haemorrhage(s)
% Break down the ncm_vessel structure

ncm_haemorrage = struct('version',[],...
                        'anchor',[]);

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
            ncm_haemorrage.(tok) = value{1};

        case {'anchor'}
            % 2D point
            [iOpen,iClosed] = find_braces(s);
            substr = strtrim(s(iOpen(1)+1:iClosed(1)-1));
            s = s(iClosed(1)+1:end);
            ncm_haemorrage.(tok) = read_2d_points(substr);

        otherwise
            error(['Unexpected token: ',tok]);
    end
end

function [points] = read_2d_points(s)
% convert list of (x,y) coordinates to a matrix of point locations
points = reshape(str2num(s),2,[])';

function [points] = read_3d_points(s)
% convert list of (x,y,w) coordinates to a matrix of point locations
points = reshape(str2num(s),3,[])';

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

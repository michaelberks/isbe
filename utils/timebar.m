function [tb_handle] = timebar(varargin)
%TIMEBAR Creates a waitbar that automatically updates you with the
%estimated time remaining
%   [tb_handle] = timebar(varargin)
%
% Inputs:
%      varargin - mostly key-value pairs to update 
%
%
% Outputs:
%      tb_handle - handle of the created/updated timebar
%
%
% Example:
%   tb = timebar('title','Processing','limit',1000); % create new timebar
%   timebar('closeall'); % close all open timebars
%   timebar(tb,'reset'); % reset timer
%   timebar(tb,'limit',100); % set max ticks to 100
%   timebar(tb,'title','Computing FFT'); % change title
%   timebar(tb,'advance'); % advance by one tick
%   timebar(tb,'advance',100); % advance by 100 ticks
%   timebar(tb,'position',[0 0 10 100]); % change position of window
%   timebar(tb,'close'); % close the timebar
%
% Notes:
%
% See also:
%
% Created: 10-Feb-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

% do nothing for linux
if ~ispc, tb_handle = 0; return; end

% get parameters from varargin
[tb_handle,f_reset,f_close,f_closeall,title,limit,n_ticks,unused] = ...
	parse_args(varargin{:});

if isempty(tb_handle)
    if (f_closeall)
        existing_tbs = findall(0,'Tag','TMWWaitbar');
        close(existing_tbs);
        return;
    end
    
	% if no handle given then create new timebar
	
	% define default position [x,y,w,h] (bottom right of screen)
	pos = [975 52 270 56.25];
	
	% move waitbar up above existing waitbars
    existing_tbs = findall(0,'Tag','TMWWaitbar');
	pos(2) = pos(2) + length(existing_tbs-1)*90; 
	
	% create waitbar at desired position
	tb_handle = waitbar(0,'Estimated time remaining = ??:??:??',...
												'Position',pos,...
												'Name','Timebar');

	% timebar.UserData stores:
	%   [limit current_progress start_time]
	
	% set some default values
	set(tb_handle,'UserData',[100 0 clock]); % def
else
	% if handle given then check it is valid
	
	if ~ishandle(tb_handle)
		warning('TIMEBAR:InvalidHandle','Invalid timebar handle');
		return;
	end
end

% set the title of the window
if ~isempty(title)
	set(tb_handle,'Name',title);
end

% if handle is a timebar then close it
if f_close && strcmp(get(tb_handle,'Tag'),'TMWWaitbar')
    close(tb_handle);
	tb_handle = 0;
end

% if reset flag set then restart clock and set progress to 0
if f_reset
	ud = get(tb_handle,'UserData');
	ud(2:end) = [0 clock];
	set(tb_handle,'UserData',ud);
	
	waitbar(0,tb_handle,'Estimated time remaining = ??:??:??');
end

% update limit of timebar
if ~isempty(limit)
	ud = get(tb_handle,'UserData');
	ud(1) = limit;
	set(tb_handle,'UserData',ud);
	
	% this assumes the limit is set early on and not reset halfway through
	% may need to account for this possibility
end

% update timebar progress and estimated time
if ~isempty(n_ticks)
	ud = get(tb_handle,'UserData');
	ud(2) = ud(2) + n_ticks; % number of ticks performed so far
	
	% estimate time remaining
	n_remaining = ud(1) - ud(2); % number of clicks left to do
	t = etime(clock,ud(3:end))/ud(2) * n_remaining;
	
	% convert to minutes and seconds
	hh = floor(t/3600); t = t-hh*3600;
	mm = floor(t/60); t = t-mm*60;
	ss = t;
	
	% update message displayed on timebar
	msg = sprintf('Estimated time remaining = %02.0f:%02.0f:%02.0f',hh,mm,ss);
	waitbar(ud(2)/ud(1),tb_handle,msg);
	set(tb_handle,'UserData',ud);
end

% pass any unused parameters as figure properties
if ~isempty(unused)
	while ~isempty(unused)
		set(tb_handle,unused{1},unused{2});
		unused(1:2) = [];
	end
end


function [handle,f_reset,f_close,f_closeall,title,limit,n_ticks,unused] = ...
    parse_args(varargin)
% get parameter values from varargin

% default values
handle = [];			% handle of the timebar
f_reset = false;        % reset flag
f_close = false;
f_closeall = false;
title = '';				% title of the window
limit = [];				% maximum number of ticks
n_ticks = [];			% number of ticks to advance by	
unused = {};			% unused parameters

% check for a handle as first parameter
if ~isempty(varargin) && isnumeric(varargin{1})
	handle = varargin{1};
	varargin(1) = [];
end

while ~isempty(varargin)
	switch lower(varargin{1})
		case 'reset',
			% reset counter
			f_reset = true;
			varargin(1) = [];
            
        case 'close',
			% close timebar
			f_close = true;
			varargin(1) = [];
			
        case 'closeall',
			% close all open timebars
			f_closeall = true;
			varargin(1) = [];

        case 'title',
			% get window title
			if length(varargin)<2
				error('Must give a value for timebar.title');
			end
			title = varargin{2};
			varargin(1:2) = [];
			
		case 'limit',
			% get maximum count
			if length(varargin)<2
				error('Must give a value for timebar.limit');
			end
			if ischar(varargin{2})
                limit = str2num(varargin{2});
            else	
                limit = varargin{2};
			end
			if isempty(limit)
				error(['Invalid value for timebar.limit: ',num2str(limit)]);
			end
			varargin(1:2) = [];
			
		case 'advance',
			% get number of ticks to advance by
			if length(varargin)<2
				% default advance of 1 tick
				n_ticks = 1;
				varargin(1) = [];
			else
				if ischar(varargin{2})
                    n_ticks = str2num(varargin{2});
                else
                    n_ticks = varargin{2};
				end
				if isempty(n_ticks)
					% use default advance of 1 tick
					n_ticks = 1;
					varargin(1) = [];
				else
					varargin(1:2) = [];
				end
			end
			
		otherwise,
			% unrecognised parameters get returned and passed to waitbar
			unused(end+1:end+2) = varargin(1:2);
			varargin(1:2) = [];
	end
end




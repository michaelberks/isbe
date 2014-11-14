function [] = display_images_shell(display_fun, num_images, curr_im, varargin)
%DISPLAY_IMAGES_SHELL *Insert a one line summary here*
%   [] = display_images_shell(display_fun, num_images)
%
% Inputs:
%      display_fun - *Insert description of input variable here*
%
%      num_images - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 11-Jan-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if nargin < 3
    curr_im = 1;
end

screen_size = get(0,'ScreenSize');
f1 = figure(...
    'windowstyle', 'normal',...
    'NumberTitle', 'off',...
    'Position', [0 30 screen_size(3), screen_size(4)-100],...
    'visible', 'off');

feval(display_fun, curr_im, f1, varargin);
set(f1,...
    'name', ['Running display function: @' display_fun ...
        '...      Image ' num2str(curr_im) ' of '  num2str(num_images)],...
    'KeyPressFcn', @key_press_Callback,...
    'visible', 'on');

% --------------------------------------------------------------------
    function key_press_Callback(hObject, eventdata) %#ok
    % Callback to...
    
    switch eventdata.Key;
        case 'leftarrow'
            if curr_im > 1
                h = waitbar(0,'Loading previous image. Please wait...');
                curr_im = curr_im - 1;
                feval(display_fun, curr_im, f1, varargin);
                set(f1, 'name', ['Running display function: ' display_fun ...
                    '. Image ' num2str(curr_im) ' of '  num2str(num_images)],...
                    'KeyPressFcn', @key_press_Callback);
                close(h);
            end
            
        case 'rightarrow'
            if curr_im < num_images
                h = waitbar(0,'Loading next image. Please wait...');
                curr_im = curr_im + 1;
                feval(display_fun, curr_im, f1, varargin);
                set(f1, 'name', ['Running display function: ' display_fun ...
                    '. Image ' num2str(curr_im) ' of '  num2str(num_images)],...
                    'KeyPressFcn', @key_press_Callback);
                close(h);
            end
    end

    end
end

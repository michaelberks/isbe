function [] = quick_apex_label(num_images, curr_im)
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
if ~exist('curr_im', 'var')
    curr_im = 1;
end

if ~exist('num_images', 'var')
    num_images = 100;
end

first_im = curr_im;
last_im = first_im + num_images - 1;

f1 = figure(...
    'windowstyle', 'normal',...
    'NumberTitle', 'off',...
    'Position', [100 100 1054, 800],...
    'visible', 'on',...
    'name', ['Running display function:' ...
        '...      Image ' num2str(curr_im) ' of '  num2str(num_images)],...
    'Colormap', gray(256),...
    'KeyPressFcn', @key_press_Callback);

a1 = axes(...
    'Parent', f1,...
    'Units', 'pixels',...
    'Position', [10 278 512 512],...
    'NextPlot', 'add',...
    'Ydir', 'reverse',...
    'DataAspectRatio', [1 1 1],...
    'PlotBoxAspectRatio', [1 1 1],...
    'XTick', [], 'YTick', []);

a2 = axes(...
    'Parent', f1,...
    'Units', 'pixels',...
    'Position', [532 278 512 512],...
    'NextPlot', 'add',...
    'Ydir', 'reverse',...
    'DataAspectRatio', [1 1 1],...
    'PlotBoxAspectRatio', [1 1 1],...
    'XTick', [], 'YTick', [],...
    'CLim', [0 1]);

uicontrol(...
    'Parent', f1,...
    'Style','pushbutton',...
    'String', 'Background',...
    'Tag','background',...
    'Callback', @negative_Callback,...
    'Position', [10, 10, 512, 258],...
    'Enable', 'on');

uicontrol(...
    'Parent', f1,...
    'Style','pushbutton',...
    'String', 'Capillary',...
    'Tag','capillary',...
    'Callback', @positive_Callback,...
    'Position', [532, 10, 512, 258],...
    'Enable', 'on');

candidate_path = 'C:\isbe\nailfold\data\wellcome_study\candidate_patches\candidate';
im_h1 = imagesc(rand(130),...
    'parent', a1);
im_h2 = imagesc(rand(130),...
    'parent', a2);
p1 = plot(a1, [1 130 130 1 1], [1 1 130 130 1], 'k', 'linewidth', 3);
p2 = plot(a2, [1 130 130 1 1], [1 1 130 130 1], 'k', 'linewidth', 3);
plot_function();

% --------------------------------------------------------------------
    function plot_function() 
        h = waitbar(0,'Loading image. Please wait...');
        candidate = load([candidate_path zerostr(curr_im,5) '.mat']);
        candidate.vessel_feature_patch(candidate.vessel_feature_patch < 80) = 80;
        set(im_h1, 'CData', candidate.vessel_feature_patch);
        set(im_h2, 'CData', candidate.vessel_prob_patch);
        set(f1, 'name', ['Candidate ' num2str(curr_im) ' (' num2str(first_im) ' to ' num2str(last_im) ')']);
                
        if isfield(candidate, 'label')
            if candidate.label
                set(p1, 'color', 'g');
                set(p2, 'color', 'g');
            else
                set(p1, 'color', 'r');
                set(p2, 'color', 'r');
            end
        else
            set(p1, 'color', 'y');
            set(p2, 'color', 'y');
        end
        close(h);
        
    end

% --------------------------------------------------------------------
    function forwards()
        if curr_im < last_im
            curr_im = curr_im + 1;
            plot_function();
        end      
    end

% --------------------------------------------------------------------
    function backwards()
        if curr_im > 1               
            curr_im = curr_im - 1;
            plot_function();               
        end     
    end

% --------------------------------------------------------------------
    function positive_Callback(hObject, eventdata) %#ok
        label = true; %#ok
        save([candidate_path zerostr(curr_im,5) '.mat'],...
            'label', '-append');
        forwards();
    end

% --------------------------------------------------------------------
    function negative_Callback(hObject, eventdata) %#ok
        label = false; %#ok
        save([candidate_path zerostr(curr_im,5) '.mat'],...
            'label', '-append');
        forwards();
    end

% --------------------------------------------------------------------
    function key_press_Callback(hObject, eventdata) %#ok
    % Callback to...
    
        switch eventdata.Key;
            case 'leftarrow'
                backwards()

            case 'rightarrow'
                forwards();
        end

    end
end

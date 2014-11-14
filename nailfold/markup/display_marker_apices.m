function [] = display_marker_apices()
%OBSERVER_STUDY_GUI *Insert a one line summary here*
%   [] = observer_study_gui()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 14-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Program code that runs
%
%--------------------------------------------------------------------------
%Must have window style set to normal
orig_window_style = get(0,'DefaultFigureWindowStyle');
if ~strcmp(orig_window_style, 'normal')
    display('Warning: changing window style to normal for function');
    set(0,'DefaultFigureWindowStyle','normal');
end

%Set constants/variables that persist for all user sessions:
color1 = [1 1 1];
color2 = [212 208 200]/255; %#ok
color3 = [0 0 0];
buff = 5;
screen_size = get(0,'ScreenSize');

%Create empty variables that exist globally and will be filled auxilliary
%functions - this set do not need to be reset for each session
if strcmpi(getenv('UserName'), 'mberks')
    mammo_dir = 'C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\';
    nailfold_dir = 'C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\';
else
    mammo_dir = 'E:\asymmetry_project\data\mammograms\2004_screening\abnormals\';
    nailfold_dir = 'E:\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\';
end

num_ims = [];
num_markers = [];
nailfold_strings = [];
marker_names = [];
im_num = [];
curr_im = [];

ui = [];
axes_pos = [];
panel_pos = [];

im1 = [];

meta_lx = [];
meta_rx = [];
meta_ly = [];
meta_ry = [];

%main program
create_main_fig
%get_pair_information;

%Reset window style
set(0,'DefaultFigureWindowStyle',orig_window_style);

%End of function

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Auxilliary functions so set up main UI figure (axes, buttons etc)
%
%--------------------------------------------------------------------------
    function create_main_fig
        %Generate main figure if it doens't already exist

        ui.main_fig = figure(...
            'Position', [0 30 screen_size(3), screen_size(4)-50],...
            'Visible','on',...
            'Name', 'Organisation maps display tool',...
            'NumberTitle', 'off',...
            'MenuBar', 'none',...
            'WindowStyle', 'normal',...
            'Color', color1,...
            'CloseRequestFcn', @quit_Callback);
        
        
        x_max = screen_size(3);
        y_max = screen_size(4)-50;
        button_h = 40;
        button_w = 100;
        text_w = 250;
        panel_h = 30;
        panel_w = x_max - 2*buff;
        axes_h = y_max - 3*buff - panel_h;
        axes_w = x_max - 2*buff;
        
        figure(ui.main_fig);
        set(ui.main_fig,...
            'Color', color3,...
            'CloseRequestFcn', @quit_Callback);

        
        axes_pos = [buff, panel_h+2*buff, axes_w, axes_h];
        panel_pos = [buff, buff, panel_w, panel_h];
        
        checkbox_w = panel_pos / num_markers;
        checkbox_h = panel_h;
        
        ui.frame = uicontrol(...
            'Style','frame',...
            'Parent', ui.main_fig,...
            'Units', 'pixels',...
            'Position', panel_pos,...
            'Visible', 'on');
        
        ui.panel = uipanel(...
            'Parent', ui.main_fig,...
            'Units', 'pixels',...
            'Position', panel_pos,...
            'Visible', 'on');
        
        ui.marker_checkboxes = zeros(num_markers,1);
        
        for i_m = 1:num_markers
            
            checkbox_pos = [(i_m-1)*checkbox_w 0 checkbox_w checkbox_h];
            
            ui.marker_checkboxes(i_m) = uicontrol(...
                'Style','checkbox',...
                'Parent', ui.panel,...
                'Units', 'pixels',...
                'Position', checkbox_pos,...
                'String', marker_names{i_m},...
                'Callback', toggle_marker_Callback,...
                'Visible', 'on');
        end
        
%         %------------------------------------------------------------------
%         ui.panel1 = uipanel(...
%             'Parent', ui.main_fig,...
%             'Units', 'pixels',...
%             'Position', [panel_pos(1)+buff y_max-buff-panel_h panel_w-2*buff panel_h],...
%             'Visible', 'on');
%            
%         ui.mammo_dir_text = uicontrol(... 
%             'Parent', ui.panel1,... 
%             'Style','text',...
%             'BackgroundColor', get(ui.panel, 'BackgroundColor'),...
%             'FontName', 'Arial',...
%             'String', 'Mammogram pairs folder:',...
%             'HorizontalAlignment', 'left',...
%             'Position', [0 button_h text_w 25]); 
% 
%         ui.mammo_dir_box = uicontrol(...
%             'Style', 'edit',...
%             'Position', [buff buff text_w button_h],...
%             'BackgroundColor', [1 1 1],...
%             'Parent', ui.panel1,...
%             'String', mammo_dir);
% 
%         ui.mammo_dir_select = uicontrol(... 
%             'Style', 'pushbutton',...
%             'Position', [2*buff+text_w buff button_w button_h],...
%             'String', 'Select',...
%             'Parent', ui.panel1,...
%             'Callback', @mammo_dir_select_Callback);
%         %---------------------------------------------------------------
%         %---------------------------------------------------------------
%         ui.panel8 = uipanel(...
%             'Parent', ui.main_fig,...
%             'Units', 'pixels',...
%             'Position', [panel_pos(1)+buff y_max-8*(buff+panel_h) panel_w-2*buff panel_h],...
%             'Visible', 'on');
%         
%         ui.min_m_slider_text = uicontrol(... 
%             'Parent', ui.panel8,... 
%             'Style','text',...
%             'BackgroundColor', get(ui.panel, 'BackgroundColor'),...
%             'FontName', 'Arial',...
%             'String', 'Select min/max values of mammogram contrast range: ' ,...
%             'HorizontalAlignment', 'left',...
%             'Position', [0 25 text_w+button_w+buff 25]); 
%         
%         ui.min_m_slider = uicontrol(... 
%             'Style', 'slider',...
%             'Position', [buff buff (text_w+button_w)/2 25],...
%             'String', 'Select',...
%             'Min', 0,...
%             'Max', 1,...
%             'Value', 0,...
%             'SliderStep', [0.01 0.1],...
%             'Parent', ui.panel8,...
%             'Enable', 'off',...
%             'Callback', @min_m_slider_Callback);
%         
%         ui.max_m_slider = uicontrol(... 
%             'Style', 'slider',...
%             'Position', [2*buff+(text_w+button_w)/2 buff (text_w+button_w)/2 25],...
%             'String', 'Select',...
%             'Min', 0,...
%             'Max', 1,...
%             'Value', 0,...
%             'SliderStep', [0.01 0.1],...
%             'Parent', ui.panel8,...
%             'Enable', 'off',...
%             'Callback', @max_m_slider_Callback);
%         %---------------------------------------------------------------
%         ui.panel9 = uipanel(...
%             'Parent', ui.main_fig,...
%             'Units', 'pixels',...
%             'Position', [panel_pos(1)+buff y_max-9*(buff+panel_h) panel_w-2*buff panel_h],...
%             'Visible', 'on');
%         
%         ui.zoom_on = uicontrol(... 
%             'Parent', ui.panel9,... 
%             'Style','togglebutton',...
%             'String','Zoom',...
%             'Tag','zoom_on',...
%             'Callback', @zoom_Callback,...
%             'Position', [buff, buff, button_w, button_h],...
%             'Enable', 'off'); 
%         
%         ui.pan_on = uicontrol(... 
%             'Parent', ui.panel9,... 
%             'Style','togglebutton',...
%             'String','Pan',...
%             'Tag','zoom_on',...
%             'Callback', @pan_Callback,...
%             'Position', [2*buff+button_w, buff, button_w, button_h],...
%             'Enable', 'off');        

        %---------------------------------------------------------------                    
        ui.axes = axes(...
            'Parent', ui.main_fig,...
            'Units', 'pixels',...
            'Visible', 'off');
        ui.region = imagesc([]);
        set(ui.axes,...
            'Position', axes_pos,...
            'Xtick', [],...
            'Ytick', [],...
            'YDir','reverse',...
            'NextPlot', 'add');        
        
        set(ui.main_fig,...
             'Colormap', gray(256));
         
        ui.distal_apices = zeros(num_markers,1);
        ui.distal_anchors = zeros(num_markers,1);
        ui.nondistal_anchors = zeros(num_markers,1);
        marker_colors = lines(num_markers);
        
        for i_m = 1:num_markers
        
            ui.distal_apices(i_m) = plot(1, 1,...
                'Parent', ui.axes,...
                'Visible', 'off',...
                'Marker', '^',...
                'MarkerEdgeColor', marker_colors(i_m,:),...
                'MarkerSize', 2,...
                'LineStyle', 'none');
            
            ui.distal_anchors(i_m) = plot(1, 1,...
                'Parent', ui.axes,...
                'Visible', 'off',...
                'Marker', 'o',...
                'MarkerEdgeColor', marker_colors(i_m,:),...
                'MarkerSize', 2,...
                'LineStyle', 'none');
            
            ui.nondistal_anchors(i_m) = plot(1, 1,...
                'Parent', ui.axes,...
                'Visible', 'off',...
                'Marker', '+',...
                'MarkerEdgeColor', marker_colors(i_m,:),...
                'MarkerSize', 2,...
                'LineStyle', 'none');
            
        end
        
    end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Auxilliary functions that control data
%
%--------------------------------------------------------------------------

    function get_nailfold_information
        nailfold_list = dir([nailfold_dir '\*.png']);
        num_ims = length(nailfold_list);
        nailfold_strings = cell(num_ims,1);
        
        
        if num_ims
            
            for i_n = 1:num_ims
                nailfold_strings{i_n} = nailfold_list(i_n).name;
            end
            
            curr_im = 1;
            
            %Update uicontrols now we have data
            set(ui.nailfold_number_selecter,...
                'Enable', 'on',...
                'String', nailfold_strings);

            %Load in the images for the first pair and upfate the pair
            %selecter
            update_curr_pair;
            
            set(ui.zoom_on, 'Enable', 'on');
            set(ui.pan_on, 'Enable', 'on');
            set(ui.meta_on, 'Enable', 'on');
            set(ui.capture_slider, 'Enable', 'on');
        else
            warndlg('No confocality maps found in this directory');
        end        
    end
%--------------------------------------------------------------------------       
    function load_mammograms
        h = waitbar(0,'Loading mammograms. Please wait...');
        
        %Work out the names of the associated mammograms
        mammo_name_l = dir([mammo_dir '\*' im_num 'L' pair_view '*']);
        mammo_name_r = dir([mammo_dir '\*' im_num 'R' pair_view '*']);
        mammo_name_l = [mammo_dir '\' mammo_name_l(1).name];
        mammo_name_r = [mammo_dir '\' mammo_name_r(1).name];
       
        %Load in the mammograms
        if strcmpi(mammo_name_r(end-3:end), '.mat')
            im1 = u_load(mammo_name_r);
        else
            im1 = imread(mammo_name_r);
        end
        if strcmpi(mammo_name_l(end-3:end), '.mat')
            im2 = u_load(mammo_name_l);
        else
            im2 = imread(mammo_name_l);
        end       
        
        %Now check for meta data (i.e. abnormality annotations)
        meta_lx = [];
        meta_rx = [];
        meta_ly = [];
        meta_ry = [];
        if exist([mammo_dir '\meta'], 'dir')
            meta_name_l = dir([mammo_dir '\meta\*' im_num 'L' pair_view '*']);
            if ~isempty(meta_name_l)
                meta_l = u_load([mammo_dir '\meta\' meta_name_l(1).name]);
                if ~isempty(meta_l)
                    meta_lx = meta_l(:,1);
                    meta_ly = meta_l(:,2);
                end
            end
            meta_name_r = dir([mammo_dir '\meta\*' im_num 'R' pair_view '*']);
            if ~isempty(meta_name_r)
                meta_r = u_load([mammo_dir '\meta\' meta_name_r(1).name]);
                if ~isempty(meta_r)
                    meta_rx = meta_r(:,1);
                    meta_ry = meta_r(:,2);
                end
            end    
        end
        close(h);
        
        update_mammo_axes_display;
    end
%--------------------------------------------------------------------------
    function load_maps
        
        h = waitbar(0,'Loading maps. Please wait...');

        try 
            map_name_r = [nailfold_dir '\' map_names_r(curr_map).name];
            
            try 
                s = load(map_name_r);
                im3 = s.out_var;
                scaling3 = s.scaling.m / 255;
                clear s;
            catch %#ok
                err = lasterror; %#ok
                display(err.message);
                
                im3 = u_load(map_name_r);
            end
            map_dist = num2str(str2double(map_name_r(end-6:end-4)));
        catch %#ok
            im3 = zeros(size(im1));
        end
        try
            map_name_l = [nailfold_dir '\' map_names_l(curr_map).name];
            try
                s = load(map_name_l);
                im4 = s.out_var;
                scaling4 = s.scaling.m / 255;
                clear s;
            catch %#ok
                err = lasterror; %#ok
                display(err.message);
                im4 = u_load(map_name_l);
            end
            map_dist = num2str(str2double(map_name_l(end-6:end-4)));
        catch %#ok
            im4 = zeros(size(im2));
        end
        
        set(ui.distance_slider_text, 'String', ['Select capture range weighting: sigma = ' map_dist]);
        update_map_axes_display;
        
        close(h);
    end
%--------------------------------------------------------------------------
    function update_map_axes_display
        
        %Smooth the maps by capture range sigma
        smooth3 = scaling3*imfilter(sum(im3,3), fspecial('gaussian', round(5*curr_capture), curr_capture));
        smooth4 = scaling4*imfilter(sum(im4,3), fspecial('gaussian', round(5*curr_capture), curr_capture));
        
        %Compute map limits for color scaling
        max_map = max([max(smooth3(:)) max(smooth4(:))]);
        map_lim = [0 2*max_map];
        
        set(ui.axes3, 'Clim', map_lim);
        set(ui.axes4, 'Clim', map_lim);
        set(ui.axes5,...
            'Clim', map_lim,...
            'XColor', 'red',...
            'YColor', 'red',...
            'Xtick', [],...
            'YTickLabel', num2str(linspace(0, max_map, length(get(ui.axes5, 'Ytick')))',3));
        
        %Make the maps visible
        set(ui.region3,...
            'Visible', 'on',...
            'CData', smooth3,...
            'XData', [1 axes_pos3(3)],...
            'YData', [1 axes_pos3(4)],...
            'ButtonDownFcn', @display_hist);
        set(ui.region4,...
            'Visible', 'on',...
            'CData', smooth4,...
            'XData', [1 axes_pos4(3)],...
            'YData', [1 axes_pos4(4)],...
            'ButtonDownFcn', @display_hist);
        set(ui.region5,...
            'Visible', 'on',...
            'CData', linspace(0, max_map, 128)',...
            'XData', [1 axes_pos5(3)],...
            'YData', [1 axes_pos5(4)]);
        
        set(ui.meta3,...
            'XData', axes_pos3(3)*meta_rx,...
            'YData', axes_pos3(4)*meta_ry);
        set(ui.meta4,...
            'XData', axes_pos4(3)*meta_lx,...
            'YData', axes_pos4(4)*meta_ly);
        
        set(ui.min_c_slider,...
            'Min', 0,...
            'Max', 0.99*max_map,...
            'Value', 0,...
            'Enable', 'on')
        set(ui.max_c_slider,...
            'Min', 0.01,...
            'Max', max_map,...
            'Value', max_map,...
            'Enable', 'on')
        
        set(ui.min_c_slider_text,...
            'String', ['Select min/max values of map contrast range: [0, ' num2str(max_map,3) ']']);
        
        
    end

    function update_mammo_axes_display
        
        %Compute the apsect ratios for these images (they may vary from
        %pair to pair)
        aspect_ratio_r = size(im1,2) / size(im1,1);
        aspect_ratio_l = size(im2,2) / size(im2,1);
        
        %Update the size and position of the axes
        axes_pos(1) = axes_pos(1) + axes_pos(3) - axes_pos(4)*aspect_ratio_r;
        axes_pos3(1) = axes_pos3(1) + axes_pos3(3) - axes_pos3(4)*aspect_ratio_r;
        axes_pos5(1) = axes_pos3(1) - 25;
        axes_pos6(1) = axes_pos(1) - 25;
        
        axes_pos(3) = axes_pos(4)*aspect_ratio_r;
        axes_pos2(3) = axes_pos2(4)*aspect_ratio_l;
        axes_pos3(3) = axes_pos3(4)*aspect_ratio_r;
        axes_pos4(3) = axes_pos4(4)*aspect_ratio_l;
        
        
        set(ui.axes,...
            'Position', axes_pos,...
            'Xlim', [0.5 axes_pos(3)+0.5],...
            'Ylim', [0.5 axes_pos(4)+0.5]);
        set(ui.axes2,...
            'Position', axes_pos2,...
            'Xlim', [0.5 axes_pos2(3)+0.5],...
            'Ylim', [0.5 axes_pos2(4)+0.5]);
        set(ui.axes3,...
            'Position', axes_pos3,...
            'Xlim', [0.5 axes_pos3(3)+0.5],...
            'Ylim', [0.5 axes_pos3(4)+0.5]);
        set(ui.axes4,...
            'Position', axes_pos4,...
            'Xlim', [0.5 axes_pos4(3)+0.5],...
            'Ylim', [0.5 axes_pos4(4)+0.5]);
        set(ui.axes5,...
            'Position', axes_pos5,...
            'Xlim', [0.5 axes_pos5(3)+0.5],...
            'Ylim', [0.5 axes_pos5(4)+0.5]);
        set(ui.axes6,...
            'Position', axes_pos6,...
            'Xlim', [0.5 axes_pos6(3)+0.5],...
            'Ylim', [0.5 axes_pos6(4)+0.5]);
       
        %Compute mammo limits for color scaling
        max_mammo = double(max([max(im1(:)) max(im2(:))]));
        min_mammo = double(min([min(im1(:)) min(im2(:))]));
        mammo_lim = [2*min_mammo - max_mammo max_mammo];
        
        set(ui.axes, 'Clim', mammo_lim);
        set(ui.axes2, 'Clim', mammo_lim);
        set(ui.axes6,...
            'Clim', mammo_lim,...
            'XColor', 'red',...
            'YColor', 'red',...
            'Xtick', [],...
            'YTickLabel', num2str(linspace(min_mammo, max_mammo, length(get(ui.axes6, 'Ytick')))',3));
        
        %Make the mammograms visible
        set(ui.region,...
            'Visible', 'on',...
            'CData', im1,...
            'XData', [1 axes_pos(3)],...
            'YData', [1 axes_pos(4)]);
        set(ui.region2,...
            'Visible', 'on',...
            'CData', im2,...
            'XData', [1 axes_pos2(3)],...
            'YData', [1 axes_pos2(4)]);
        set(ui.region6,...
            'Visible', 'on',...
            'CData', linspace(min_mammo, max_mammo, 128)',...
            'XData', [1 axes_pos6(3)],...
            'YData', [1 axes_pos6(4)]);
        
        set(ui.meta1,...
            'XData', axes_pos(3)*meta_rx,...
            'YData', axes_pos(4)*meta_ry);
        set(ui.meta2,...
            'XData', axes_pos2(3)*meta_lx,...
            'YData', axes_pos2(4)*meta_ly);
        
        set(ui.min_m_slider,...
            'Min', 0,...
            'Max', 0.99*max_mammo,...
            'Value', 0,...
            'Enable', 'on')
        set(ui.max_m_slider,...
            'Min', 0.01,...
            'Max', max_mammo,...
            'Value', max_mammo,...
            'Enable', 'on')
        
        set(ui.min_m_slider_text,...
            'String', ['Select min/max values of mammo contrast range: [' num2str(min_mammo,3) ', ' num2str(max_mammo,3) ']']);
        
        
    end
%----------------------------------------------------------------------
    function update_nailfold
        
        set(ui.nailfold_number_text, 'String', ['Select mammogram pair:' num2str(curr_im) ' of ' num2str(num_ims)]);
        set(ui.nailfold_number_selecter, 'Value', curr_im);
        
        if curr_im == num_ims
            set(ui.next_pair, 'Enable', 'off');
        else
            set(ui.next_pair, 'Enable', 'on');
        end
        if curr_im == 1
            set(ui.previous_pair, 'Enable', 'off');
        else
            set(ui.previous_pair, 'Enable', 'on');
        end


        %Load in the nailfold
        load_mammograms;
        load_markers;
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% UI Callbacks
%
%--------------------------------------------------------------------------
% --------------------------------------------------------------------
    function mammo_dir_select_Callback(hObject, eventdata) %#ok
    % Callback to...
        set(ui.mammo_dir_select, 'Enable', 'off');
        temp_dir = ...
            uigetdir(mammo_dir, 'Select the directory containing the mammograms');
        if temp_dir
            mammo_dir = temp_dir;
            set(ui.mammo_dir_box, 'string', mammo_dir);
        end
        set(ui.mammo_dir_select, 'Enable', 'on');
    end
%--------------------------------------------------------------------------
    function map_dir_select_Callback(hObject, eventdata) %#ok
    % Callback to...
        set(ui.map_dir_select, 'Enable', 'off');
        temp_dir = ...
            uigetdir(nailfold_dir, 'Select the directory containing the organisation maps');
        if temp_dir
            nailfold_dir = temp_dir;
            set(ui.map_dir_box, 'string', nailfold_dir);
        end
        set(ui.map_dir_select, 'Enable', 'on');
    end
%--------------------------------------------------------------------------
    function update_Callback(hObject, eventdata) %#ok
    % Callback to...
        get_pair_information;
    end    
% --------------------------------------------------------------------
    function quit_Callback(hObject, eventdata) %#ok
    % Callback executed if the user tries to quit
        
        %Check if they're ok to quit
        delete(ui.main_fig);
    end

% --------------------------------------------------------------------
    function distance_slider_Callback(hObject, eventdata) %#ok
    % Callback...
        
    	curr_map = round(get(ui.distance_slider, 'value'));
        load_maps;

    end

% --------------------------------------------------------------------
    function capture_slider_Callback(hObject, eventdata) %#ok
    % Callback...
        curr_capture = 2^(get(ui.capture_slider, 'value')-1);
        set(ui.capture_slider_text, 'String', ['Select capture range weighting: sigma = ' num2str(curr_capture,1)]);
        update_map_axes_display;
    end

% --------------------------------------------------------------------
    function previous_pair_Callback(hObject, eventdata) %#ok
    % Callback to the "zoom" button allowing the user to zoom in/out the 
    % mammogram
        curr_im = curr_im - 1;
        update_curr_pair;
    end

% --------------------------------------------------------------------
    function next_pair_Callback(hObject, eventdata) %#ok
    % Callback to the "pan" button allowing the user to pan around the 
    % mammogram
        curr_im = curr_im + 1;
        update_curr_pair;
    end
% --------------------------------------------------------------------
    function pair_number_selecter_Callback(hObject, eventdata) %#ok
    % Callback...
        curr_im = get(ui.nailfold_number_selecter, 'value');
        update_curr_pair;
    end
% --------------------------------------------------------------------
    function min_c_slider_Callback(hObject, eventdata) %#ok
    % Callback...
        min_slider = get(ui.min_c_slider, 'value');
        max_slider = max(min_slider+0.01, get(ui.max_c_slider, 'value'));
        
        set(ui.max_c_slider, 'value', max_slider);
        
        new_max_c = 2*max_slider - min_slider;
        new_min_c = min_slider;
        
        %Compute map limits for color scaling
        map_lim = [new_min_c new_max_c];
        
        %Compute mammo limits for color scaling
        set(ui.axes3, 'Clim', map_lim);
        set(ui.axes4, 'Clim', map_lim);
        
        set(ui.axes5,...
            'YTickLabel', num2str(linspace(min_slider, max_slider, length(get(ui.axes5, 'Ytick')))',3));
        %
        set(ui.min_c_slider_text,...
            'String', ['Select min/max value of map contrast range: [' num2str(min_slider,3) ', ' num2str(max_slider,3) ']']);
    end
% --------------------------------------------------------------------
    function max_c_slider_Callback(hObject, eventdata) %#ok
    % Callback...
        max_slider = get(ui.max_c_slider, 'value');
        min_slider = min(max_slider-0.01, get(ui.min_c_slider, 'value'));
        
        set(ui.min_c_slider, 'value', min_slider);
        
        new_max_c = 2*max_slider - min_slider;
        new_min_c = min_slider;
        
        %Compute map limits for color scaling
        map_lim = [new_min_c new_max_c];
        
        %Compute mammo limits for color scaling
        set(ui.axes3, 'Clim', map_lim);
        set(ui.axes4, 'Clim', map_lim);
        
        set(ui.axes5,...
            'YTickLabel', num2str(linspace(min_slider, max_slider, length(get(ui.axes5, 'Ytick')))',3));
        %
        set(ui.min_c_slider_text,...
            'String', ['Select min/max value of map contrast range: [' num2str(min_slider,3) ', ' num2str(max_slider,3) ']']);
    end
% --------------------------------------------------------------------
    function min_m_slider_Callback(hObject, eventdata) %#ok
    % Callback...
        min_slider = get(ui.min_m_slider, 'value');
        max_slider = max(min_slider+0.01, get(ui.max_m_slider, 'value'));
        
        set(ui.max_m_slider, 'value', max_slider);
        
        new_max_m = max_slider;
        new_min_m = 2*min_slider - max_slider;
        
        %Compute mammo limits for color scaling
        mammo_lim = [new_min_m new_max_m];
        
        %Compute mammo limits for color scaling
        set(ui.axes, 'Clim', mammo_lim);
        set(ui.axes2, 'Clim', mammo_lim);
        
        set(ui.axes6,...
            'YTickLabel', num2str(linspace(min_slider, max_slider, length(get(ui.axes6, 'Ytick')))',3));
        
        %
        set(ui.min_m_slider_text,...
            'String', ['Select min/max value of mammogram contrast range: [' num2str(min_slider,3) ', ' num2str(max_slider,3) ']']);
    end
% --------------------------------------------------------------------
    function max_m_slider_Callback(hObject, eventdata) %#ok
    % Callback...
        max_slider = get(ui.max_m_slider, 'value');
        min_slider = min(max_slider-0.01, get(ui.min_m_slider, 'value'));
        
        set(ui.min_m_slider, 'value', min_slider);
        
        new_max_m = max_slider;
        new_min_m = 2*min_slider - max_slider;
        
        %Compute mammo limits for color scaling
        mammo_lim = [new_min_m new_max_m];
        
        %Compute mammo limits for color scaling
        set(ui.axes, 'Clim', mammo_lim);
        set(ui.axes2, 'Clim', mammo_lim);
        
        set(ui.axes6,...
            'YTickLabel', num2str(linspace(min_slider, max_slider, length(get(ui.axes6, 'Ytick')))',3));
        
        %
        set(ui.min_m_slider_text,...
            'String', ['Select min/max value of mammogram contrast range: [' num2str(min_slider,3) ', ' num2str(max_slider,3) ']']);
    end
% --------------------------------------------------------------------
    function zoom_Callback(hObject, eventdata) %#ok
    % Callback to the "Zoom" button on the main control panel, allowing the
    % user to zoom in/out of any of 3 ROI figures (zooming occurs
    % simultaneously in all 3 figures)
        if get(ui.zoom_on, 'Value')
            set(ui.pan_on, 'Value', 0);
            axes(ui.axes); zoom on;
            axes(ui.axes2); zoom on;
            axes(ui.axes3); zoom on;
            axes(ui.axes4); zoom on;
          
        else
            axes(ui.axes); zoom off;
            axes(ui.axes2); zoom off;
            axes(ui.axes3); zoom off;
            axes(ui.axes4); zoom off;
            
            set(ui.region3, 'buttondownfcn', @display_hist);
            set(ui.region4, 'buttondownfcn', @display_hist);
        end
    end

    % --------------------------------------------------------------------
    function pan_Callback(hObject, eventdata) %#ok
    % Callback to the "Pan" button on the main control panel, allowing the
    % user to pan around any of 3 ROI figures (panning occurs
    % simultaneously in all 3 figures)
        if get(ui.pan_on, 'Value')
            set(ui.zoom_on, 'Value', 0);
            axes(ui.axes); pan on;
            axes(ui.axes2); pan on;
            axes(ui.axes3); pan on;
            axes(ui.axes4); pan on;
          
        else
            axes(ui.axes); pan off;
            axes(ui.axes2); pan off;
            axes(ui.axes3); pan off;
            axes(ui.axes4); pan off;
            
            set(ui.region3, 'buttondownfcn', @display_hist);
            set(ui.region4, 'buttondownfcn', @display_hist);
        end
    end
% --------------------------------------------------------------------
    function meta_Callback(hObject, eventdata) %#ok
    % Callback to the "Pan" button on the main control panel, allowing the
    % user to pan around any of 3 ROI figures (panning occurs
    % simultaneously in all 3 figures)
        if get(ui.meta_on, 'Value')
            set(ui.meta1, 'visible', 'on');
            set(ui.meta2, 'visible', 'on');
            set(ui.meta3, 'visible', 'on');
            set(ui.meta4, 'visible', 'on');
          
        else
            set(ui.meta1, 'visible', 'off');
            set(ui.meta2, 'visible', 'off');
            set(ui.meta3, 'visible', 'off');
            set(ui.meta4, 'visible', 'off');
        end
    end

% --------------------------------------------------------------------
    function display_hist(hObject, eventdata) %#ok
    %
        g_width = round(curr_capture);
        g_filt = fspecial('gaussian', 6*g_width + 1, curr_capture);
        phi = linspace(0, pi, 100);
        max_val = 0.1*max(get(ui.axes5, 'Clim'));
        
        if gca == ui.axes3
            XYZ = get(ui.axes3, 'CurrentPoint');
            x_lim = get(ui.axes3, 'Xlim');
            y_lim = get(ui.axes3, 'Ylim');
            x = round( size(im3, 2) * XYZ(1,1) / x_lim(2));
            y = round( size(im3, 1) * XYZ(1,2) / y_lim(2) );
            
            num_angles = size(im3,3);
            roi = scaling3*double(reshape(im3(y-3*g_width:y+3*g_width, x-3*g_width:x+3*g_width,:), [], num_angles));
            
        elseif gca == ui.axes4
            XYZ = get(ui.axes4, 'CurrentPoint');
            x_lim = get(ui.axes4, 'Xlim');
            y_lim = get(ui.axes4, 'Ylim');
            x = round( size(im4, 2) * XYZ(1,1) / x_lim(2) );
            y = round( size(im4, 1) * XYZ(1,2) / y_lim(2) );
            
            num_angles = size(im4,3);
            roi = scaling4*double(reshape(im4(y-3*g_width:y+3*g_width, x-3*g_width:x+3*g_width,:), [], num_angles));
           
            
        else
            display('Unknown axis caller');
        end

        set(ui.popup_fig, 'visible', 'on');
        axes(ui.popup_axes(ui.popup_idx));
        hold off;
        ui.popup_idx = mod(ui.popup_idx,2)+1;
        
        rad_vals = g_filt(:)' * roi;
        [theta rho] = polar_bar([rad_vals 0 0], pi*(-1:num_angles)/num_angles);

        %Plot radial axes marker
        plot([-max_val max_val], [0 0], 'k--'); hold on;
        plot([0 0], [0 max_val], 'k--');
        for ii = 2:2:10
            plot(ii*max_val*cos(phi)/10, ii*max_val*sin(phi)/10, 'k--');
        end

        plot(rho.*cos(theta), rho.*sin(theta)); axis equal;
        axis([-max_val max_val min(rho.*sin(theta)) max_val]);
        set(gca, 'xticklabel', []);
        title(['Orientation histogram for x = ' num2str(x), ' y = ' num2str(y)]);
        xlabel(['Sum = ' num2str(sum(rad_vals),3), ', S.D. = ' num2str(std(rad_vals),3) ', Entropy = ' num2str(mb_entropy(rad_vals),3)]);  
    end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---------------------- END OF FUNCTION -----------------------------------
%--------------------------------------------------------------------------
end
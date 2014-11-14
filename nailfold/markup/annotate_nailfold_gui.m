function annotate_nailfold_gui(varargin)

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
color2 = [212 208 200]/255; %#ok
color3 = [0 0 0];
buff = 5;
screen_size = get(0,'ScreenSize');
sigma = 3;

%Create empty variables that exist globally and will be filled auxilliary
%functions - this set do not need to be reset for each session
    
if strcmpi(getenv('UserName'), 'mberks') || strcmpi(getenv('UserName'), 'Michael Berks')
    nailfold_dir = 'C:\isbe\nailfold\images\';
    anno_dir = 'C:\isbe\nailfold\images\';
    
else
    try
        fid = fopen('V:\Desktop\mb_nailfold_software\software\default_locations.txt');
        def_txt = textscan(fid, '%s %s', 'delimiter', '=');
        fclose(fid);

        nailfold_dir = ...
            def_txt{2}{strncmpi(def_txt{1}, 'nailfold_dir', length('nailfold_dir'))};
        anno_dir = ...
            def_txt{2}{strncmpi(def_txt{1}, 'annotations_dir', length('annotations_dir'))};

    catch
        display('Problem opening default file locations');
        nailfold_dir = '';
        anno_dir = '';
    end
end

curr_vessel = 0;
curr_idx = 0;

ui = [];
axes_pos = [];

colorbar_pos = [];
panel_pos = [];

nailfold = [];
vessels = [];
x_vessels = [];
y_vessels = [];
scaling = 0;
axes_offset = [];
savename = [];

%main program
create_main_fig


%Reset window style
set(0,'DefaultFigureWindowStyle',orig_window_style);

%End of function

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Auxilliary functions so set up main UI figure (axes, buttons etc)
%
%--------------------------------------------------------------------------
    function create_main_fig
        %-----------------------------------------------------------------
        %
        %-----------------------------------------------------------------
        
        %get screen size of this computer
        x_max = screen_size(3);
        y_max = screen_size(4)-50;
        
        %Generate main figure
        ui.main_fig = figure(...
            'Position', [0 30 x_max, y_max-25],...
            'Visible','on',...
            'Name', 'Nailfod annotation tool',...
            'NumberTitle', 'off',...
            'MenuBar', 'none',...
            'WindowStyle', 'normal',...
            'Color', color3,...
            'CloseRequestFcn', @quit_Callback);  
        
        % Make the menu bar to load/save in files
        %------------------------------------------------------------------
        
        ui.file_menu = uimenu('Label','File');
        
        ui.file_open = uimenu(ui.file_menu,'Label','Open',...
            'Callback', @file_open_Callback);
                
        ui.save = uimenu(ui.file_menu,'Label','Save',...
            'Callback', @file_save_Callback,...
            'Enable', 'off');
                            
        ui.load = uimenu(ui.file_menu,'Label','Load vessels',...
            'Callback', @load_Callback,...
            'Enable', 'off');                    
                    
        ui.close = uimenu(ui.file_menu,'Label', 'Close file',...
            'Callback', @close_Callback,...
            'Enable', 'off');
        ui.quit = uimenu(ui.file_menu,'Label', 'Quit',...
            'Callback', @quit_Callback,...
            'Separator','on');                    
        
        % Make the panel of buttons and controls
        %------------------------------------------------------------------              
       
        ui.panel = uipanel('Units', 'pixels');
        
        %Define button size and text width
        button_h = 40;
        button_w = 100;
        
        %------------------------------------------------------------------
        %position the buttons
        %------------------------------------------------------------------
        % Row 1, counting upwards from bottom
        %------------------------------------
        rr = 1; cc = 1;
        %--------------------
        ui.min_m_slider_text = uicontrol(... 
            'Parent', ui.panel,... 
            'Style','text',...
            'BackgroundColor', get(ui.panel, 'BackgroundColor'),...
            'FontName', 'Arial',...
            'String', 'Select min/max values of nailfold contrast range: ' ,...
            'HorizontalAlignment', 'left',...
            'Position', [cc*buff+(cc-1)*button_w rr*buff+(rr-1)*button_h+18 4*button_w 20]); 
        
        ui.min_m_slider = uicontrol(... 
            'Style', 'slider',...
            'Position', [cc*buff+(cc-1)*button_w rr*buff+(rr-1)*button_h 2*button_w 25],...
            'String', 'Select',...
            'Min', 0,...
            'Max', 1,...
            'Value', 0,...
            'SliderStep', [0.01 0.1],...
            'Parent', ui.panel,...
            'Enable', 'off',...
            'Callback', @min_m_slider_Callback);
        cc = cc + 2;
        %--------------------
        ui.max_m_slider = uicontrol(... 
            'Style', 'slider',...
            'Position', [cc*buff+(cc-1)*button_w rr*buff+(rr-1)*button_h 2*button_w 25],...
            'String', 'Select',...
            'Min', 0,...
            'Max', 1,...
            'Value', 0,...
            'SliderStep', [0.01 0.1],...
            'Parent', ui.panel,...
            'Enable', 'off',...
            'Callback', @max_m_slider_Callback);
        cc = cc + 2;
        %--------------------
        %--------------------
        ui.auto = uicontrol(...
            'Parent', ui.panel,...
            'Style','pushbutton',...
            'String', 'Turn auto mode off',...
            'Tag','show_all',...
            'Callback', @auto_Callback,...
            'Position', [cc*buff+(cc-1)*button_w rr*buff+(rr-1)*button_h button_w button_h],...
            'Enable', 'off');
        %cc = cc + 1;
        %------------------------------------------------------------------
        rr = rr + 1; cc = 1;
        ui.zoom_on = uicontrol(...
            'Parent', ui.panel,...
            'Style','togglebutton',...
            'String','Zoom',...
            'Tag','zoom_on',...
            'Callback', @zoom_Callback,...
            'Position', [cc*buff+(cc-1)*button_w rr*buff+(rr-1)*button_h button_w button_h],...
            'Enable', 'off');
        cc = cc + 1;
        %------------------
        ui.pan_on = uicontrol(...
            'Parent', ui.panel,...
            'Style','togglebutton',...
            'String','Pan',...
            'Tag','pan_on',...
            'Callback', @pan_Callback,...
            'Position', [cc*buff+(cc-1)*button_w rr*buff+(rr-1)*button_h button_w button_h],...
            'Enable', 'off');
        cc = cc + 1;
        %--------------------
        ui.show_vessel = uicontrol(...
            'Parent', ui.panel,...
            'Style','pushbutton',...
            'String','Hide vessel',...
            'Tag','show_outline',...
            'Callback', @show_vessel_Callback,...
            'Position', [cc*buff+(cc-1)*button_w rr*buff+(rr-1)*button_h button_w button_h],...
            'Enable', 'off');
        cc = cc + 1;
        %--------------------
        ui.show_all = uicontrol(...
            'Parent', ui.panel,...
            'Style','pushbutton',...
            'String','Hide all vessels',...
            'Tag','show_all',...
            'Callback', @show_all_Callback,...
            'Position', [cc*buff+(cc-1)*button_w rr*buff+(rr-1)*button_h button_w button_h],...
            'Enable', 'off');
        cc = cc + 1;
        %--------------------
        ui.show_lines = uicontrol(...
            'Parent', ui.panel,...
            'Style','pushbutton',...
            'String','Hide lines',...
            'Tag','show_all',...
            'Callback', @show_lines_Callback,...
            'Position', [cc*buff+(cc-1)*button_w rr*buff+(rr-1)*button_h button_w button_h],...
            'Enable', 'off');
        %------------------------
        % Now finalise the panel position
        panel_h = (rr+1)*buff + rr*button_h;
        panel_w = x_max - 2*buff;
        panel_pos = [buff, buff, panel_w, panel_h];
        
        set(ui.panel, 'Position', panel_pos);
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        
        % Workout what size axes we can have given the panel size
        axes_w = x_max - 3*buff - 20;
        axes_h = y_max - 8*buff - panel_h;
        
        axes_pos = [2*buff+20  2*buff+panel_h axes_w axes_h];
        colorbar_pos = [buff axes_pos(2) 20 axes_h];
        
        % Create the axes to store the loaded nailfold and annotations
        ui.axes = axes(...
            'Parent', ui.main_fig,...
            'Units', 'pixels',...
            'Visible', 'off',...
            'Drawmode', 'fast');
        ui.nailfold = imagesc([]);
        set(ui.axes,...
            'Position', axes_pos,...
            'Xtick', [],...
            'Ytick', [],...
            'YDir','reverse',...
            'NextPlot', 'add');
        
        ui.colorbar_axes = axes(...
            'Parent', ui.main_fig,...
            'Units', 'pixels');
        ui.colorbar = imagesc([]);
        set(ui.colorbar_axes,...
            'Xtick', [],...
            'Yticklabel', [],...
            'YAxisLocation', 'right',...
            'Position', colorbar_pos); 

        set(ui.main_fig,...
            'Colormap', gray(256));
        ui.show_vessel_flag = 1;
        ui.show_all_flag = 1;
        ui.show_lines_flag = 1;
        ui.auto_flag = 1;
        ui.not_saved = 0;
        
    end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Auxilliary functions that control data
%
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --------------------------------------------------------------------
    function file_open_Callback(hObject, eventdata) %#ok %#ok
    % hObject    handle to menu_file_open_orig (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
        
        if ui.not_saved
            selection = questdlg('The vessel annotations have changed since you last saved. Do you still wish to open a new image?',...
                             'Save not up to date',...
                             'Yes','No','Yes');
            switch selection,
                case 'No'
                    return
            end    
        end
        
        %Select filename of nailfold
        [filename pathname] = uigetfile([nailfold_dir '*.bmp'],'File Selector');
        
        if ~filename
            return;
        end
        ui.not_saved = 0;
        
        if length(vessels) %#ok
            %First, clear the axes of exiting lines, pts etc
            delete_vessels;
        end
        if ~isempty(x_vessels);
            delete(ui.gauss_lines);
        end
        
        %cd(pathname);
        ui.file_orig_text = strcat(pathname, filename);       
        h = waitbar(0,'Loading nailfold. Please wait...');
        
        %Set up the save name
        savename = [filename(1:end-4) '_vessels.mat'];
        
        %load in nailfold and compute min max of grey range
        nailfold = imread(ui.file_orig_text);
        if size(nailfold,3) == 3;
            nailfold(:,:,2:3) = [];
        end
%         nailfold = imresize(nailfold, 0.5, 'bilinear');
%         %Contrast enhance the nailfold
%         mask = nailfold < 255;
%         nailfold(~mask) = 0;
%         win_size = 63;
%         nailfold = double(nailfold);
%         
%         local_n = imfilter(double(mask), ones(win_size));
%         local_sum = imfilter(nailfold, ones(win_size));
%         local_sum2 = imfilter(nailfold.^2, ones(win_size));
%         %
%         local_mean = local_sum ./ local_n;
%         local_mean2 = local_sum2 ./ local_n;
%         local_mean(~mask) = 0;
%         local_mean2(~mask) = 0;
% 
%         local_std = sqrt(local_mean2 - local_mean.^2);
%         nailfold = (nailfold - local_mean) ./ local_std;

        max_nailfold = double(max(nailfold(:)));
        min_nailfold = double(min(nailfold(:)));        
        
        %Update axes position give size of loaded nailfold
        %Compute the apsect ratios for these images (they may vary from
        %pair to pair)
        aspect_ratio = size(nailfold,1) / size(nailfold,2);
        
        %Update the size and position and color limits of the axes
        axes_h = axes_pos(3)*aspect_ratio;
        axes_offset = (axes_pos(4)-axes_h) / 2;
        axes_end = axes_offset + axes_h - 1;
        
        set(ui.axes,...
            'Position', axes_pos,...
            'Xlim', [0.5 axes_pos(3)+0.5],...
            'Ylim', [0.5 axes_pos(4)+0.5],...
            'Clim', [min_nailfold max_nailfold]);            
            
        %Make the nailfold visible and colorbar visible, and set the
        %buttondown function
        set(ui.main_fig,...
            'WindowButtonDownFcn', @plot_vessel_Callback,...
            'Name', ['Nailfod annotation tool: ' ui.file_orig_text]);%,...'Pointer','crosshair'
            
        set(ui.nailfold,...
            'Visible', 'on',...
            'CData', nailfold,...
            'XData', [1 axes_pos(3)],...
            'YData', [axes_offset axes_end]); %'ButtonDownFcn', @plot_vessel_Callback

        set(ui.colorbar,...
            'Visible', 'on',...
            'CData', linspace(min_nailfold, max_nailfold, 128)',...
            'XData', [1 colorbar_pos(3)],...
            'YData', [1 colorbar_pos(4)]);
        
        set(ui.colorbar_axes,...
            'Clim', [min_nailfold max_nailfold],...
            'XColor', 'red',...
            'YColor', 'red',...
            'Position', colorbar_pos,...
            'Xtick', [],...
            'Xlim', [0.5 colorbar_pos(3)+0.5],...
            'Ylim', [0.5 colorbar_pos(4)+0.5],...
            'YTickLabel', num2str(linspace(min_nailfold, max_nailfold, length(get(ui.colorbar_axes, 'Ytick')))',3));
        
        set(ui.min_m_slider,...
            'Min', 0,...
            'Max', 0.99*max_nailfold,...
            'Value', 0,...
            'Enable', 'on')
        set(ui.max_m_slider,...
            'Min', 0.01,...
            'Max', max_nailfold,...
            'Value', max_nailfold,...
            'Enable', 'on')
        
        set(ui.min_m_slider_text,...
            'String', ['Select min/max values of mammo contrast range: [' num2str(min_nailfold,3) ', ' num2str(max_nailfold,3) ']']);
        
        %Update widgets
        set(ui.save, 'Enable', 'on');
        set(ui.load, 'Enable', 'on');
        set(ui.close, 'Enable', 'on');
        set(ui.zoom_on, 'Enable', 'on');
        set(ui.pan_on, 'Enable', 'on');
%         set(ui.add_vessel, 'Enable', 'on');
        set(ui.auto, 'Enable', 'on');
        set(ui.show_lines, 'Enable', 'on');
        set(ui.min_m_slider, 'Enable', 'on');
        set(ui.max_m_slider, 'Enable', 'on');
        
        %Do some processing on the nailfold
        [nailfold_strength, nailfold_orientation] = ...
            gaussian_2nd_derivative_line2(nailfold, sigma);
        nailfold_nms = bwareaopen(...
            mb_non_maximal_supp(nailfold_strength, nailfold_orientation) > 0, 10);
        [y_vessels x_vessels] = find(nailfold_nms);
        
        scaling = axes_h / size(nailfold,1);
        x_vessels = x_vessels*scaling;
        y_vessels = y_vessels*scaling + axes_offset - 1;
        
        %Plot these vessels
        ui.gauss_lines = ...
            plot(ui.axes, x_vessels, y_vessels, 'y.', 'markersize', 2, 'HitTest', 'off');
        
        %create a pointer manager so we get crosshairs within the figure
        %axes
        ui.pointer = nan(16);
        ui.pointer([1 9],4:6) = 1; 
        ui.pointer([2 8],3) = 1; ui.pointer([2 8],7) = 1; 
        ui.pointer([3 7], 2) = 1; ui.pointer([3 7],8) = 1; 
        ui.pointer(4:6, [1 9]) = 1; 
        ui.pointer(5,5) = 1;
        ui.enterFcn = @(figHandle, currentPoint)...
        set(figHandle, 'Pointer','custom','PointerShapeCData', ui.pointer,...
            'PointerShapeHotSpot',[5 5]);
        iptSetPointerBehavior(ui.axes, ui.enterFcn);
        iptPointerManager(ui.main_fig);
        
        close(h);
    end

    % --------------------------------------------------------------------
    function file_save_Callback(hObject, eventdata) %#ok %#ok
    % Callback to the "file->save" function used to save the annotations
    % as data in a ".mat" file
        [filename pathname] = uiputfile([anno_dir savename],'File Selector');
        if filename
            
            %Rescale the data to save (but keep copy of original)
            vessels_orig = vessels;
            for ii = 1:length(vessels)
                vessels{ii}(:,1) = vessels{ii}(:,1)/scaling;
                vessels{ii}(:,2) = (vessels{ii}(:,2) + 1 - axes_offset)/scaling;
            end
            
            save(strcat(pathname, filename), 'vessels');
            ui.not_saved = 0;
            vessels = vessels_orig;
        end
    end

    % --------------------------------------------------------------------
    function load_Callback(hObject, eventdata) %#ok  %#ok
    % Callback to the "file->load" function used to load and display a
    % previously saved annotation. This creates new ROI figures containing
    % the annotations, replacing the existing ones if necessary
        if ~isempty(vessels)
            selection = questdlg('Loading an annotation will delete any existing vessels. Do you still want to proceed?',...
                             'Warning',...
                             'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    
                    %First, clear the axes of exiting lines, pts etc
                    delete_vessels
                    
                    %Update widgets 
                    set(ui.show_vessel, 'Enable', 'off');
                    set(ui.show_all, 'Enable', 'off');
                    
                    refreshdata(ui.axes, 'caller');
                case 'No'
                    return
            end
        end
        
        %Now select the new file containing the existing annotation
        [filename pathname] = uigetfile({'*.mat; *.txt'},'File Selector', anno_dir);
        if filename
            if strcmpi(filename(end-3:end), '.mat')
                %load in the vessels
                vessels = u_load(strcat(pathname, filename));
                
            elseif strcmpi(filename(end-3:end), '.txt')
                vessels = read_vessels_from(strcat(pathname, filename));
            else
                warndlg(['Filetype ' filename(end-3:end) ' not recognised.']);
                return;
            end

            set(ui.save, 'Enable', 'on');
            set(ui.zoom_on, 'Enable', 'on');
            set(ui.pan_on, 'Enable', 'on');

            %If the new vessels aren't empty then plot them
            if ~isempty(vessels)
                set(ui.show_vessel, 'Enable', 'on');
                set(ui.show_all, 'Enable', 'on');

                curr_idx = 1;
                for ii = 1:length(vessels)
                    vessels{ii}(:,1) = vessels{ii}(:,1)*scaling;
                    vessels{ii}(:,2) = vessels{ii}(:,2)*scaling + axes_offset - 1;
                    curr_vessel = ii;
                    plot_vessel(ii == 1);
                    refreshdata(ui.axes, 'caller');
                end
                curr_vessel = 1;
            end

            %Refresh the display
            refreshdata(ui.axes, 'caller');
        end
    end

    % --------------------------------------------------------------------
    function close_Callback(hObject, eventdata) %#ok %#ok
    % Callback to the "file->close" function, close the current nailfold
    % and annotation but keeps the program open
        if ui.not_saved
            selection = questdlg('The vessel annotations have changed since you last saved. Do you still wish to close?',...
                             'Save not up to date',...
                             'Yes','No','Yes');
            switch selection,
                case 'No'
                    return
            end    
        end
        
        %First, clear the axes of exiting lines, pts etc
        delete_vessels
        
        %Now clear the nailfold (data and image)
        set(ui.nailfold,'Visible', 'off');
        delete(ui.gauss_lines);
        nailfold = [];
        x_vessels = [];
        y_vessels = [];
        
        %Update widgets ***check this***
        set(ui.save, 'Enable', 'off');
        set(ui.load, 'Enable', 'off');
        set(ui.close, 'Enable', 'off');
        set(ui.zoom_on, 'Enable', 'off');
        set(ui.pan_on, 'Enable', 'off');
%         set(ui.add_vessel, 'Enable', 'off');
        set(ui.auto, 'Enable', 'off');
        set(ui.show_lines, 'Enable', 'off');
        set(ui.min_m_slider, 'Enable', 'off');
        set(ui.max_m_slider, 'Enable', 'off');
        
        %
         set(ui.main_fig,...
            'WindowButtonDownFcn', '',...
            'Name', 'Nailfod annotation tool');
            
    end
    
    % --------------------------------------------------------------------
    function quit_Callback(hObject, eventData) %#ok
    % Callback to the "file->quit" function that exits the program
        if ui.not_saved
            selection = questdlg('The vessel annotations have changed since you last saved. Do you still wish to close?',...
                             'Save not up to date',...
                             'Yes','No','Yes');
            switch selection,
                case 'No'
                    return
            end    
        end
        delete(ui.main_fig);    
    end

   
    % --------------------------------------------------------------------
    function zoom_Callback(hObject, eventdata) %#ok %#ok
    % Callback to the "Zoom" button on the main control panel, allowing the
    % user to zoom in/out of any of 3 ROI figures (zooming occurs
    % simultaneously in all 3 figures)
        if get(ui.zoom_on, 'Value')
            set(ui.pan_on, 'Value', 0);
            axes(ui.axes);  zoom on       
        else
            axes(ui.axes);  zoom off
        end
    end

    % --------------------------------------------------------------------
    function pan_Callback(hObject, eventdata) %#ok %#ok
    % Callback to the "Pan" button on the main control panel, allowing the
    % user to pan around any of 3 ROI figures (panning occurs
    % simultaneously in all 3 figures)
        if get(ui.pan_on, 'Value')
            set(ui.zoom_on, 'Value', 0);
            axes(ui.axes); pan on
        else
            axes(ui.axes); pan off
        end
    end
% --------------------------------------------------------------------
    function min_m_slider_Callback(hObject, eventdata) %#ok
    % Callback...
        min_slider = get(ui.min_m_slider, 'value');
        max_slider = max(min_slider+0.01, get(ui.max_m_slider, 'value'));
        
        set(ui.max_m_slider, 'value', max_slider);
        
        %Compute mammo limits for color scaling
        set(ui.axes, 'Clim', [min_slider max_slider]);
        set(ui.colorbar_axes,...
            'YTickLabel', num2str(linspace(min_slider, max_slider, length(get(ui.colorbar_axes, 'Ytick')))',3));
        
        %
        set(ui.min_m_slider_text,...
            'String', ['Select min/max value of nailfold contrast range: [' num2str(min_slider,3) ', ' num2str(max_slider,3) ']']);
    end
% --------------------------------------------------------------------
    function max_m_slider_Callback(hObject, eventdata) %#ok
    % Callback...
        max_slider = get(ui.max_m_slider, 'value');
        min_slider = min(max_slider-0.01, get(ui.min_m_slider, 'value'));
        
        set(ui.min_m_slider, 'value', min_slider);
        
        %Compute mammo limits for color scaling
        set(ui.axes, 'Clim', [min_slider max_slider]);
        set(ui.colorbar_axes,...
            'YTickLabel', num2str(linspace(min_slider, max_slider, length(get(ui.colorbar_axes, 'Ytick')))',3));
        
        %
        set(ui.min_m_slider_text,...
            'String', ['Select min/max value of nailfold contrast range: [' num2str(min_slider,3) ', ' num2str(max_slider,3) ']']);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vessel callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --------------------------------------------------------------------
    function plot_vessel_Callback(hObject, eventdata) %#ok
    % Callback controlling the actions following various mouse-clicks in
    % the nailfold imaging (except when zooming or panning are enabled)
    
        %Workout the selection type
        selection_type = get(ui.main_fig, 'SelectionType');
        
        %Get location of current point
        XYZ = get(ui.axes, 'CurrentPoint');
        
        switch selection_type,
            
            %If it's a normal click and a vessel exists we plot the point
            %and set continous mode on if the button is held down
            case 'normal'
                
                if curr_idx
                    %plot_pt(XYZ(1,1), XYZ(1,2));
                    set(hObject,'WindowButtonMotionFcn',@plot_continuous);
                    set(hObject,'WindowButtonUpFcn',@continuous_stop);
                end
                
            %If it's a right click we attempt to delete the vessel nearest
            %to the click (asking the user for confirmation)
            case 'alt'
                if ~isempty(vessels)
                    select_vessel(XYZ(1,1), XYZ(1,2));
                    delete_vessel;
                end
            
            %If it's a shift-click we select the nearest vessel as the
            %current vessel
            case 'extend'
                if ~isempty(vessels)
                    select_vessel(XYZ(1,1), XYZ(1,2));
                end
                
            %If it's a double click we create a new vessel and start
            %drawing it    
            case 'open'
                %Turn on the red outline for the old vessel if one exists
                if curr_vessel && curr_vessel <= length(vessels) && ui.show_all_flag
                    set(ui.vessel_outlines(curr_vessel), 'Visible', 'on');
                end
                
                vessels{end + 1} = [];
                curr_vessel = length(vessels);
                curr_idx = 0;
                set(ui.show_vessel, 'Enable', 'on');
                
                plot_pt(XYZ(1,1), XYZ(1,2));
                set(hObject,'WindowButtonMotionFcn',@plot_continuous);
                set(hObject,'WindowButtonUpFcn',@continuous_stop);

            otherwise
                %plot_pt(XYZ(1,1), XYZ(1,2));
        end

        
    end

    % --------------------------------------------------------------------
    function select_vessel(x, y)

        if curr_vessel && curr_vessel <= length(vessels) && ui.show_all_flag
            %Turn the red outline of the original vessel on
            set(ui.vessel_outlines(curr_vessel), 'Visible', 'on');
        end
        
        %Find the nearest vessel to the selected point
        min_dist = inf;
        for ii = 1:length(vessels)
            dists = (vessels{ii}(:,1) - x).^2 + (vessels{ii}(:,2)-y).^2;
            mm = min(dists);
            if  mm < min_dist
                curr_vessel = ii;
                min_dist = mm;
            end
        end
        
        %Turn the red outline of the new vessel off
        set(ui.vessel_outlines(curr_vessel), 'Visible', 'off');
        
        %Set the current index to the end point of the new vessel
        curr_idx = size(vessels{curr_vessel},1);
        
        %Refresh the display
        refreshdata(ui.axes, 'caller');
        set(ui.curr_outline, 'Visible', 'on');
        set(ui.vessel_pts, 'Visible', 'on');
        set(ui.curr_pt, 'Visible', 'on');
    end
            
    % --------------------------------------------------------------------
    function show_vessel_Callback(hObject, eventdata) %#ok
    % Callback to the "show vessel" button, shows/hides the current
    % vessel
        if ui.show_vessel_flag
            set(ui.show_vessel, 'String', 'Show vessel');
            set(ui.curr_outline, 'Visible', 'off');
            ui.show_vessel_flag = 0;
            
        else
            set(ui.show_vessel, 'String', 'Hide vessel');
            set(ui.curr_outline, 'Visible', 'on');
            ui.show_vessel_flag = 1;
        end
    end 

    % --------------------------------------------------------------------
    function show_all_Callback(hObject, eventdata) %#ok
    % Callback to the "Show all" button, shows/hides all the vessels
    
        if ui.show_all_flag
            set(ui.show_all, 'String', 'Show all vessels');
            set(ui.vessel_outlines, 'Visible', 'off');
            ui.show_all_flag = 0;
        else
            set(ui.show_all, 'String', 'Hide all vessels');
            set(ui.vessel_outlines, 'Visible', 'on');
            set(ui.vessel_outlines(curr_vessel), 'Visible', 'off');
            ui.show_all_flag = 1;
        end
    end
    
    % --------------------------------------------------------------------
    function delete_vessel
        
        selection = questdlg('Delete current vessel?',...
                         'Delete current vessel',...
                         'Yes','No','No');
        switch selection,
            case 'Yes',
                
                %Remove the data for the vessel, its handle and the plot
                vessels{curr_vessel} = [0 0];
                
                %If this is the only vessel, delete everything...
                if isempty(vessels) == 1
                    delete_vessels;
                else
                    curr_idx = 0;
                    set(ui.show_vessel, 'Enable', 'off');
                    set(ui.show_vessel, 'String', 'Hide vessel');        
                    
                    ui.show_vessel_flag = 0;                  
                    set(ui.vessel_pts, 'visible', 'off');
                    set(ui.curr_pt, 'visible', 'off');
                    set(ui.curr_outline, 'Visible', 'off');
                end
                
            case 'No',
                return
        end
    end

    % --------------------------------------------------------------------
    function show_lines_Callback(hObject, eventdata) %#ok
    % Callback to the "show vessel" button, shows/hides the current
    % vessel
        if ui.show_lines_flag
            set(ui.show_lines, 'String', 'Show lines');
            set(ui.gauss_lines, 'Visible', 'off');
            ui.show_lines_flag = 0;
            
        else
            set(ui.show_lines, 'String', 'Hide lines');
            set(ui.gauss_lines, 'Visible', 'on');
            ui.show_lines_flag = 1;
        end
    end

    % --------------------------------------------------------------------
    function auto_Callback(hObject, eventdata) %#ok
    % Callback to the "show vessel" button, shows/hides the current
    % vessel
        if ui.auto_flag
            set(ui.auto, 'String', 'Turn Auto mode on');
            ui.auto_flag = 0;       
        else
            set(ui.auto, 'String', 'Turn Auto mode off');
            ui.auto_flag = 1;
        end
    end 

    
    % --------------------------------------------------------------------
    function plot_pt(x, y)
    % Auxiliary function that adds a new point to the current vessel - and if
    % it's the 1st point calls plot_vessel to initialise the plot ui
    
        if ui.auto_flag
            %snap the [x y] point to nearest potential vessel
            neighbours = (abs(x_vessels - x) < 10) & (abs(y_vessels - y) < 10);
            
            if any(neighbours)
                xn = x_vessels(neighbours);
                yn = y_vessels(neighbours);

                dists = (xn - x).^2 + (yn - y).^2;
                [dummy min_idx] = min(dists);
                x = xn(min_idx);
                y = yn(min_idx);
            end
        end
        
        ui.not_saved = 1;
        vessels{curr_vessel}(curr_idx+2:end+1,:) = vessels{curr_vessel}(curr_idx+1:end,:);
        vessels{curr_vessel}(curr_idx+1, :) = [x,y];
        curr_idx = curr_idx + 1;
        
        if size(vessels{curr_vessel},1)==1
            plot_vessel;
            set(ui.vessel_outlines(curr_vessel), 'Visible', 'off');
            set(ui.vessel_pts, 'visible', 'on');
            set(ui.curr_pt, 'visible', 'on');
        end 
        refreshdata(ui.axes, 'caller');
    end

    % --------------------------------------------------------------------
    function plot_continuous(hObject, eventdata) %#ok
       cp = get(ui.axes,'CurrentPoint');
       plot_pt(cp(1,1), cp(1,2));
    end

    % --------------------------------------------------------------------
    function continuous_stop(hObject, eventdata) %#ok
      set(hObject,'WindowButtonMotionFcn','');
      set(hObject,'WindowButtonUpFcn','');
    end

    % --------------------------------------------------------------------
    function plot_vessel(plot_pts)
    % Auxiliary function that...
    
    if nargin < 1
        plot_pts = length(vessels) == 1;
    end
    
        %Create plot ui for the vessel to be shown with all vessels
        ui.vessel_outlines(curr_vessel) =...
            plot(vessels{curr_vessel}(:,1), vessels{curr_vessel}(:,2),...
            'r', 'LineWidth', 1,...
            'XDataSource', ['vessels{', num2str(curr_vessel), '}(:,1)'],...
            'YDataSource', ['vessels{', num2str(curr_vessel), '}(:,2)'],...
            'HitTest', 'off', 'Visible', 'off');
        
        %If this is the first vessel we need to create plot ui for the
        %current vessel points
        if plot_pts
                
            ui.vessel_pts =...
                plot(ui.axes,...
                vessels{curr_vessel}(:,1),...
                vessels{curr_vessel}(:,2), 'mx',...
                'XDataSource', 'vessels{curr_vessel}(:,1)',...
                'YDataSource', 'vessels{curr_vessel}(:,2)',...
                'HitTest', 'off');
            
            ui.curr_outline =...
                plot(ui.axes,...
                vessels{curr_vessel}(:,1),...
                vessels{curr_vessel}(:,2), 'g',...
                'XDataSource', 'vessels{curr_vessel}(:,1)',...
                'YDataSource', 'vessels{curr_vessel}(:,2)',...
                'HitTest', 'off');

            ui.curr_pt =...
                plot(vessels{curr_vessel}(curr_idx,1), vessels{curr_vessel}(curr_idx,1), 'mo',...
                    'XDataSource', 'vessels{curr_vessel}(curr_idx,1)',...
                    'YDataSource', 'vessels{curr_vessel}(curr_idx,2)',...
                    'HitTest', 'off');
            set(ui.show_all, 'Enable', 'on');
        end
        set(ui.curr_outline, 'Visible', 'on');
        set(ui.vessel_pts, 'Visible', 'on');
        set(ui.curr_pt, 'Visible', 'on');
        
        
    end

    % --------------------------------------------------------------------
    
    function delete_vessels
    % Auxiliary function that...
    
        if curr_vessel
    
            %First, clear the axes of exiting lines, pts etc
            delete(ui.vessel_outlines);
            delete(ui.vessel_pts);
            delete(ui.curr_outline);
            delete(ui.curr_pt);

            %Now initialise data and plots to be empty    
            vessels = {};
            ui.vessel_outlines = [];
            ui.vessel_pts = [];
            ui.curr_pt = [];
            ui.curr_outline = [];

            curr_vessel = 0;
            curr_idx = 0;
        end
        
    end
    
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
end
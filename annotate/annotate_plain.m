function annotate_plain(main_fig)
%ANNOTATE_PLAIN
% 
%  overview: the GUI is designed to allow the user to annotate spciulated
%  masses in mammograms. It contains the basic annotation tools without any 
%  of the additional enhanced ROI features in annotate_gui  
%
%  design: the user initially opens the mammogram to be annotated, which is
%  then displayed in the main window. The user annotates the mass by
%  placing points in on the displayed mammogram
%
%  mass outline: the mass outline is defined by a sequence of points, with
%  the outline taken as the polygon connecting each point. For consistency
%  points should be place clockwise
%
%  spicules: spicules associated with masses are defined by a sequence of
%  points marking a line along the length of the spicule. For consistency, 
%  the first point should be taken as the start of the spicule nearest the 
%  centre of the mass
%
%  nipple: in addition the user should mark the position of nipple in the
%  main mammogram
%
%  function architecture: the main GUI function contains three groups of
%  child functions: 1) initialising functions that create the figures, and
%  user controls for the program 2) Callback functions to the user
%  controls, called when any control is activated by the user 3) Auxiliary
%  functions called by the Callbacks that manage the data defining the 
%  annotations 
%
%  data architecture: the data maintained by the function is divided into 2
%  groups 1) data objects defining the annotations (mass_outline, spicules,
%  nipple) 2) data objects defining the graphical output. E.g. the mass
%  outline is stored once in the object mass_outline, which is then
%  associated with 3 graphic objects: the connecting lines, the points and 
%  the circle marking the current point in the main figure window
%
%  See also: ANNOTATE_GUI, ANNOTATE_USER_GUIDE


%Must have window style set to normal
orig_window_style = get(0,'DefaultFigureWindowStyle');
if ~strcmp(orig_window_style, 'normal')
    display('Warning: changing window style to normal for function');
    set(0,'DefaultFigureWindowStyle','normal');
end


handles = {};
plots = {};
mass.border = [];

spicules = {};
cs = 0;
idx = 0;

nipple = [];
initialise_gui;
initialise_figs;

    function initialise_gui
    % creates and positions the user controls on the main annotation panel                   
        handles.gui_fig = figure('Visible','off',...
                                'Name', 'Annotate Masses',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@quit_Callback, ...
                                'WindowStyle', 'normal');
                            
        handles.file_menu = uimenu('Label','File');
        
        handles.file_open = uimenu(handles.file_menu,'Label','Open',...
                                'Callback', @file_open_Callback);
                
        handles.save = uimenu(handles.file_menu,'Label','Save',...
                                'Callback', @file_save_Callback,...
                                'Enable', 'off');
                            
        handles.load = uimenu(handles.file_menu,'Label','Load mass border',...
                                'Callback', @load_Callback,...
                                'Enable', 'off');                    
                    
        handles.close = uimenu(handles.file_menu,'Label', 'Close file',...
                                'Callback', @close_Callback,...
                                'Enable', 'off');
        handles.quit = uimenu(handles.file_menu,'Label', 'Quit',...
                                'Callback', @quit_Callback,...
                                'Separator','on');                    
                            
        %top row of GUI                    
        handles.file_text = uicontrol(...
                                'Style','edit',...
                                'String','No file opened',...
                                'Tag','file_open',...
                                'Position', [5,350, 270, 30]);
       
        handles.image_panel = uipanel(...
                                'Title', 'Image controls',...
                                'Units', 'pixels',...
                                'Position', [5,240,270,105]);
                            
       handles.nipple_button = uicontrol(...
                                'Parent', handles.image_panel,...
                                'Style','togglebutton',...
                                'String','Select nipple ',...
                                'Tag','pan_main',...
                                'Callback', @nipple_Callback,...
                                'Position', [5,50, 80, 40],...
                                'Enable', 'off');
                                
        handles.zoom_main = uicontrol(...
                                'Parent', handles.image_panel,...
                                'Style','togglebutton',...
                                'String','Zoom',...
                                'Tag','zoom_main',...
                                'Callback', @zoom_main_Callback,...
                                'Position', [95,50, 80, 40],...
                                'Enable', 'off');
                            
        handles.pan_main = uicontrol(...
                                'Parent', handles.image_panel,...
                                'Style','togglebutton',...
                                'String','Pan',...
                                'Tag','pan_main',...
                                'Callback', @pan_main_Callback,...
                                'Position', [185,50, 80, 40],...
                                'Enable', 'off');
        
        handles.switch_modes = uibuttongroup(...
                                'Parent', handles.image_panel,...
                                'Units', 'pixels',...
                                'Position', [5, 5, 260, 40],...
                                'SelectionChangeFcn', @switch_modes_Callback,...
                                'SelectedObject', []);
        handles.switch_mode_text = uicontrol(...
                                'Parent', handles.switch_modes,...
                                'Style', 'text',...
                                'Position', [5, 5, 80, 20],...
                                'String', 'Switch modes:');
        handles.switch_button_mass = uicontrol(...
                                'Parent', handles.switch_modes,...
                                'Style', 'Radio',...
                                'String', 'Masses',...
                                'Position', [95, 5, 60, 30],...
                                'Enable', 'off');
        handles.switch_button_spicule = uicontrol(...
                                'Parent', handles.switch_modes,...
                                'Style', 'Radio',...
                                'String', 'Spicules',...
                                'Position', [175, 5, 60, 30],...
                                'Enable', 'off');
                                
        handles.mass_panel = uipanel(...
                                'Title', 'Mass controls',...
                                'Units', 'pixels',...
                                'Position', [5,170,270,70]);
                            
        handles.show_mass = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','togglebutton',...
                                'String','Show mass',...
                                'Tag','show_mass',...
                                'Callback', @show_mass_Callback,...
                                'Position', [5,5, 80, 40],...
                                'Enable', 'off');
        
        handles.delete_point_m = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','pushbutton',...
                                'String','Delete point',...
                                'Tag','delete',...
                                'Callback', @delete_point_m_Callback,...
                                'Position', [95,5, 80, 40],...
                                'Enable', 'off');
        handles.delete_mass = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','pushbutton',...
                                'String','Delete mass',...
                                'Tag','delete_mass',...
                                'Callback', @delete_mass_Callback,...
                                'Position', [185,5, 80, 40],...
                                'Enable', 'off');
                            
        handles.spicules_panel = uipanel(...
                                'Title', 'Spicule controls',...
                                'Units', 'pixels',...
                                'Position', [5,5,270,160]);
        %1st row of spicules panel
        handles.add_spicule = uicontrol(...
                                'Parent', handles.spicules_panel,...
                                'Style','pushbutton',...
                                'String','Add spicule',...
                                'Tag','add_spicule',...
                                'Callback', @add_spicule_Callback,...
                                'Position', [5,95, 80, 40],...
                                'Enable', 'off');
        
        handles.show_spicule = uicontrol(...
                                'Parent', handles.spicules_panel,...
                                'Style','togglebutton',...
                                'String','Show spicule',...
                                'Tag','show_outline',...
                                'Callback', @show_spicule_Callback,...
                                'Position', [95,95, 80, 40],...
                                'Enable', 'off');
        handles.show_all = uicontrol(...
                                'Parent', handles.spicules_panel,...
                                'Style','togglebutton',...
                                'String','Show all',...
                                'Tag','show_all',...
                                'Callback', @show_all_Callback,...
                                'Position', [185,95, 80, 40],...
                                'Enable', 'off');
        %2nd row of spicules panel                    
        handles.switch_spicule = uicontrol(...
                                'Parent', handles.spicules_panel,...
                                'Style','pushbutton',...
                                'String','Switch spicule',...
                                'Tag','switch_spicule',...
                                'Callback', @switch_spicule_Callback,...
                                'Position', [5,50, 80, 40],...
                                'Enable', 'off');               
        handles.delete_point_s = uicontrol(...
                                'Parent', handles.spicules_panel,...
                                'Style','pushbutton',...
                                'String','Delete point',...
                                'Tag','delete',...
                                'Callback', @delete_point_s_Callback,...
                                'Position', [95,50, 80, 40],...
                                'Enable', 'off');
        handles.delete_spicule = uicontrol(...
                                'Parent', handles.spicules_panel,...
                                'Style','pushbutton',...
                                'String','Delete spicule',...
                                'Tag','delete_spicule',...
                                'Callback', @delete_spicule_Callback,...
                                'Position', [185,50, 80, 40],...
                                'Enable', 'off');                         
        %3rd row of spicules panel                    
        handles.delete_all = uicontrol(...
                                'Parent', handles.spicules_panel,...
                                'Style','pushbutton',...
                                'String','Delete all',...
                                'Tag','delete_all',...
                                'Callback', @delete_all_Callback,...
                                'Position', [95,5, 80, 40],...
                                'Enable', 'off');
        handles.not_saved = 0;
        
        set(handles.gui_fig,'Position',[800,300,280,385],'Visible','on');
    end
    
    function initialise_figs
    % Initialises the figure that will contain the mammogram to be
    % annotated
        handles.fig1 = figure(  'Position',[100,0,600,800],...
                                'Visible','off',...
                                'Name', 'Original mammogram',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_fig_Warning, ...
                                'WindowStyle', 'normal');
        handles.axes1 = axes;
        set(handles.fig1,'CurrentAxes', handles.axes1)
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main image callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % --------------------------------------------------------------------
    function zoom_main_Callback(hObject, eventdata)
    % Callback to the "zoom" button allowing the user to zoom in/out the 
    % mammogram 
        if get(handles.zoom_main, 'Value')
            set(handles.pan_main, 'Value', 0);
            axes(handles.axes1);
            zoom on
        else
            axes(handles.axes1);
            zoom off
        end
    end

    % --------------------------------------------------------------------
    function pan_main_Callback(hObject, eventdata)
    % Callback to the "pan" button allowing the user to pan around the 
    % mammogram
        if get(handles.pan_main, 'Value')
            set(handles.zoom_main, 'Value', 0);
            axes(handles.axes1);
            pan on
        else
            axes(handles.axes1);
            pan off
        end
    end

    % --------------------------------------------------------------------
    function nipple_Callback(hObject, eventdata)
    % Callback to the "select nipple" allowing the user to mark the 
    % position of the nipple on the mammogram 
        if get(handles.nipple_button, 'Value')
            set(handles.zoom_main, 'Value', 0);
            set(handles.pan_main, 'Value', 0);
            axes(handles.axes1);
            zoom off
            pan off
            set(handles.fig1, 'WindowButtonDownFcn', @get_nipple,...
                            'WindowButtonUpFcn', []);
        else
            set(handles.fig1, 'WindowButtonDownFcn', @get_R1C1,...
                            'WindowButtonUpFcn', @get_R2C2);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --------------------------------------------------------------------
    function file_open_Callback(hObject, eventdata)
    % Callback to the "file->open" function used to load the mammogram into
    % the program
        
        [filename pathname] = uigetfile('*.bmp','File Selector');
        if filename
            cd(pathname);
            handles.file_orig_text = strcat(pathname, filename);
            set(handles.file_text, 'String', handles.file_orig_text);

            im_orig = imread(handles.file_orig_text);

            set(handles.fig1,   'Visible', 'on',...
                'Name', ['Original mammogram - ', filename]);
            axes(handles.axes1)
            handles.im1 = image(im_orig); axis image; hold on; colormap(gray(256));
            set(handles.im1, 'ButtonDownFcn', @plot_mass_Callback1);

            set(handles.file_open, 'Enable', 'off');
            set(handles.load, 'Enable', 'on');
            set(handles.close, 'Enable', 'on');
            set(handles.save, 'Enable', 'on');
            set(handles.zoom_main, 'Enable', 'on');
            set(handles.pan_main, 'Enable', 'on');
            set(handles.switch_button_mass, 'Enable', 'on');
            set(handles.switch_button_spicule, 'Enable', 'on');

            handles.nipple_point = 0;
            clear im_orig;
        end
    end

    % --------------------------------------------------------------------
    function file_save_Callback(hObject, eventdata)
    % Callback to the "file->save" function used to save the annotations
    % as data in a ".mat" file
        im_orig = imread(handles.file_orig_text);
        [rows cols] = size(im_orig);
        
        if ~isempty(spicules)            
            for ii = 1:length(spicules)
                rm1(ii) = min(spicules(ii).outline(:,2));
                rm2(ii) = max(spicules(ii).outline(:,2));
                cm1(ii) = min(spicules(ii).outline(:,1));
                cm2(ii) = max(spicules(ii).outline(:,1));
            end
            R1 = max(1, round(min([mass.border(:,2); rm1']) - 200));
            R2 = min(rows, round(max([mass.border(:,2); rm2']) + 200));
            C1 = max(1, round(min([mass.border(:,1); cm1']) - 200));
            C2 = min(cols, round(max([mass.border(:,1); cm2']) + 200));
            
            for ii = 1:length(spicules)
                mass_spicules(ii).outline(:,1) = spicules(ii).outline(:,1) - C1 + 1;
                mass_spicules(ii).outline(:,2) = spicules(ii).outline(:,2) - R1 + 1;
            end
        else
            display('no spicules');
            R1 = max(1, round(min(mass.border(:,2)) - 200));
            R2 = min(rows, round(max(mass.border(:,2)) + 200));
            C1 = max(1, round(min(mass.border(:,1)) - 200));
            C2 = min(cols, round(max(mass.border(:,1)) + 200));
            mass_spicules = {};
        end

        mass_ROI = im_orig(R1:R2, C1:C2); clear im_orig;
        mass_outline(:,1) = mass.border(:,1) - C1 + 1; %off set coordinates to
        mass_outline(:,2) = mass.border(:,2) - R1 + 1; %match ROI

        if ~isempty(nipple)
            nipple(1) = nipple(1) - C1 + 1;
            nipple(2) = nipple(2) - R1 + 1;
        end


        [filename pathname] = uiputfile('*.mat','File Selector');
        if filename
            save(strcat(pathname, filename),...
                'mass_outline', 'mass_ROI', 'mass_spicules',...
                'R1', 'R2', 'C1', 'C2', 'nipple');
            cd(pathname);
            handles.not_saved = 0;
        end
        
    end

    function load_Callback(hObject, eventdata)
    % Callback to the "file->load" function used to load and display a
    % previously saved annotation on the current mammogram
        if ~isempty(mass.border) || ~isempty(spicules)
            selection = questdlg('Loading an annotation will create a new ROI and delete any existing mass and spicules. Do you still want to proceed?',...
                             'Warning',...
                             'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    mass.border = [];
                    spicules = {};
                    plots = {};
                    cs = 0;
                    idx = 0;
                    
                    set(handles.show_spicule, 'Enable', 'off');
                    set(handles.switch_spicule, 'Enable', 'off');
                    set(handles.show_all, 'Enable', 'off');
                    set(handles.delete_point_s, 'Enable', 'off');
                    set(handles.delete_all, 'Enable', 'off');
                    set(handles.delete_spicule, 'Enable', 'off');
                
                    set(handles.show_mass, 'Enable', 'off');
                    set(handles.delete_mass, 'Enable', 'off');
                    set(handles.delete_point_m, 'Enable', 'off');
                    set(handles.show_mass, 'Value', 0)

                    refresh_points; display('5');
                case 'No'
                    return
            end
        end
        [filename pathname] = uigetfile('*.mat','File Selector');
        if filename
            S = load(strcat(pathname, filename));

            mass.border = S.mass_outline;
            spicules = S.mass_spicules;
            R1 = S.R1;  
            C1 = S.C1;
            clear S;

            set(handles.save, 'Enable', 'on');
            set(handles.zoom_main, 'Enable', 'on');
            set(handles.pan_main, 'Enable', 'on');
            set(handles.switch_button_mass, 'Enable', 'on');
            set(handles.switch_button_spicule, 'Enable', 'on');

            if ~isempty(spicules)
                for ii = 1:length(spicules)
                    spicules(ii).outline(:,1) = spicules(ii).outline(:,1) + C1 - 1;
                    spicules(ii).outline(:,2) = spicules(ii).outline(:,2) + R1 - 1;
                end

                set(handles.show_all, 'Enable', 'on');
                set(handles.show_all, 'Value', 1)
                idx = 1;
                cs = 1;
                handles.curr_point = spicules(1).outline(1,:);
                plot_points_s
                set(plots.p_s1, 'Visible', 'off');
                set(plots.cp_s1, 'Visible', 'off');


                for cs = 1:length(spicules)
                    plot_lines_s
                end
            end
            if ~isempty(mass.border)
                set(handles.show_mass, 'Enable', 'on');
                set(handles.delete_mass, 'Enable', 'on');
                set(handles.delete_point_m, 'Enable', 'on');
                set(handles.show_mass, 'Value', 1)

                mass.border(:,1) = mass.border(:,1) + C1 - 1; %off set coordinates to
                mass.border(:,2) = mass.border(:,2) + R1 - 1;
                handles.curr_point = mass.border(1,:);
                plot_points_m
            end

            refresh_points
        end
    end

    % --------------------------------------------------------------------
    function close_Callback(hObject, eventdata)
    % Callback to the "file->close" function, close the current mammogram
    % and annotation but keeps the program open
        if handles.not_saved
            selection = questdlg('The mass outline has changed since you last saved. Do you still wish to close?',...
                             'Save not up to date',...
                             'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    delete(handles.fig1); clear handles.fig1;
                    delete(handles.gui_fig); clear handles.gui_fig;
                    clear all;
                    
                    handles = {};
                    plots = {};
                    mass.border = [];
                    spicules = {};
                    cs = 0;
                    idx = 0;
                    
                    initialise_gui;
                    initialise_figs;
                case 'No'
                    return
            end
        else
            delete(handles.fig1); clear handles.fig1;
            delete(handles.gui_fig); clear handles.gui_fig;
            clear all;
            
            handles = {};
            plots = {};
            mass.border = [];
            spicules = {};
            cs = 0;
            idx = 0;

            initialise_gui;
            initialise_figs;
        end
    end
    
    % --------------------------------------------------------------------
    function quit_Callback(hObject, eventData)
    % Callback to the "file->quit" function that exits the program
        if handles.not_saved
            selection = questdlg('The mass outline has changed since you last saved. Do you still wish to quit?',...
                             'Save not up to date',...
                             'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    delete(handles.fig1); clear handles.fig1;
                    delete(handles.gui_fig); clear handles.gui_fig;
                    clear all;
                case 'No'
                    return
            end
        else
            delete(handles.fig1); clear handles.fig1;
            delete(handles.gui_fig); clear handles.gui_fig;
            clear all;
        end
        set(0,'DefaultFigureWindowStyle', orig_window_style);
    end

    % --------------------------------------------------------------------
    function get_nipple(hObject, eventdata)
    % called when the user clicks on the mammogram in "select nipple" mode
        XYZ = get(handles.axes1, 'CurrentPoint');
        nipple = round(XYZ(1,1:2));
        if handles.nipple_point
            refreshdata(handles.axes1, 'caller');
        else
            handles.nipple_point = ...
            plot(nipple(1), nipple(2), 'r*', 'MarkerSize', 12,...
            'XDataSource', 'nipple(1)',...
            'YDataSource', 'nipple(2)',...
            'HitTest', 'off');
        end
    end

    function refresh_points
        %get(handles.axes1, 'Children')
        refreshdata(handles.axes1, 'caller');
    end

    % --------------------------------------------------------------------
    function close_fig_Warning(hObject, eventData)
    % Sets the close request function so that the mammogram figure can only
    % be closed by the main user panel
        warndlg('Individual figures may not be closed. If you would like to close this file, please use the button on the control bar',...
        'Cannot close window')
    end

    % --------------------------------------------------------------------
    function switch_modes_Callback(hObject, eventData)
    % Callback to the button that toggles between masses and spicules mode
        selection = get(get(hObject,'SelectedObject'),'String');
        switch selection,
            case 'Masses',
                set(handles.im1, 'ButtonDownFcn', @plot_mass_Callback1);
                if not(isempty(mass.border))
                    set(handles.show_mass, 'Enable', 'on');
                    set(handles.delete_mass, 'Enable', 'on');
                    set(handles.delete_point_m, 'Enable', 'on');
                    
                    set(plots.p_m1, 'Visible', 'on');
                    
                    set(plots.cp_m1, 'Visible', 'on');
                end
                %turn some spicule buttons off
                set(handles.show_spicule, 'Enable', 'off');
                set(handles.switch_spicule, 'Enable', 'off');
                set(handles.delete_all, 'Enable', 'off');
                set(handles.delete_spicule, 'Enable', 'off');
                set(handles.add_spicule, 'Enable', 'off');
                
                if cs
                    set(plots.p_s1, 'Visible', 'off');
                    set(plots.cp_s1, 'Visible', 'off');
                end
            case 'Spicules'
                set(handles.im1, 'ButtonDownFcn', @plot_spicule_Callback1);                
                if cs
                    set(handles.show_spicule, 'Enable', 'on');
                    set(handles.switch_spicule, 'Enable', 'on');
                    set(handles.show_all, 'Enable', 'on');        
                    set(handles.delete_all, 'Enable', 'on');
                    set(handles.delete_spicule, 'Enable', 'on');
                    
                    set(plots.p_s1, 'Visible', 'on');
                    set(plots.cp_s1, 'Visible', 'on');
                end
                if ~isempty(mass.border)
                    set(plots.p_m1, 'Visible', 'off');
                    set(plots.cp_m1, 'Visible', 'off');
                end
                set(handles.add_spicule, 'Enable', 'on'); 
                set(handles.delete_mass, 'Enable', 'off');
                set(handles.delete_point_m, 'Enable', 'off');
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % --------------------------------------------------------------------
    function show_mass_Callback(hObject, eventdata)
    % Callback to the "Show mass" button that allows the user to show or
    % hide the mass outline
        if get(handles.show_mass, 'Value')
            set(plots.l_m1, 'Visible', 'on');
        else
            set(plots.l_m1, 'Visible', 'off');
        end
    end 

    % --------------------------------------------------------------------
    function delete_point_m_Callback(hObject, eventdata)
    % Callback to the "Delete point" button that delete's the current point
    % from the mass outline
        mass.border(end, :) = [];
        if isempty(mass.border)
            delete(plots.p_m1); clear plots.p_m1;             

            delete(plots.l_m1); clear plots.l_m1;

            delete(plots.cp_m1); clear plots.cp_m1;

            set(handles.show_mass, 'Enable', 'off');
            set(handles.delete_mass, 'Enable', 'off');
            set(handles.delete_point_m, 'Enable', 'off');
            set(handles.show_mass, 'Value', 0)
        else
            handles.curr_point = mass.border(end,:);
            refresh_points
        end
    end

    function delete_mass_Callback(hObject, eventdata)
    % Callback to the "Delete mass" button that deletes all the points from
    % the current mass outline
        selection = questdlg('Delete all points?',...
                         'Delete all points',...
                         'Yes','No','No');
        switch selection,
            case 'Yes',
                mass.border = [];
                
                delete(plots.p_m1); clear plots.p_m1;
                delete(plots.l_m1); clear plots.l_m1;
                delete(plots.cp_m1); clear plots.cp_m1;
               
                set(handles.show_mass, 'Enable', 'off');
                set(handles.delete_mass, 'Enable', 'off');
                set(handles.delete_point_m, 'Enable', 'off');
                set(handles.show_mass, 'Value', 0)
                
                refresh_points;
            case 'No',
                return
        end
    end 
    
    % --------------------------------------------------------------------
    function plot_mass_Callback1(hObject, eventdata)
    % Callback that handles a mouse-click in the mammogram when in
    % mass mode (and not zooming/panning) 
        l_or_r = get(handles.fig1, 'SelectionType');
        XYZ = get(handles.axes1, 'CurrentPoint');
        switch l_or_r,
            case 'normal',
                plot_point_m(XYZ(1,1), XYZ(1,2));
            case 'alt'
                highlight_point_m(XYZ(1,1), XYZ(1,2));
            case 'extend'
                move_point_m(XYZ(1,1), XYZ(1,2));
            otherwise
        end
    end

    % --------------------------------------------------------------------
    function plot_point_m(x, y)
    % Auxiliary function that adds a new point to the mass outline - and if
    % it's the 1st point calls plot_points_m to initialise the plot handles
        handles.not_saved = 1;
        handles.curr_point = [x y];
        if isempty(mass.border)
            mass.border(1,:) = [x,y];
            set(handles.show_mass, 'Enable', 'on');
            set(handles.delete_mass, 'Enable', 'on');
            set(handles.delete_point_m, 'Enable', 'on');
            plot_points_m
        else
            mass.border = [mass.border; [x,y]];
            refresh_points
        end
        
    end

    function highlight_point_m(x, y)
    % Auxiliary function selects the closest point in the mass outline to
    % the new right mouse-click
        d = (mass.border(:,1) - x).^2 + (mass.border(:,2) - y).^2;
        idx = find(d == min(d));
        
        if length(mass.border(:,1)) - idx
            mass.border =...
            [mass.border(idx+1:end,:); mass.border(1:idx,:)];
            handles.curr_point = mass.border(end, :);
        end
        
        refresh_points
    end
    
    % --------------------------------------------------------------------
    function move_point_m(x, y)
    % Auxiliary function that moves the current point in the mass outline
    % to the new mouse-click
        mass.border(end, :) = [x y];
        handles.curr_point = [x y];
        refresh_points;
    end

    % --------------------------------------------------------------------
    function plot_points_m
    % Auxiliary function that initialises the plot handles for 1) the points 
    % in the mass outline; 2) the mass outline; 3) the current point in the
    % mass outline
        %axes1
        axes(handles.axes1);
        plots.p_m1 =...
            plot(mass.border(:,1),...
            mass.border(:,2), 'mx',...
            'XDataSource', 'mass.border(:,1)',...
            'YDataSource', 'mass.border(:,2)',...
            'HitTest', 'off');
        
        plots.l_m1 =...
        plot([mass.border(:,1); mass.border(1,1)],...
            [mass.border(:,2); mass.border(1,2)],...
            'g', 'LineWidth', 1.5,...
            'XDataSource', '[mass.border(:,1); mass.border(1,1)]',...
            'YDataSource', '[mass.border(:,2); mass.border(1,2)]',...
            'HitTest', 'off');
        plots.cp_m1 =...
        plot(handles.curr_point(1), handles.curr_point(2), 'mo',...
            'XDataSource', 'handles.curr_point(1)',...
            'YDataSource', 'handles.curr_point(2)',...
            'HitTest', 'off');
        
        if not(get(handles.show_mass, 'Value'));
            set(plots.l_m1, 'Visible', 'off');
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spicule callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --------------------------------------------------------------------
    function add_spicule_Callback(hObject, eventdata)
    % Callback to the "Add spicule" button to create a new empty spicule
    % object
        spicules(end + 1).outline = [];
        cs = length(spicules);
        set(handles.show_spicule, 'Value', 0);
        set(handles.add_spicule, 'Enable', 'off');
        set(handles.delete_point_s, 'Enable', 'off');
        set(handles.switch_spicule, 'Enable', 'off');
        set(handles.show_spicule, 'Enable', 'off');
        set(handles.show_all, 'Enable', 'off');        
        set(handles.delete_all, 'Enable', 'off');
        set(handles.delete_spicule, 'Enable', 'off');
    end

    % --------------------------------------------------------------------
    function switch_spicule_Callback(hObject, eventdata)
    % Callback to the "Switch spicule" button, allows the user to select 
    % which spicule is the current spicule 
        if cs == length(spicules);
        	cs = 1;
        else
            cs = cs + 1;
        end
        handles.curr_point = spicules(cs).outline(end,:);
        refresh_points;
    end

    % --------------------------------------------------------------------
    function show_spicule_Callback(hObject, eventdata)
    % Callback to the "show spicule" button, shows/hides the current
    % spicule
        if get(handles.show_spicule, 'Value')
            set(plots.l_s1(cs), 'Visible', 'on');
        else
            set(plots.l_s1(cs), 'Visible', 'off');
        end
    end 

    function show_all_Callback(hObject, eventdata)
    % Callback to the "Show all" button, shows/hides all the spicules
        v = get(handles.show_all, 'Value');
        set(handles.show_spicule, 'Value', v);
    
        if v
            set(plots.l_s1(:), 'Visible', 'on');
        else
            set(plots.l_s1(:), 'Visible', 'off');
        end
    end

    % --------------------------------------------------------------------
    function delete_point_s_Callback(hObject, eventdata)
    % Callback to the "Delete point" point button on the spicules panel, 
    % deletes the current point from the current spicule     
        spicules(cs).outline(idx,:) = [];
        idx = idx - 1;
        if isempty(spicules(cs).outline)
            delete_spicule_Callback(hObject, eventdata);
        elseif idx == 0
            idx = 1;
            handles.curr_point = spicules(cs).outline(1,:);
        else
            handles.curr_point = spicules(cs).outline(idx,:);
        end
        refresh_points
    end

    function delete_spicule_Callback(hObject, eventdata)
    % Callback to the "Delete spicule" button, deletes the current spicule
        selection = questdlg('Delete current spicule?',...
                         'Delete current spicule',...
                         'Yes','No','No');
        switch selection,
            case 'Yes',
                spicules(cs) = [];
                delete(plots.l_s1(end)); plots.l_s1(end) = [];
                
                if cs == 1
                    if length(spicules)
                        cs = length(spicules);
                        handles.curr_point = spicules(cs).outline(end,:);
                    else
                        %last spicule deleted, delete all points
                        delete(plots.p_s1); clear plots.p_s1;                
                        delete(plots.l_s1); clear plots.l_s1;
                        delete(plots.cp_s1); clear plots.cp_s1;
                        
                        set(handles.show_spicule, 'Enable', 'off');
                        set(handles.switch_spicule, 'Enable', 'off');
                        set(handles.show_all, 'Enable', 'off');
                        set(handles.delete_point_s, 'Enable', 'off');
                        set(handles.delete_all, 'Enable', 'off');
                        set(handles.delete_spicule, 'Enable', 'off');
                        
                        cs = 0;
                        idx = 0;
                    end
                else
                    cs = cs - 1;                
                	handles.curr_point = spicules(cs).outline(end,:);
                end
                refresh_points
                
            case 'No',
                return
        end
    end

    function delete_all_Callback(hObject, eventdata)
    % Callback to the "Delete all" button, deletes all the spicules
        selection = questdlg('Delete all spicules?',...
                         'Delete all spicules',...
                         'Yes','No','No');
        switch selection,
            case 'Yes',
                spicules = [];
                
                delete(plots.p_s1); clear plots.p_s1;                
                delete(plots.l_s1); clear plots.l_s1;
                delete(plots.cp_s1); clear plots.cp_s1;
                
                set(handles.show_spicule, 'Enable', 'off');
                set(handles.switch_spicule, 'Enable', 'off');
                set(handles.show_all, 'Enable', 'off');
                set(handles.delete_point_s, 'Enable', 'off');
                set(handles.delete_all, 'Enable', 'off');
                set(handles.delete_spicule, 'Enable', 'off');
                
                cs = 0;
                idx = 0;
                refresh_points
                
            case 'No',
                return
        end
    end

    % --------------------------------------------------------------------
    function plot_spicule_Callback1(hObject, eventdata)
    % Callback that handles a mouse-click in the mammogram when in
    % spicule mode (and not zooming/panning)
    
        if cs
            l_or_r = get(handles.fig1, 'SelectionType');
            XYZ = get(handles.axes1, 'CurrentPoint');
            switch l_or_r,
                case 'normal',
                    plot_point_s(XYZ(1,1), XYZ(1,2));
                case 'alt'
                    highlight_point_s(XYZ(1,1), XYZ(1,2));
                case 'extend'
                    move_point_s(XYZ(1,1), XYZ(1,2));
                otherwise
            end
        else
            warndlg('Please add a spicule first');
        end
        
    end

    % --------------------------------------------------------------------
    function plot_point_s(x, y)
    % Auxiliary function that adds a new point to the current spicule - and if
    % it's the 1st point calls plot_points_s to initialise the plot handles    
        handles.not_saved = 1;
        handles.curr_point = [x y];
        
        if isempty(spicules(cs).outline)
            spicules(cs).outline(1,:) = [x,y];
            set(handles.add_spicule, 'Enable', 'on');
            set(handles.delete_point_s, 'Enable', 'on');
            set(handles.switch_spicule, 'Enable', 'on');
            set(handles.show_spicule, 'Enable', 'on');
            set(handles.show_all, 'Enable', 'on');        
            set(handles.delete_all, 'Enable', 'on');
            set(handles.delete_spicule, 'Enable', 'on');
            
            if length(spicules) == 1
                plot_points_s;
            end
            plot_lines_s;
            refresh_points
            idx = 1;
        else
            idx = idx + 1;
            spicules(cs).outline(idx+1:end+1,:) = spicules(cs).outline(idx:end,:);
            spicules(cs).outline(idx, :) = [x,y];
            refresh_points
        end        
    end

    function highlight_point_s(x, y)
    % Auxiliary function selects the closest point in the current spicule to
    % the new right mouse-click
        d = (spicules(cs).outline(:,1) - x).^2 + (spicules(cs).outline(:,2) - y).^2;
        idx = find(d == min(d));
        handles.curr_point = spicules(cs).outline(idx, :);
        refresh_points
    end
    
    % --------------------------------------------------------------------
    function move_point_s(x, y)
    % Auxiliary function that moves the current point in the current
    % spicule to the new mouse-click
        spicules(cs).outline(idx, :) = [x y];
        handles.curr_point = [x y];
        refresh_points;
    end

    % --------------------------------------------------------------------
    function plot_points_s
    % Auxiliary function that initialises the plot handles for 1) the points 
    % in the current spicule; 2) the current point in the current spicule
        %axes1
        axes(handles.axes1);
        plots.p_s1 =...
            plot(spicules(cs).outline(:,1),...
            spicules(cs).outline(:,2), 'mx',...
            'XDataSource', 'spicules(cs).outline(:,1)',...
            'YDataSource', 'spicules(cs).outline(:,2)',...
            'HitTest', 'off');
        
        plots.cp_s1 =...
        plot(handles.curr_point(1), handles.curr_point(2), 'mo',...
            'XDataSource', 'handles.curr_point(1)',...
            'YDataSource', 'handles.curr_point(2)',...
            'HitTest', 'off');
    end

    % --------------------------------------------------------------------
    function plot_lines_s
    % Auxiliary function that initialises the plot handles for the current 
    % spicule
        %axes1
        axes(handles.axes1);
        plots.l_s1(cs) =...
        plot(spicules(cs).outline(:,1), spicules(cs).outline(:,2),...
            'g', 'LineWidth', 1.5,...
            'XDataSource', ['spicules(', num2str(cs), ').outline(:,1)'],...
            'YDataSource', ['spicules(', num2str(cs), ').outline(:,2)'],...
            'HitTest', 'off');
        
        if not(get(handles.show_all, 'Value'));
            set(plots.l_s1(cs), 'Visible', 'off');
        end
    end        

end


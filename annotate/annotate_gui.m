function annotate_gui(varargin)
%  ANNOTATE_GUI - call to open SPAM 
%        (the Software Package for Annotating Mammograms)
%  
%  
%  overview: the GUI is designed to allow the user to annotate spciulated
%  masses in mammograms. It contains a more advanced set of features than
%  annotate_plain, designed to make annotating easier and more accurate
%
%  design: the user initially opens the mammogram to be annotated, which is
%  then displayed in the main window. The user then selects a rectangular
%  ROI containing the mass to be annotated. This ROI is displayed in 3
%  windows displaying the ROI original grayscale (256), contrast enhanced
%  grayscale and colour enhanced. The user then annotates the mass by
%  placing points in any of the 3 windows (each point will appear
%  simultaneously in all 3)
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
%  associated with 9 graphic objects: the connecting lines, the points and 
%  the circle marking the current point in each of 3 ROI windows
%
%  See also: ANNOTATE_PLAIN, ANNOTATE_USER_GUIDE

handles = {};
plots = {};
mass_outline = [];
x_contour = [];
y_contour = [];
spicules = {};
cs = 0;
idx = 0;

rows = 0;
cols = 0;
r1 = 0;
r2 = 0;
c1 = 0;
c2 = 0;
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
                                'Enable', 'off',...
                                'Position', [5,400, 270, 30]);
       
       handles.image_panel = uipanel(...
                                'Title', 'Image controls',...
                                'Units', 'pixels',...
                                'Position', [5,290,270,105]);
                            
       handles.create_roi = uicontrol(...
                                'Parent', handles.image_panel,...
                                'Style','pushbutton',...
                                'String','Create ROI',...
                                'Tag','create_ROI',...
                                'Callback', @create_ROI_Callback,...
                                'Position', [5,50, 80, 40],...
                                'Enable', 'off');
                                
        handles.zoom_on = uicontrol(...
                                'Parent', handles.image_panel,...
                                'Style','togglebutton',...
                                'String','Zoom',...
                                'Tag','zoom_on',...
                                'Callback', @zoom_Callback,...
                                'Position', [95,50, 80, 40],...
                                'Enable', 'off');
                            
        handles.pan_on = uicontrol(...
                                'Parent', handles.image_panel,...
                                'Style','togglebutton',...
                                'String','Pan',...
                                'Tag','pan_on',...
                                'Callback', @pan_Callback,...
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
                                'Position', [5,170,270,120]);
        
        %1st row                    
        handles.show_mass = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','togglebutton',...
                                'String','Show mass',...
                                'Tag','show_mass',...
                                'Callback', @show_mass_Callback,...
                                'Position', [5,55, 80, 40],...
                                'Enable', 'off');
        
        handles.delete_point_m = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','pushbutton',...
                                'String','Delete point',...
                                'Tag','delete',...
                                'Callback', @delete_point_m_Callback,...
                                'Position', [95,55, 80, 40],...
                                'Enable', 'off');
        handles.delete_mass = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','pushbutton',...
                                'String','Delete mass',...
                                'Tag','delete_mass',...
                                'Callback', @delete_mass_Callback,...
                                'Position', [185,55, 80, 40],...
                                'Enable', 'off');
        %2nd row
        handles.contour_level = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','slider',...
                                'Tag','contour_level',...
                                'Min', 0,...
                                'Max', 255,...
                                'SliderStep', [1/255, 10/255],...
                                'Callback', @contour_level_Callback,...
                                'Position', [5,25, 80, 20],...
                                'Enable', 'on');
        handles.contour_tag = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','text',...
                                'String','0',...
                                'Tag','contour_tag',...
                                'Position', [5,5, 80, 20]);                     
        handles.add_contour = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','pushbutton',...
                                'String','Add contour',...
                                'Tag','add_contour',...
                                'Callback', @add_contour_Callback,...
                                'Position', [95,5, 80, 40],...
                                'Enable', 'off');
        handles.show_contour = uicontrol(...
                                'Parent', handles.mass_panel,...
                                'Style','togglebutton',...
                                'String','Show contour',...
                                'Tag','show_contour',...
                                'Callback', @show_contour_Callback,...
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
        
        set(handles.gui_fig,'Position',[800,300,280,435],'Visible','on');
    end
    
    function initialise_figs
    % initialise the figures that will contain the main mammogram and the
    % ROI to be annotated
        handles.fig1 = figure(  'Position',[100,0,600,800],...
                                'Visible','off',...
                                'Name', 'Original mammogram',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_fig_Warning, ...
                                'WindowStyle', 'normal');
        handles.zoom_main = uicontrol(...
                                'Style','togglebutton',...
                                'String','Zoom',...
                                'Tag','zoom_main',...
                                'Callback', @zoom_main_Callback,...
                                'Position', [215,10, 80, 40],...
                                'Enable', 'on');
        handles.pan_main = uicontrol(...
                                'Style','togglebutton',...
                                'String','Pan',...
                                'Tag','pan_main',...
                                'Callback', @pan_main_Callback,...
                                'Position', [305,10, 80, 40],...
                                'Enable', 'on');
        handles.nipple_button = uicontrol(...
                                'Style','togglebutton',...
                                'String','Select nipple ',...
                                'Tag','pan_main',...
                                'Callback', @nipple_Callback,...
                                'Position', [390,10, 80, 40],...
                                'Enable', 'on');
                            
        handles.axes1 = axes;
        set(handles.fig1,'CurrentAxes', handles.axes1)
        
        handles.fig2 = figure(  'Visible','off',...
                                'Name', 'Mass ROI - grayscale',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_fig_Warning, ...
                                'WindowStyle', 'normal');
        handles.axes2 = axes;
        set(handles.fig2,'CurrentAxes', handles.axes2)
        
        handles.fig3 = figure(  'Visible','off',...
                                'Name', 'Mass ROI - color',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_fig_Warning, ...
                                'WindowStyle', 'normal');
        handles.axes3 = axes;                    
        set(handles.fig3,'CurrentAxes', handles.axes3)
        
        handles.fig4 = figure(  'Visible','off',...
                                'Name', 'Mass ROI - enhanced',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_fig_Warning, ...
                                'WindowStyle', 'normal');
        handles.axes4 = axes;                    
        set(handles.fig4,'CurrentAxes', handles.axes4);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main image callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % --------------------------------------------------------------------
    function zoom_main_Callback(hObject, eventdata)
    % Call back to the "Zoom" button on the main mammogram window allowing
    % the user to zoom in/out when selecting the target ROI
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
    % Call back to the "Pan" button on the main mammogram window allowing
    % the user to pan around when selecting the target ROI
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
            [rows cols] = size(im_orig);
            
            set(handles.fig1,   'Visible', 'on',...
                                'WindowButtonDownFcn', @get_R1C1,...
                                'WindowButtonUpFcn', @get_R2C2,...
                                'Name', ['Original mammogram - ', filename]);
            axes(handles.axes1)
            handles.im1 = image(im_orig); axis image; hold on;
            colormap(gray(256));

            set(handles.file_open, 'Enable', 'off');
            set(handles.load, 'Enable', 'on');
            set(handles.close, 'Enable', 'on');

            handles.rec = 0;
            handles.nipple_point = 0;
            handles.roi = [];
            clear im_orig;
        end
    end

    % --------------------------------------------------------------------
    function file_save_Callback(hObject, eventdata)
    % Callback to the "file->save" function used to save the annotations
    % as data in a ".mat" file
        [filename pathname] = uiputfile('*.mat','File Selector');
        if filename
            
            im_orig = imread(handles.file_orig_text);

            if ~isempty(spicules)            
                for ii = 1:length(spicules)
                    rm1(ii) = min(spicules(ii).outline(:,2));
                    rm2(ii) = max(spicules(ii).outline(:,2));
                    cm1(ii) = min(spicules(ii).outline(:,1));
                    cm2(ii) = max(spicules(ii).outline(:,1));
                end
                R1 = max(1, round(min([mass_outline(:,2); rm1']) - 200 + r1));
                R2 = min(rows, round(max([mass_outline(:,2); rm2']) + 200 + r1));
                C1 = max(1, round(min([mass_outline(:,1); cm1']) - 200 + c1));
                C2 = min(cols, round(max([mass_outline(:,1); cm2']) + 200 + c1));

                for ii = 1:length(spicules)
                    mass_spicules(ii).outline(:,1) = spicules(ii).outline(:,1) - C1 + c1;
                    mass_spicules(ii).outline(:,2) = spicules(ii).outline(:,2) - R1 + r1;
                end
            else
                R1 = max(1, round(min(mass_outline(:,2)) - 200 + r1));
                R2 = min(rows, round(max(mass_outline(:,2)) + 200 + r1));
                C1 = max(1, round(min(mass_outline(:,1)) - 200 + c1));
                C2 = min(cols, round(max(mass_outline(:,1)) + 200 + c1));
                mass_spicules = {};
            end

            mass_ROI = im_orig(R1:R2, C1:C2); clear im_orig;
            mass_outline(:,1) = mass_outline(:,1) - C1 + c1; %off set coordinates to
            mass_outline(:,2) = mass_outline(:,2) - R1 + r1; %match ROI

            if ~isempty(nipple)
                nipple(1) = nipple(1) - C1 + 1;
                nipple(2) = nipple(2) - R1 + 1;
            end
            mass.mass_outline = mass_outline;
            mass.mass_ROI = mass_ROI;
            mass.mass_spicules = mass_spicules;
            mass.R1 = R1;
            mass.R2 = R2;
            mass.C1 = C1;
            mass.C2 = C2;
            mass.nipple = nipple;
            mass.name = handles.file_orig_text;
            save(strcat(pathname, filename), 'mass');
            clear mass;
            cd(pathname);
            handles.not_saved = 0;

            mass_outline(:,1) = mass_outline(:,1) + C1 - c1; %off set coordinates to
            mass_outline(:,2) = mass_outline(:,2) + R1 - r1;
        end
    end

    function load_Callback(hObject, eventdata)
    % Callback to the "file->load" function used to load and display a
    % previously saved annotation. This creates new ROI figures containing
    % the annotations, replacing the existing ones if necessary
        if ~isempty(handles.roi)
            selection = questdlg('Loading an annotation will create a new ROI and delete any existing mass and spicules. Do you still want to proceed?',...
                             'Warning',...
                             'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    mass_outline = [];
                    spicules = {};
                    plots = {};
                    cs = 0;
                    idx = 0;
                    
                    delete(handles.rec);
                    if isfield(handles, 'nipple_point')
                        delete(handles.nipple_point);
                        clear('handles.nipple_point');
                    end
                    delete(get(handles.axes2, 'Children'));
                    delete(get(handles.axes3, 'Children'));
                    delete(get(handles.axes4, 'Children'));
                    
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
                    set(handles.show_contour, 'Enable', 'off');
                    set(handles.show_contour, 'Value', 1);
                    
                    refresh_points;
                case 'No'
                    return
            end
        end
        [filename pathname] = uigetfile('*.mat','File Selector');
        if filename
            S = load(strcat(pathname, filename));
            mass = S.mass;
            clear S;
            
            mass_outline = mass.mass_outline;
            spicules = mass.mass_spicules;
            handles.r1 = mass.R1; r1 = handles.r1;  
            handles.r2 = min(rows, mass.R2); r2 = handles.r2;
            handles.c1 = mass.C1; c1 = handles.c1;
            handles.c2 = min(cols, mass.C2); c2 = handles.c2;
            nipple = mass.nipple;
            clear mass;
            
            im_orig = imread(handles.file_orig_text);
            handles.roi = im_orig(r1:r2, c1:c2);
            %handles.linop = [];
            clear im_orig;

            axes(handles.axes1);
            handles.rec = ...
            plot([handles.c1 handles.c2 handles.c2 handles.c1 handles.c1],...
            [handles.r1 handles.r1 handles.r2 handles.r2 handles.r1],...
            'XDataSource', '[handles.c1 handles.c2 handles.c2 handles.c1 handles.c1]',...
            'YDataSource', '[handles.r1 handles.r1 handles.r2 handles.r2 handles.r1]',...
            'HitTest', 'off');
            
            if ~isempty(nipple)
                nipple = nipple + [handles.c1 handles.r1] - 1;
                handles.nipple_point = ...
                plot(nipple(1), nipple(2), 'r*', 'MarkerSize', 12,...
                'XDataSource', 'nipple(1)',...
                'YDataSource', 'nipple(2)',...
                'HitTest', 'off');
            end
        
            set(handles.create_roi, 'Enable', 'on');
            set(handles.save, 'Enable', 'on');
            set(handles.zoom_on, 'Enable', 'on');
            set(handles.pan_on, 'Enable', 'on');
            set(handles.add_contour, 'Enable', 'on');
            set(handles.switch_button_mass, 'Enable', 'on');
            set(handles.switch_button_spicule, 'Enable', 'on');

            set(handles.fig2, 'Position',[100,200,360,480],...
                              'Visible', 'on',...
                              'Name', [handles.file_orig_text, ' - gray']);
            axes(handles.axes2)
            handles.im2 = image(handles.roi); axis image; hold on;
            colormap(gray(256));
            set(handles.im2, 'ButtonDownFcn', @plot_mass_Callback2);

            set(handles.fig3, 'Position',[500,200,360,480],...
                              'Visible', 'on',...
                              'Name', [handles.file_orig_text, ' - gray enhanced']);
            axes(handles.axes3)
            handles.im3 = imagesc(handles.roi); axis image; hold on;
            colormap(gray(256));
            set(handles.im3, 'ButtonDownFcn', @plot_mass_Callback3);

            set(handles.fig4, 'Position',[900,200,360,480],...
                              'Visible', 'on',...
                              'Name', [handles.file_orig_text, ' - colour']);
            axes(handles.axes4)
            handles.im4 = imagesc(handles.roi, 'CData', handles.roi);
            colormap(jet(256));
            axis image; hold on; %{3}%
            set(handles.im4, 'ButtonDownFcn', @plot_mass_Callback4);

            linkaxes([handles.axes2, handles.axes3, handles.axes4]);

            if ~isempty(spicules)
                set(handles.show_all, 'Enable', 'on');
                set(handles.show_all, 'Value', 1)
                idx = 1;
                cs = 1;
                handles.curr_point = spicules(1).outline(1,:);
                plot_points_s
                set(plots.p_s2, 'Visible', 'off');
                set(plots.p_s3, 'Visible', 'off');
                set(plots.p_s4, 'Visible', 'off');
                set(plots.cp_s2, 'Visible', 'off');
                set(plots.cp_s3, 'Visible', 'off');
                set(plots.cp_s4, 'Visible', 'off');

                for cs = 1:length(spicules)
                    plot_lines_s
                end
            end
            if ~isempty(mass_outline)
                set(handles.show_mass, 'Enable', 'on');
                set(handles.delete_mass, 'Enable', 'on');
                set(handles.delete_point_m, 'Enable', 'on');
                set(handles.show_mass, 'Value', 1)

                handles.curr_point = mass_outline(1,:);
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
                    delete(handles.fig2); clear handles.fig2;
                    delete(handles.fig3); clear handles.fig3;
                    delete(handles.fig4); clear handles.fig4;
                    delete(handles.gui_fig); clear handles.gui_fig;
                    clear all;
                    
                    handles = {};
                    plots = {};
                    mass_outline = [];
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
            delete(handles.fig2); clear handles.fig2;
            delete(handles.fig3); clear handles.fig3;
            delete(handles.fig4); clear handles.fig4;
            delete(handles.gui_fig); clear handles.gui_fig;
            clear all;
            
            handles = {};
            plots = {};
            mass_outline = [];
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
                    delete(handles.fig2); clear handles.fig2;
                    delete(handles.fig3); clear handles.fig3;
                    delete(handles.fig4); clear handles.fig4;
                    delete(handles.gui_fig); clear handles.gui_fig;
                    clear all;
                case 'No'
                    return
            end
        else
            delete(handles.fig1); clear handles.fig1;
            delete(handles.fig2); clear handles.fig2;
            delete(handles.fig3); clear handles.fig3;
            delete(handles.fig4); clear handles.fig4;
            delete(handles.gui_fig); clear handles.gui_fig;
            clear all;
        end
    end

    % --------------------------------------------------------------------
    function create_ROI_Callback(hObject, eventdata)
    % Callback to the "Create ROI" button, displays the user selected ROI
    % in the 3 ROI figures, allowing the user to begin annotating the mass
        
        if ~isempty(handles.roi)
            selection = questdlg('Creating a new ROI will delete any mass and spicules. Do you still want to proceed?',...
                         'Warning',...
                         'Yes','No','No');
            switch selection,
                case 'Yes',
                    mass_outline = [];
                    spicules = {};
                    plots = {};
                    cs = 0;
                    idx = 0;
                    
                    delete(get(handles.axes2, 'Children'));
                    delete(get(handles.axes3, 'Children'));
                    delete(get(handles.axes4, 'Children'));
                    
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
                    set(handles.show_contour, 'Enable', 'off');
                    set(handles.show_contour, 'Value', 0);
                    refresh_points;
                case 'No',
                    return
            end
        end
        r1 = handles.r1;    r2 = handles.r2;
        c1 = handles.c1;    c2 = handles.c2;
        
        im_orig = imread(handles.file_orig_text);
        handles.roi = im_orig(r1:r2, c1:c2); 
        clear im_orig; 
        %handles.linop = [];
        
        set(handles.save, 'Enable', 'on');
        set(handles.zoom_on, 'Enable', 'on');
        set(handles.pan_on, 'Enable', 'on');
        set(handles.add_contour, 'Enable', 'on');
        set(handles.switch_button_mass, 'Enable', 'on');
        set(handles.switch_button_spicule, 'Enable', 'on');

        set(handles.fig2, 'Position',[100,200,360,480],...
                          'Visible', 'on',...
                          'Name', [handles.file_orig_text, ' - gray']);
        axes(handles.axes2)
        handles.im2 = image(handles.roi); axis image; hold on;
        colormap(gray(256));
        set(handles.im2, 'ButtonDownFcn', @plot_mass_Callback2);

        set(handles.fig3, 'Position',[500,200,360,480],...
                          'Visible', 'on',...
                          'Name', [handles.file_orig_text, ' - gray enhanced']);
        axes(handles.axes3)
        handles.im3 = imagesc(handles.roi); axis image; hold on;
        colormap(gray(256));
        set(handles.im3, 'ButtonDownFcn', @plot_mass_Callback3);

        set(handles.fig4, 'Position',[900,200,360,480],...
                              'Visible', 'on',...
                              'Name', [handles.file_orig_text, ' - colour']);
        axes(handles.axes4)
        handles.im4 = imagesc(handles.roi, 'CData', handles.roi);
        axis image; hold on; colormap(jet(256));
        set(handles.im4, 'ButtonDownFcn', @plot_mass_Callback4);

        linkaxes([handles.axes2, handles.axes3, handles.axes4]);       
        refresh_points
    end

    % --------------------------------------------------------------------
    function zoom_Callback(hObject, eventdata)
    % Callback to the "Zoom" button on the main control panel, allowing the
    % user to zoom in/out of any of 3 ROI figures (zooming occurs
    % simultaneously in all 3 figures)
        if get(handles.zoom_on, 'Value')
            set(handles.pan_on, 'Value', 0);
            set(handles.switch_button_mass, 'Enable', 'off');
            set(handles.switch_button_spicule, 'Enable', 'off');
            
            axes(handles.axes2);
            zoom on

            axes(handles.axes3);
            zoom on

            axes(handles.axes4);
            zoom on           
        else
            axes(handles.axes2);
            zoom off

            axes(handles.axes3);
            zoom off

            axes(handles.axes4);
            zoom off
            set(handles.switch_button_mass, 'Enable', 'on');
            set(handles.switch_button_spicule, 'Enable', 'on');
        end
    end

    % --------------------------------------------------------------------
    function pan_Callback(hObject, eventdata)
    % Callback to the "Pan" button on the main control panel, allowing the
    % user to pan around any of 3 ROI figures (panning occurs
    % simultaneously in all 3 figures)
        if get(handles.pan_on, 'Value')
            set(handles.zoom_on, 'Value', 0);
            set(handles.switch_button_mass, 'Enable', 'off');
            set(handles.switch_button_spicule, 'Enable', 'off');
            axes(handles.axes2);
            pan on

            axes(handles.axes3);
            pan on

            axes(handles.axes4);
            pan on
        else
            axes(handles.axes2);
            pan off

            axes(handles.axes3);
            pan off

            axes(handles.axes4);
            pan off
            set(handles.switch_button_mass, 'Enable', 'on');
            set(handles.switch_button_spicule, 'Enable', 'on');
        end
    end

    % --------------------------------------------------------------------
    function get_R1C1(hObject, eventdata)
    % Auxiliary function used to determine the top left corner of the
    % rectangular ROI selected in the fullview mammogram
        XYZ = get(handles.axes1, 'CurrentPoint');
        handles.r1 = max(1,round(XYZ(1,2)));
        handles.c1 = max(1,round(XYZ(1,1)));
    end

     % --------------------------------------------------------------------
    function get_R2C2(hObject, eventdata)
    % Auxiliary function used to determine the bottom right corner of the
    % rectangular ROI selected in the fullview mammogram
        XYZ = get(handles.axes1, 'CurrentPoint');
        handles.r2 = min(rows, round(XYZ(1,2)));
        handles.c2 = min(cols, round(XYZ(1,1)));
        
        if handles.rec
            refreshdata(handles.axes1, 'caller');
        else
            handles.rec = ...
            plot([handles.c1 handles.c2 handles.c2 handles.c1 handles.c1],...
            [handles.r1 handles.r1 handles.r2 handles.r2 handles.r1],...
            'XDataSource', '[handles.c1 handles.c2 handles.c2 handles.c1 handles.c1]',...
            'YDataSource', '[handles.r1 handles.r1 handles.r2 handles.r2 handles.r1]',...
            'HitTest', 'off');
            set(handles.create_roi, 'Enable', 'on');
        end
    end

    % --------------------------------------------------------------------
    function get_nipple(hObject, eventdata)
    % Auxiliary function called when the user clicks on the mammogram in 
    % "select nipple" mode
    
    % hObject    handle to plot_point_Callback (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
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
        %get(handles.axes2, 'Children')
        refreshdata(handles.axes2, 'caller');
        refreshdata(handles.axes3, 'caller');
        refreshdata(handles.axes4, 'caller');
    end

    % --------------------------------------------------------------------
    function close_fig_Warning(hObject, eventData)
    % Sets the close request function so that the mammogram figures can only
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
                set(handles.im4, 'CData', handles.roi);
                set(handles.im2, 'ButtonDownFcn', @plot_mass_Callback2);
                set(handles.im3, 'ButtonDownFcn', @plot_mass_Callback3);
                set(handles.im4, 'ButtonDownFcn', @plot_mass_Callback4);
                if not(isempty(mass_outline))
                    set(handles.show_mass, 'Enable', 'on');
                    set(handles.delete_mass, 'Enable', 'on');
                    set(handles.delete_point_m, 'Enable', 'on');
                    
                    set(plots.p_m2, 'Visible', 'on');
                    set(plots.p_m3, 'Visible', 'on');
                    set(plots.p_m4, 'Visible', 'on');
                    
                    set(plots.cp_m2, 'Visible', 'on');
                    set(plots.cp_m3, 'Visible', 'on');
                    set(plots.cp_m4, 'Visible', 'on');
                end
                %turn some spicule buttons off
                set(handles.show_spicule, 'Enable', 'off');
                set(handles.switch_spicule, 'Enable', 'off');
                set(handles.delete_all, 'Enable', 'off');
                set(handles.delete_spicule, 'Enable', 'off');
                set(handles.add_spicule, 'Enable', 'off');
                set(handles.delete_point_s, 'Enable', 'off');
                
                if cs
                    set(plots.p_s2, 'Visible', 'off');
                    set(plots.p_s3, 'Visible', 'off');
                    set(plots.p_s4, 'Visible', 'off');
                    set(plots.cp_s2, 'Visible', 'off');
                    set(plots.cp_s3, 'Visible', 'off');
                    set(plots.cp_s4, 'Visible', 'off');
                end
            case 'Spicules'
                %if isempty(handles.linop)
                %    handles.linop = line_operator_conv(handles.roi, 12, 5, 21);
                %end
                %set(handles.im4, 'CData', handles.linop);
                set(handles.im2, 'ButtonDownFcn', @plot_spicule_Callback2);
                set(handles.im3, 'ButtonDownFcn', @plot_spicule_Callback3);
                set(handles.im4, 'ButtonDownFcn', @plot_spicule_Callback4);
                
                if cs
                    set(handles.show_spicule, 'Enable', 'on');
                    set(handles.switch_spicule, 'Enable', 'on');
                    set(handles.show_all, 'Enable', 'on');        
                    set(handles.delete_all, 'Enable', 'on');
                    set(handles.delete_spicule, 'Enable', 'on');
                    
                    set(plots.p_s2, 'Visible', 'on');
                    set(plots.p_s3, 'Visible', 'on');
                    set(plots.p_s4, 'Visible', 'on');
                    
                    set(plots.cp_s2, 'Visible', 'on');
                    set(plots.cp_s3, 'Visible', 'on');
                    set(plots.cp_s4, 'Visible', 'on');
                end
                if ~isempty(mass_outline)
                    
                    set(plots.p_m2, 'Visible', 'off');
                    set(plots.p_m3, 'Visible', 'off');
                    set(plots.p_m4, 'Visible', 'off');
                    
                    set(plots.cp_m2, 'Visible', 'off');
                    set(plots.cp_m3, 'Visible', 'off');
                    set(plots.cp_m4, 'Visible', 'off');
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
            set(plots.l_m2, 'Visible', 'on');
            set(plots.l_m3, 'Visible', 'on');
            set(plots.l_m4, 'Visible', 'on');
        else
            set(plots.l_m2, 'Visible', 'off');
            set(plots.l_m3, 'Visible', 'off');
            set(plots.l_m4, 'Visible', 'off');
        end
    end 

    % --------------------------------------------------------------------
    function delete_point_m_Callback(hObject, eventdata)
    % Callback to the "Delete point" button that delete's the current point
    % from the mass outline
        mass_outline(end, :) = [];
        if isempty(mass_outline)
            delete(plots.p_m2); clear plots.p_m2;
            delete(plots.p_m3); clear plots.p_m3;
            delete(plots.p_m4); clear plots.p_m4;              

            delete(plots.l_m2); clear plots.l_m2;
            delete(plots.l_m3); clear plots.l_m3;
            delete(plots.l_m4); clear plots.l_m4;

            delete(plots.cp_m2); clear plots.cp_m2;
            delete(plots.cp_m3); clear plots.cp_m3;
            delete(plots.cp_m4); clear plots.cp_m4;

            set(handles.show_mass, 'Enable', 'off');
            set(handles.delete_mass, 'Enable', 'off');
            set(handles.delete_point_m, 'Enable', 'off');
            set(handles.show_mass, 'Value', 0)
        else
            handles.curr_point = mass_outline(end,:);
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
                mass_outline = [];
                
                delete(plots.p_m2); clear plots.p_m2;
                delete(plots.p_m3); clear plots.p_m3;
                delete(plots.p_m4); clear plots.p_m4;              
        
                delete(plots.l_m2); clear plots.l_m2;
                delete(plots.l_m3); clear plots.l_m3;
                delete(plots.l_m4); clear plots.l_m4;
                   
                delete(plots.cp_m2); clear plots.cp_m2;
                delete(plots.cp_m3); clear plots.cp_m3;
                delete(plots.cp_m4); clear plots.cp_m4;
 
                set(handles.show_mass, 'Enable', 'off');
                set(handles.delete_mass, 'Enable', 'off');
                set(handles.delete_point_m, 'Enable', 'off');
                set(handles.show_mass, 'Value', 0)
                
                refresh_points;
            case 'No',
                return
        end
    end

    function contour_level_Callback(hObject, eventdata)
    % Callback to update the displayed value of the contour slide bar
        set(handles.contour_tag, 'String',...
            get(handles.contour_level, 'Value'));
    end

    function add_contour_Callback(hObject, eventdata)
    % Callback to update the displayed value of the contour slide bar
        lev = get(handles.contour_level, 'Value');
        cm = contourc(double(handles.roi), [lev lev]);
        if isempty(cm); return; end;
        
        pts_vec = [];
        idx_vec = [];
        idx = 1;
        while idx < size(cm, 2)
            idx_vec(end+1) = idx;
            pts_vec(end+1) = cm(2, idx);
            idx = idx + cm(2, idx) + 1;
        end

        [m idx] = max(pts_vec);
        idx = idx_vec(idx);
        clear pts_vec idx_vec;
        x_contour = cm(1, idx+1:5:idx+m);
        y_contour = cm(2, idx+1:5:idx+m);
        clear cm;
        
        if isfield(plots, 'contour2')
            refresh_points
        else
            axes(handles.axes2);
            plots.contour2 =...
                plot(x_contour, y_contour,...
                'k.', 'MarkerSize', 0.1,...
                'XDataSource', 'x_contour',...
                'YDataSource', 'y_contour',...
                'HitTest', 'off');
            axes(handles.axes3);
            plots.contour3 =...
                plot(x_contour, y_contour,...
                'k.', 'MarkerSize', 0.1,...
                'XDataSource', 'x_contour',...
                'YDataSource', 'y_contour',...
                'HitTest', 'off');
            axes(handles.axes4);
            plots.contour4 =...
                plot(x_contour, y_contour,...
                'm.', 'MarkerSize', 0.1,...
                'XDataSource', 'x_contour',...
                'YDataSource', 'y_contour',...
                'HitTest', 'off');
        end
        set(handles.show_contour, 'Enable', 'on');
        set(handles.show_contour, 'Value', 1);
        
    end 

    function show_contour_Callback(hObject, eventdata)
    % Callback to update the displayed value of the contour slide bar
        if get(handles.show_contour, 'Value')
            set(plots.contour2, 'Visible', 'on');
            set(plots.contour3, 'Visible', 'on');
            set(plots.contour4, 'Visible', 'on');
        else
            set(plots.contour2, 'Visible', 'off');
            set(plots.contour3, 'Visible', 'off');
            set(plots.contour4, 'Visible', 'off');
        end
    end 
    
    % --------------------------------------------------------------------
    function plot_mass_Callback2(hObject, eventdata)
    % Callback that handles a mouse-click in the 1st ROI figure when in
    % mass mode (and not zooming/panning)
        l_or_r = get(handles.fig2, 'SelectionType');
        XYZ = get(handles.axes2, 'CurrentPoint');
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
    function plot_mass_Callback3(hObject, eventdata)
    % Callback that handles a mouse-click in the 2nd ROI figure when in
    % mass mode (and not zooming/panning)
        l_or_r = get(handles.fig3, 'SelectionType');
        XYZ = get(handles.axes3, 'CurrentPoint');
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
    function plot_mass_Callback4(hObject, eventdata)
    % Callback that handles a mouse-click in the 3rd ROI figure when in
    % mass mode (and not zooming/panning)
        l_or_r = get(handles.fig4, 'SelectionType');
        XYZ = get(handles.axes4, 'CurrentPoint');
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
        if isempty(mass_outline)
            mass_outline(1,:) = [x,y];
            set(handles.show_mass, 'Enable', 'on');
            set(handles.delete_mass, 'Enable', 'on');
            set(handles.delete_point_m, 'Enable', 'on');
            plot_points_m
        else
            mass_outline = [mass_outline; [x,y]];
            refresh_points
        end
        
    end

    function highlight_point_m(x, y)
    % Auxiliary function selects the closest point in the mass outline to
    % the new right mouse-click
        d = (mass_outline(:,1) - x).^2 + (mass_outline(:,2) - y).^2;
        idx = find(d == min(d));
        display([num2str(idx)]);
        if length(mass_outline(:,1)) - idx
            mass_outline =...
            [mass_outline(idx+1:end,:); mass_outline(1:idx,:)];
            handles.curr_point = mass_outline(end, :);
        end
        
        refresh_points
    end
    
    % --------------------------------------------------------------------
    function move_point_m(x, y)
    % Auxiliary function that moves the current point in the mass outline
    % to the new mouse-click
        mass_outline(end, :) = [x y];
        handles.curr_point = [x y];
        refresh_points;
    end

    % --------------------------------------------------------------------
    function plot_points_m
    % Auxiliary function that initialises the plot handles for 1) the points 
    % in the mass outline; 2) the mass outline; 3) the current point in the
    % mass outline in each of the 3 ROI figures
        %axes2
        axes(handles.axes2);
        plots.p_m2 =...
            plot(mass_outline(:,1),...
            mass_outline(:,2), 'mx',...
            'XDataSource', 'mass_outline(:,1)',...
            'YDataSource', 'mass_outline(:,2)',...
            'HitTest', 'off');
        
        plots.l_m2 =...
        plot([mass_outline(:,1); mass_outline(1,1)],...
            [mass_outline(:,2); mass_outline(1,2)],...
            'g', 'LineWidth', 1.5,...
            'XDataSource', '[mass_outline(:,1); mass_outline(1,1)]',...
            'YDataSource', '[mass_outline(:,2); mass_outline(1,2)]',...
            'HitTest', 'off');
        plots.cp_m2 =...
        plot(handles.curr_point(1), handles.curr_point(2), 'mo',...
            'XDataSource', 'handles.curr_point(1)',...
            'YDataSource', 'handles.curr_point(2)',...
            'HitTest', 'off');
        
        %axes3
        axes(handles.axes3);
        plots.p_m3 =...
        plot(mass_outline(:,1), mass_outline(:,2), 'mx',...
            'XDataSource', 'mass_outline(:,1)',...
            'YDataSource', 'mass_outline(:,2)',...
            'HitTest', 'off');
        
        plots.l_m3 = ...
        plot([mass_outline(:,1); mass_outline(1,1)],...
            [mass_outline(:,2); mass_outline(1,2)],...
            'm', 'LineWidth', 1.5,...
            'XDataSource', '[mass_outline(:,1); mass_outline(1,1)]',...
            'YDataSource', '[mass_outline(:,2); mass_outline(1,2)]',...
            'HitTest', 'off');
        plots.cp_m3 =...
        plot(handles.curr_point(1), handles.curr_point(2), 'mo',...
            'XDataSource', 'handles.curr_point(1)',...
            'YDataSource', 'handles.curr_point(2)',...
            'HitTest', 'off');
        
        %axes 4
        axes(handles.axes4);
        plots.p_m4 =...
        plot(mass_outline(:,1), mass_outline(:,2), 'kx',...
            'XDataSource', 'mass_outline(:,1)',...
            'YDataSource', 'mass_outline(:,2)',...
            'HitTest', 'off');        
        plots.l_m4 =...
        plot([mass_outline(:,1); mass_outline(1,1)],...
            [mass_outline(:,2); mass_outline(1,2)],...
            'k', 'LineWidth', 1.5,...
            'XDataSource', '[mass_outline(:,1); mass_outline(1,1)]',...
            'YDataSource', '[mass_outline(:,2); mass_outline(1,2)]',...
            'HitTest', 'off');
        plots.cp_m4 =...
        plot(handles.curr_point(1), handles.curr_point(2), 'ko',...
            'XDataSource', 'handles.curr_point(1)',...
            'YDataSource', 'handles.curr_point(2)',...
            'HitTest', 'off');
        
        if not(get(handles.show_mass, 'Value'));
            set(plots.l_m2, 'Visible', 'off');
            set(plots.l_m3, 'Visible', 'off');
            set(plots.l_m4, 'Visible', 'off');
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
            set(plots.l_s2(cs), 'Visible', 'on');
            set(plots.l_s3(cs), 'Visible', 'on');
            set(plots.l_s4(cs), 'Visible', 'on');
        else
            set(plots.l_s2(cs), 'Visible', 'off');
            set(plots.l_s3(cs), 'Visible', 'off');
            set(plots.l_s4(cs), 'Visible', 'off');
        end
    end 

    function show_all_Callback(hObject, eventdata)
    % Callback to the "Show all" button, shows/hides all the spicules
        v = get(handles.show_all, 'Value');
        set(handles.show_spicule, 'Value', v);
    
        if v
            set(plots.l_s2(:), 'Visible', 'on');
            set(plots.l_s3(:), 'Visible', 'on');
            set(plots.l_s4(:), 'Visible', 'on');
        else
            set(plots.l_s2(:), 'Visible', 'off');
            set(plots.l_s3(:), 'Visible', 'off');
            set(plots.l_s4(:), 'Visible', 'off');
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
                delete(plots.l_s2(end)); plots.l_s2(end) = [];
                delete(plots.l_s3(end)); plots.l_s3(end) = [];
                delete(plots.l_s4(end)); plots.l_s4(end) = [];
                
                if cs == 1
                    if length(spicules)
                        cs = length(spicules);
                        handles.curr_point = spicules(cs).outline(end,:);
                    else
                        %last spicule deleted, delete all points
                        delete(plots.p_s2); clear plots.p_s2;                
                        delete(plots.p_s3); clear plots.p_s3;
                        delete(plots.p_s4); clear plots.p_s4;

                        delete(plots.l_s2); clear plots.l_s2;
                        delete(plots.l_s3); clear plots.l_s3;
                        delete(plots.l_s4); clear plots.l_s4;

                        delete(plots.cp_s2); clear plots.cp_s2;
                        delete(plots.cp_s3); clear plots.cp_s3;
                        delete(plots.cp_s4); clear plots.cp_s4;
                        
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
                
                delete(plots.p_s2); clear plots.p_s2;                
                delete(plots.p_s3); clear plots.p_s3;
                delete(plots.p_s4); clear plots.p_s4;
                
                delete(plots.l_s2); clear plots.l_s2;
                delete(plots.l_s3); clear plots.l_s3;
                delete(plots.l_s4); clear plots.l_s4;
                
                delete(plots.cp_s2); clear plots.cp_s2;
                delete(plots.cp_s3); clear plots.cp_s3;
                delete(plots.cp_s4); clear plots.cp_s4;
                
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
    function plot_spicule_Callback2(hObject, eventdata)
    % Callback that handles a mouse-click in the 1st ROI figure when in
    % spicule mode (and not zooming/panning)
    
        if cs
            l_or_r = get(handles.fig2, 'SelectionType');
            XYZ = get(handles.axes2, 'CurrentPoint');
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
    function plot_spicule_Callback3(hObject, eventdata)
    % Callback that handles a mouse-click in the 2nd ROI figure when in
    % spicule mode (and not zooming/panning)
    
        if cs
            l_or_r = get(handles.fig3, 'SelectionType');
            XYZ = get(handles.axes3, 'CurrentPoint');
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
    function plot_spicule_Callback4(hObject, eventdata)
    % Callback that handles a mouse-click in the 3rd ROI figure when in
    % spicule mode (and not zooming/panning)
    
        if cs
            l_or_r = get(handles.fig4, 'SelectionType');
            XYZ = get(handles.axes4, 'CurrentPoint');
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
    % in each of the 3 ROI figures
    
        %axes2
        axes(handles.axes2);
        plots.p_s2 =...
            plot(spicules(cs).outline(:,1),...
            spicules(cs).outline(:,2), 'mx',...
            'XDataSource', 'spicules(cs).outline(:,1)',...
            'YDataSource', 'spicules(cs).outline(:,2)',...
            'HitTest', 'off');
        
        plots.cp_s2 =...
        plot(handles.curr_point(1), handles.curr_point(2), 'mo',...
            'XDataSource', 'handles.curr_point(1)',...
            'YDataSource', 'handles.curr_point(2)',...
            'HitTest', 'off');
        
        %axes3
        axes(handles.axes3);
        plots.p_s3 =...
        plot(spicules(cs).outline(:,1), spicules(cs).outline(:,2), 'mx',...
            'XDataSource', 'spicules(cs).outline(:,1)',...
            'YDataSource', 'spicules(cs).outline(:,2)',...
            'HitTest', 'off');
        
        plots.cp_s3 =...
        plot(handles.curr_point(1), handles.curr_point(2), 'mo',...
            'XDataSource', 'handles.curr_point(1)',...
            'YDataSource', 'handles.curr_point(2)',...
            'HitTest', 'off');
        
        %axes 4
        axes(handles.axes4);
        plots.p_s4 =...
        plot(spicules(cs).outline(:,1), spicules(cs).outline(:,2), 'kx',...
            'XDataSource', 'spicules(cs).outline(:,1)',...
            'YDataSource', 'spicules(cs).outline(:,2)',...
            'HitTest', 'off');        
        
        plots.cp_s4 =...
        plot(handles.curr_point(1), handles.curr_point(2), 'ko',...
            'XDataSource', 'handles.curr_point(1)',...
            'YDataSource', 'handles.curr_point(2)',...
            'HitTest', 'off');
    end

    % --------------------------------------------------------------------
    function plot_lines_s
    % Auxiliary function that initialises the plot handles for the current 
    % spicule in each of the 3 ROI figures
    
        %axes2
        axes(handles.axes2);
        plots.l_s2(cs) =...
        plot(spicules(cs).outline(:,1), spicules(cs).outline(:,2),...
            'g', 'LineWidth', 1.5,...
            'XDataSource', ['spicules(', num2str(cs), ').outline(:,1)'],...
            'YDataSource', ['spicules(', num2str(cs), ').outline(:,2)'],...
            'HitTest', 'off');
        
        %axes3
        axes(handles.axes3);
        plots.l_s3(cs) = ...
        plot(spicules(cs).outline(:,1), spicules(cs).outline(:,2),...
            'm', 'LineWidth', 1.5,...
            'XDataSource', ['spicules(', num2str(cs), ').outline(:,1)'],...
            'YDataSource', ['spicules(', num2str(cs), ').outline(:,2)'],...
            'HitTest', 'off');
        
        %axes4
        axes(handles.axes4);
        plots.l_s4(cs) =...
        plot(spicules(cs).outline(:,1), spicules(cs).outline(:,2),...
            'k', 'LineWidth', 1.5,...
            'XDataSource', ['spicules(', num2str(cs), ').outline(:,1)'],...
            'YDataSource', ['spicules(', num2str(cs), ').outline(:,2)'],...
            'HitTest', 'off');
        
        if not(get(handles.show_all, 'Value'));
            set(plots.l_s2(cs), 'Visible', 'off');
            set(plots.l_s3(cs), 'Visible', 'off');
            set(plots.l_s4(cs), 'Visible', 'off');
        end
    end

end


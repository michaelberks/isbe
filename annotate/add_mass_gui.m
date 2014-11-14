function add_mass_gui(varargin)
%ANNOTATE_PLAIN
% 
%  overview: the GUI is designed to allow the user to add a synthesised mass
%  to a mammogram.  
%
%  design: the user initially opens the mammogram which is
%  then displayed in the main window. The user then selects a mass, and
%  adds this to the mammogram. The resulting mass region is then depicted
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
%  See also:


%Must have window style set to normal
orig_window_style = get(0,'DefaultFigureWindowStyle');
if ~strcmp(orig_window_style, 'normal')
    display('Warning: changing window style to normal for function');
    set(0,'DefaultFigureWindowStyle','normal');
end


mammogram_fig = []; %handle to mammogram figure
mammogram_im = []; %handle to mammogram image object
mammogram_axes = []; %handle to axes in mammogram figure

mass_fig = []; %handle to mass figure
mass_im = []; %handle to mass image object

new_region_fig1 = []; %handle to 1st new region figure
new_region_im1 = []; %handle to 1st new region image object
new_region_fig2 = []; %handle to 2nd new region figure
new_region_im2 = []; %handle to 2nd new region image object

mass = []; %structure to store mass

handles = []; %structure to store all other GUI handles (buttons etc)
handles.not_saved = 0;

initialise_figs;
    
    function initialise_figs
    % Initialises the figure that will contain the mammogram to be
    % annotated
        mammogram_fig = figure(  'Position',[100,0,600,800],...
                                'Visible','on',...
                                'Name', 'No mammogram loaded',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'WindowStyle', 'normal',...
                                'CloseRequestFcn',@quit_Callback);
                  
        handles.file_menu = uimenu('Label','File');
        
        handles.open_mammogram = uimenu(handles.file_menu,'Label','Open mammogram',...
                                'Callback', @open_mammogram_Callback,...
                                'Enable', 'on');
                            
        handles.open_mass = uimenu(handles.file_menu,'Label','Open mass',...
                                'Callback', @open_mass_Callback,...
                                'Enable', 'on');                    
                
        handles.save = uimenu(handles.file_menu,'Label','Save region',...
                                'Callback', @save_Callback,...
                                'Enable', 'off');
                            
        handles.load = uimenu(handles.file_menu,'Label','Load region',...
                                'Callback', @load_Callback,...
                                'Enable', 'on');                    
                    
        handles.close_mammogram = uimenu(handles.file_menu,'Label', 'Close mammogram',...
                                'Callback', @close_mammogram_Callback,...
                                'Enable', 'off');
                            
        handles.close_mass = uimenu(handles.file_menu,'Label', 'Close mass',...
                                'Callback', @close_mammogram_Callback,...
                                'Enable', 'off');
                            
        handles.quit = uimenu(handles.file_menu,'Label', 'Quit',...
                                'Callback', @quit_Callback,...
                                'Separator','on');
                            
        handles.zoom_main = uicontrol(...
                                'Style','togglebutton',...
                                'String','Zoom',...
                                'Tag','zoom_main',...
                                'Callback', @zoom_main_Callback,...
                                'Position', [215,10, 80, 40],...
                                'Enable', 'off');
        handles.pan_main = uicontrol(...
                                'Style','togglebutton',...
                                'String','Pan',...
                                'Tag','pan_main',...
                                'Callback', @pan_main_Callback,...
                                'Position', [305,10, 80, 40],...
                                'Enable', 'off');
        handles.add_mass_button = uicontrol(...
                                'Style','togglebutton',...
                                'String','Add mass ',...
                                'Tag','add_mass',...
                                'Callback', @add_mass_Callback,...
                                'Position', [390,10, 80, 40],...
                                'Enable', 'off');
          
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main image callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % --------------------------------------------------------------------
    function zoom_main_Callback(hObject, eventdata) %#ok
    % Callback to the "zoom" button allowing the user to zoom in/out the 
    % mammogram 
        if get(handles.zoom_main, 'Value')
            set(handles.pan_main, 'Value', 0);
            set(handles.add_mass_button, 'Value', 0);
            axes(mammogram_axes);
            zoom on
        else
            axes(mammogram_axes);
            zoom off
        end
    end

    % --------------------------------------------------------------------
    function pan_main_Callback(hObject, eventdata) %#ok
    % Callback to the "pan" button allowing the user to pan around the 
    % mammogram
        if get(handles.pan_main, 'Value')
            set(handles.zoom_main, 'Value', 0);
            set(handles.add_mass_button, 'Value', 0);
            axes(mammogram_axes);
            pan on
        else
            axes(mammogram_axes);
            pan off
        end
    end

    % --------------------------------------------------------------------
    function add_mass_Callback(hObject, eventdata) %#ok
    % Callback to the "select nipple" allowing the user to mark the 
    % position of the nipple on the mammogram 
        if get(handles.add_mass_button, 'Value')
            set(handles.zoom_main, 'Value', 0);
            set(handles.pan_main, 'Value', 0);
            axes(mammogram_axes);
            zoom off
            pan off
            set(mammogram_fig, 'WindowButtonDownFcn', @add_mass,...
                            'WindowButtonUpFcn', []);
        else
            set(mammogram_fig, 'WindowButtonDownFcn', [],...
                            'WindowButtonUpFcn', []);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --------------------------------------------------------------------
    function open_mammogram_Callback(hObject, eventdata) %#ok
    % Callback to the "file->open" function used to load the mammogram into
    % the program
        
        %Warn user about unsaved data and perform clean up if we can
        %proceed - else return
        if ~unsaved_region
            return
        end
        
        [filename pathname] = uigetfile('*.bmp','File Selector');
        if filename
            cd(pathname);
            handles.mammogram_name = strcat(pathname, filename);

            mammogram = imread(handles.mammogram_name);
            
            %if mammogram_axes haven't been created yet, do so now
            if isempty(mammogram_axes)
                mammogram_axes = axes;
                set(mammogram_fig, 'CurrentAxes', mammogram_axes);
            else
                %image already exists so delete old image handle first
                delete(mammogram_im);
                set(mammogram_axes, 'NextPlot', 'Replace');
            end
            
            
            axes(mammogram_axes)
            mammogram_im = image(mammogram); axis image; hold on; colormap(gray(256));
            clear mammogram; %Clear the mammogram array from memory
            
            %Set NextPlot to add to allow region to be plotted
            set(mammogram_axes, 'NextPlot', 'Add');
            
            %Put the name of the mammogram in the figure title
            set(mammogram_fig, 'Name', ['Mammogram - ', filename]);
            
            %Enable the user to close the mammogram, zoom/pan
            set(handles.close_mammogram, 'Enable', 'on');
            set(handles.zoom_main, 'Enable', 'on');
            set(handles.pan_main, 'Enable', 'on');
            
            %if mass also exists we can now add mass
            if ~isempty(mass)
                set(handles.add_mass_button, 'Enable', 'on');
            end
            
        end
    end

% --------------------------------------------------------------------
    function open_mass_Callback(hObject, eventdata) %#ok
    % Callback to the "file->open" function used to load the mammogram into
    % the program
        
        %Warn user about unsaved data and perform clean up if we can
        %proceed - else return
        if ~unsaved_region
            return
        end
        
        [filename pathname] = uigetfile('*.mat','File Selector');
        if filename
            cd(pathname);

            mass = u_load([pathname, filename]);
            
            %if mass_fig hasn't been created yet, do so now
            if isempty(mass_im)
                mass_fig = figure(  'Visible','on',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_mass_Callback);
            else
                %image already exists so delete old image handle first
                delete(mass_im);
            end
            %Display the mass
            figure(mass_fig);
            mass_im = imagesc(mass.subtract_ROI); axis image; hold on; colormap(gray(256));
            %Put the name of the mass in the figure title
            set(mass_fig, 'Name', ['Mass - ', filename]);
            %Enable the user to close the mass
            set(handles.close_mass, 'Enable', 'on');
            
            %if mammogram also exists we can now add mass
            if ~isempty(mammogram_axes)
                set(handles.add_mass_button, 'Enable', 'on');
            end
            
        end
    end

    % --------------------------------------------------------------------
    function save_Callback(hObject, eventdata) %#ok
    % Callback to the "file->save" function used to save the annotations
    % as data in a ".mat" file
        
        [filename pathname] = uiputfile('*.mat','File Selector');
        if filename
            save(strcat(pathname, filename), 'mass');
            cd(pathname);
            handles.not_saved = 0;
        end
        
    end

    function load_Callback(hObject, eventdata) %#ok
    % Callback to the "file->load" function used to load and display a
    % previously saved annotation on the current mammogram
       
        %Warn user about unsaved data and perform clean up if we can
        %proceed - else return
        if ~unsaved_region
            return
        end
        
        [filename pathname] = uigetfile('*.mat','File Selector');
        if filename
            
            %if mammogram_axes haven't been created yet, do so now
            if isempty(mammogram_axes)
                mammogram_axes = axes;
                set(mammogram_fig, 'CurrentAxes', mammogram_axes);
            else
                %image already exists so delete old image handle first
                delete(mammogram_im);
            end
            
            %if mass_fig hasn't been created yet, do so now
            if isempty(mass_im)
                mass_fig = figure(  'Visible','on',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_mass_Callback);
            else
                %image already exists so delete old image handle first
                delete(mass_im);
            end
            
            %load the mass
            mass = u_load(strcat(pathname, filename));
            
            %display the mass
            figure(mass_fig);
            mass_im = imagesc(mass.subtract_ROI); axis image; hold on; colormap(gray(256));
            
            %Put the name of the mass in the figure title
            set(mass_fig, 'Name', ['Mass - ', filename]);

            %load the mammogram
            handles.mammogram_name = mass.mammogram_name;
            mammogram = imread(handles.mammogram_name);
            
            %Display the mammogram
            axes(mammogram_axes)
            mammogram_im = image(mammogram); axis image; hold on; colormap(gray(256));
            clear mammogram; %Clear the mammogram array from memory
            
            %Put the name of the mammogram in the figure title
            set(mammogram_fig, 'Name', ['Mammogram - ', handles.mammogram_name]);
            
            %Create the new region figures;
            new_region_fig1 = figure(  'Visible','on',...
                                'Name', 'Mass ROI - enhanced',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_fig_Warning);
        
            new_region_fig2 = figure(  'Visible','on',...
                                'Name', 'Mass ROI - original',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_fig_Warning);  
            
            %Display the new regions
            figure(new_region_fig1);
            new_region_im1 = imagesc(double(uint8(mass.subtract_ROI)) + mass.background_ROI);
            axis image; colormap(gray(256));
            figure(new_region_fig2);
            new_region_im2 = image(double(uint8(mass.subtract_ROI)) + mass.background_ROI);
            axis image; colormap(gray(256));
            
            %Plot the region in mammogram image
            axes(mammogram_axes);
            handles.region_outline = ...
            plot([mass.C1 mass.C2 mass.C2 mass.C1 mass.C1],...
                [mass.R1 mass.R1 mass.R2 mass.R2 mass.R1], 'r:', ...
                'HitTest', 'off');
            
            %Enable the user to close the mass/mammogram, zoom/pan, add
            %mass
            set(handles.close_mass, 'Enable', 'on');
            set(handles.close_mammogram, 'Enable', 'on');
            set(handles.zoom_main, 'Enable', 'on');
            set(handles.pan_main, 'Enable', 'on');
            set(handles.add_mass_button, 'Enable', 'on');
            
        end
    end

    % --------------------------------------------------------------------
    function close_mammogram_Callback(hObject, eventdata) %#ok
    % Callback to the "file->close" function, close the current mammogram
    % and annotation but keeps the program open
    
        %Warn user about unsaved data and perform clean up if we can
        %proceed - else return
        if ~unsaved_region
            return
        end
        
        %Delete the mammogram axes and image object and reset to []
        delete(mammogram_axes);
        mammogram_axes = [];
        mammogram_im = [];
        
        %Disable the user from closing mammogram, zoom/pan and add mass;
        set(handles.close_mammogram, 'Enable', 'off');
        set(handles.zoom_main, 'Enable', 'off');
        set(handles.pan_main, 'Enable', 'off');
        set(handles.add_button, 'Enable', 'off');
        set(mammogram_fig, 'Name', 'No mammogram loaded');
    end

    % --------------------------------------------------------------------
    function close_mass_Callback(hObject, eventdata) %#ok
    % Callback to the "file->close" function, close the current mammogram
    % and annotation but keeps the program open
    
        %Warn user about unsaved data and perform clean up if we can
        %proceed - else return
        if ~unsaved_region
            return
        end
        
        %Delete the mass, fig and image object and reset to []
        delete(mass_fig); %this will delete child object mass_im
        mass_fig = [];
        mass_im = [];
        mass = [];
        
        %Disable the user from closing mass and add mass;
        set(handles.close_mass, 'Enable', 'off');
        set(handles.add_button, 'Enable', 'off');
    end
    
    % --------------------------------------------------------------------
    function quit_Callback(hObject, eventData) %#ok
    % Callback to the "file->quit" function that exits the program
        
        %Warn user about unsaved data and perform clean up if we can
        %proceed - else return
        if ~unsaved_region
            return
        end
        
        delete(mammogram_fig);
        if ~isempty(mass_fig)
            delete(mass_fig);
        end
        if ~isempty(new_region_fig1)
            delete(new_region_fig1);
            delete(new_region_fig2);
        end
        clear all;
        
        %reset default window style
        set(0,'DefaultFigureWindowStyle', orig_window_style);
    end

    % --------------------------------------------------------------------
    function add_mass(hObject, eventdata) %#ok
    % called when the user clicks on the mammogram in "select nipple" mode
        XYZ = get(mammogram_axes, 'CurrentPoint');
        mass_centre = round(XYZ(1,1:2));
        
        [rows_mass cols_mass] = size(mass.subtract_ROI);
        mammogram = double(imread(handles.mammogram_name));
        [rows_mammo cols_mammo] = size(mammogram);
        
        r1 = mass_centre(2) - round(rows_mass/2);
        c1 = mass_centre(1) - round(cols_mass/2);
        r2 = r1 + rows_mass - 1;
        c2 = c1 + cols_mass - 1;
        
        %Check region will fit inside mammogram
        if r1 > 0 && c1 > 0 && r2 <= rows_mammo && c2 <= cols_mammo
            
            %save the background region and limits in mass
            mass.background_ROI = mammogram(r1:r2,c1:c2);
            mammogram(r1:r2,c1:c2) = mammogram(r1:r2,c1:c2) + double(uint8(mass.subtract_ROI));
            
            mass.R1 = r1; clear r1;
            mass.R2 = r2; clear r2;
            mass.C1 = c1; clear c1;
            mass.C2 = c2; clear c2;
            
            %save the mammogram_name in mass
            mass.mammogram_name = handles.mammogram_name;
            
            %delete the old mammogram image, and add the new one
            delete(mammogram_im);
            mammogram_im = image(mammogram); clear mammogram;
            
            %if new region figures don't exist yet then create them
            if isempty(new_region_fig1)
                %Create the new region figures;
                new_region_fig1 = figure(  'Visible','on',...
                                'Name', 'Mass ROI - enhanced',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_fig_Warning);
        
                new_region_fig2 = figure(  'Visible','on',...
                                'Name', 'Mass ROI - original',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'CloseRequestFcn',@close_fig_Warning);
            else
                %they already exist so delete image handle
                delete(new_region_im1);
                delete(new_region_im2);
            end
            
            %Display the new regions
            figure(new_region_fig1);
            new_region_im1 = imagesc(double(uint8(mass.subtract_ROI)) + mass.background_ROI);
            axis image; colormap(gray(256));
            figure(new_region_fig2);
            new_region_im2 = image(double(uint8(mass.subtract_ROI)) + mass.background_ROI);
            axis image; colormap(gray(256));
            
            handles.not_saved = 1;
            set(handles.save, 'Enable', 'on');
            
            %Plot the region outline on the mammogram
            axes(mammogram_axes);
            if isfield(handles, 'region_outline') && ishandle(handles.region_outline)
                delete(handles.region_outline);
            end
            handles.region_outline = ...
            plot([mass.C1 mass.C2 mass.C2 mass.C1 mass.C1],...
                [mass.R1 mass.R1 mass.R2 mass.R2 mass.R1], 'r:', ...
                'HitTest', 'off');
            if isfield(handles, 'region_centre') && ishandle(handles.region_centre)
                delete(handles.region_centre);
            end
            handles.region_centre = ...
            plot(mass_centre(1), mass_centre(2), 'r+', 'HitTest', 'off');
                
        else
            warndlg(['Mass won''t fit using selected centre'...
                ' Please select another position'],...
                'Mass won''t fit!')
        end
    end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxillary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --------------------------------------------------------------------
    function close_fig_Warning(hObject, eventData) %#ok
    % Sets the close request function so that the mammogram figure can only
    % be closed by the main user panel
        warndlg(['Individual figures may not be closed.'...
            'If you would like to close this file,'...
            'please use the button on the control bar'],...
        'Cannot close window')
    end

    % --------------------------------------------------------------------
    function [proceed] = unsaved_region()
    % Sets the close request function so that the mammogram figure can only
    % be closed by the main user panel
        proceed = 1;
        if handles.not_saved
            selection = questdlg(['The current region is unsaved and will be lost.'...
                                ' Do you still want to proceed?'],...
                             'Warning',...
                             'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    %Need to delete new region figures and reset to [];
                    delete(new_region_fig1);
                    delete(new_region_fig2);
                    new_region_fig1 = [];
                    new_region_im1 = [];
                    new_region_fig2 = [];
                    new_region_im2 = [];
                    
                    %Delete region outline from mammogram
                    delete(handles.region_outline);
                    
                    %Delete background from mass
                    mass.background_ROI = [];
                    handles.not_saved = 0;
                    
                case 'No'
                    proceed = 0;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


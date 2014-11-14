function synthetic_texture_roc_gui(varargin)
%SYNTHETIC_TEXTURE_ROC_GUI
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


main_fig = []; %handle to mammogram figure

handles = []; %structure to store all other GUI handles (buttons etc)
curr_image = 1;


%image_set = load(image_set);
number_of_images = 2;
image_set.answers = randsample(2, number_of_images, true);
images_text = ['Images: ', num2str(curr_image), ' of ', num2str(number_of_images)];
user_answers = zeros(number_of_images, 1);

initialise_figs;
    
    function initialise_figs
    % Initialises the figure that will contain the mammogram to be
    % annotated
        main_fig = figure(  'Position',[0,0,1200,600],...
                                'Visible','on',...
                                'Name', 'No mammogram loaded',...
                                'NumberTitle', 'off',...
                                'MenuBar', 'none',...
                                'WindowStyle', 'normal',...
                                'CloseRequestFcn',@quit_Callback);
        movegui(main_fig, 'northwest');
        
        handles.file_menu = uimenu('Label','File');
        
        handles.load_set = uimenu(handles.file_menu,'Label','Load image set',...
                                'Callback', @load_Callback,...
                                'Enable', 'on');
                            
        handles.quit = uimenu(handles.file_menu,'Label', 'Quit',...
                                'Callback', @quit_Callback,...
                                'Separator','on');
                            
        handles.previous_image = uicontrol(...
                                'Style','togglebutton',...
                                'String','Previous',...
                                'Tag','zoom_main',...
                                'Callback', @previous_image_Callback,...
                                'Position', [200,10, 100, 40],...
                                'Enable', 'off');
        handles.next_image = uicontrol(...
                                'Style','togglebutton',...
                                'String','Next',...
                                'Tag','pan_main',...
                                'Callback', @next_image_Callback,...
                                'Position', [305,10, 100, 40],...
                                'Enable', 'on');
        handles.submit_answers = uicontrol(...
                                'Style','togglebutton',...
                                'String','Submit answers',...
                                'Tag','add_mass',...
                                'Callback', @submit_answers_Callback,...
                                'Position', [410,10, 100, 40],...
                                'Enable', 'off');
                            
        handles.select_left_image = uicontrol(...
                                'Style','radio',...
                                'String',[],...
                                'Tag','select_left_image',...
                                'Callback', @select_left_image_Callback,...
                                'Position', [200,60, 40, 40],...
                                'Enable', 'on');
        handles.select_right_image = uicontrol(...
                                'Style','radio',...
                                'String',[],...
                                'Tag','select_right_image',...
                                'Callback', @select_right_image_Callback,...
                                'Position', [410,60, 40, 40],...
                                'Enable', 'on');
        handles.images_text = uicontrol(...
                                'Style','text',...
                                'String', images_text,...
                                'Tag','file_open',...
                                'Position', [245,60, 160, 40]);
                            
        handles.axes_left = axes(...
                                'Units', 'pixels',...
                                'Position', [10,110, 256, 256],...
                                'Visible', 'on');
                            
        handles.axes_right = axes(...
                                'Units', 'pixels',...
                                'Position', [630,110, 256, 256],...
                                'Visible', 'on');
                            
          
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --------------------------------------------------------------------
    function load_Callback(hObject, eventdata) %#ok
    % Callback to the "file->load" function used to load and display a
    % previously saved annotation on the current mammogram
        
    end

    % --------------------------------------------------------------------
    function quit_Callback(hObject, eventdata) %#ok
    % Callback to the "file->load" function used to load and display a
    % previously saved annotation on the current mammogram
        delete(main_fig);
        clear all;
    end

    % --------------------------------------------------------------------
    function previous_image_Callback(hObject, eventdata) %#ok
    % Callback to the "zoom" button allowing the user to zoom in/out the 
    % mammogram
    
        if curr_image == number_of_images
            set(handles.next_image, 'Enable', 'on');
        end
        curr_image = curr_image - 1;
        update_main_fig;
        if curr_image == 1
            set(handles.previous_image, 'Enable', 'off');
        end
        set(handles.previous_image, 'Value', 0);
    end

    % --------------------------------------------------------------------
    function next_image_Callback(hObject, eventdata) %#ok
    % Callback to the "pan" button allowing the user to pan around the 
    % mammogram
        if curr_image == 1
            set(handles.previous_image, 'Enable', 'on');
        end
        
        curr_image = curr_image + 1;
        update_main_fig;
        if curr_image == number_of_images;
            set(handles.next_image, 'Enable', 'off');
        end
        
        set(handles.next_image, 'Value', 0);
    end

    % --------------------------------------------------------------------
    function submit_answers_Callback(hObject, eventdata) %#ok
    % Callback to the "select nipple" allowing the user to mark the 
    % position of the nipple on the mammogram 
        correct = user_answers == image_set.answers;
        score = sum(correct);
        delete(handles.previous_image);
        delete(handles.next_image);
        delete(handles.select_left_image);
        delete(handles.select_right_image);
        delete(handles.submit_answers);
        delete(handles.images_text);
        
        figure(main_fig);
        handles.results_text = uicontrol(...
                                'Style','text',...
                                'String', images_text,...
                                'Tag','file_open',...
                                'Position', [245,200, 160, 160]);
                            
        results_text = {['You got ', num2str(score), ' out of ',...
            num2str(number_of_images), ' correct'], ...
            [num2str(round(100*score/number_of_images)), '%']};
        
        [results_text, newpos] = textwrap( handles.results_text, results_text);
        set(handles.results_text, 'String', results_text, 'Position', newpos);
        
    end

    % --------------------------------------------------------------------
    function select_left_image_Callback(hObject, eventdata) %#ok
        
        if get(handles.select_left_image, 'Value')
            set(handles.select_right_image, 'Value', 0);
            user_answers(curr_image) = 1;
        else
            user_answers(curr_image) = 0;
        end
        
        if all(user_answers)
            set(handles.submit_answers, 'Enable', 'on');
        else
            set(handles.submit_answers, 'Enable', 'off');
        end
    end

    % --------------------------------------------------------------------
    function select_right_image_Callback(hObject, eventdata) %#ok
        
        if get(handles.select_right_image, 'Value')
            set(handles.select_left_image, 'Value', 0);
            user_answers(curr_image) = 2;
        else
            user_answers(curr_image) = 0;
        end
        
        if all(user_answers)
            set(handles.submit_answers, 'Enable', 'on');
        else
            set(handles.submit_answers, 'Enable', 'off');
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxillary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --------------------------------------------------------------------
    function update_main_fig
    % Sets the close request function so that the mammogram figure can only
    % be closed by the main user panel
        
        set(handles.select_left_image, 'Value', 0);
        set(handles.select_right_image, 'Value', 0);
        
        if user_answers(curr_image) == 1
            set(handles.select_left_image, 'Value', 1);
        elseif user_answers(curr_image) == 2
            set(handles.select_right_image, 'Value', 1);
        end
        images_text = ['Images: ', num2str(curr_image), ' of ', num2str(number_of_images)];
        set(handles.images_text, 'String', images_text);
    end

    % --------------------------------------------------------------------
    function [proceed] = unsaved_answers()
    % Sets the close request function so that the mammogram figure can only
    % be closed by the main user panel
        proceed = 1;
        if handles.not_saved
            selection = questdlg(['Unsaved data.'...
                                ' Do you still want to proceed?'],...
                             'Warning',...
                             'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    %Need to delete new region figures and reset to [];
                    
                case 'No'
                    proceed = 0;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [breast_border, breast_air] = adjust_breast_border(main_fig, mammogram, segmentation)

handles = {};
plots = {};

mammogram = double(mammogram);
[breast_border breast_air] = segment_breast_resize(size(mammogram), segmentation);
idx = size(breast_border,1);

initialise_gui;
uiwait(main_fig);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% functions to control the bahaviour of the GUI
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    function initialise_gui
    % creates and positions the user controls on the main annotation panel                   
        handles.main_fig = main_fig;
                            
        handles.control_panel = uipanel(...
                                'Position', [.75,0,.25,1]);
                            
        handles.show_border = uicontrol(...
                                'Parent', handles.control_panel,...
                                'Style','togglebutton',...
                                'String','Show border',...
                                'Tag','show_border',...
                                'Callback', @show_border_Callback,...
                                'Units', 'normalized',...
                                'Position', [.01,.81, .98, .19],...
                                'Value', 1,...
                                'Enable', 'on');
        
        handles.delete_point_m = uicontrol(...
                                'Parent', handles.control_panel,...
                                'Style','pushbutton',...
                                'String','Delete point',...
                                'Tag','delete',...
                                'Callback', @delete_point_m_Callback,...
                                'Units', 'normalized',...
                                'Position', [.01,.61, .98, .19],...
                                'Enable', 'on');
        handles.delete_border = uicontrol(...
                                'Parent', handles.control_panel,...
                                'Style','pushbutton',...
                                'String','Delete border',...
                                'Tag','delete_border',...
                                'Callback', @delete_border_Callback,...
                                'Units', 'normalized',...
                                'Position', [.01,.41, .98, .19],...
                                'Enable', 'on');                    
                            
        handles.zoom_main = uicontrol(...
                                'Parent', handles.control_panel,...
                                'Style','togglebutton',...
                                'String','Zoom',...
                                'Tag','zoom_main',...
                                'Callback', @zoom_main_Callback,...
                                'Units', 'normalized',...
                                'Position', [.01,.21, .49, .19],...
                                'Enable', 'on'); 
                            
        handles.pan_main = uicontrol(...
                                'Parent', handles.control_panel,...
                                'Style','togglebutton',...
                                'String','Pan',...
                                'Tag','pan_main',...
                                'Callback', @pan_main_Callback,...
                                'Units', 'normalized',...
                                'Position', [.5,.21, .49, .19],...
                                'Enable', 'on');
                            
        handles.continue = uicontrol(...
                                'Parent', handles.control_panel,...
                                'Style','togglebutton',...
                                'String','Save and finish',...
                                'Tag','continue',...
                                'Callback', @continue_Callback,...
                                'Units', 'normalized',...
                                'Position', [.01,.01, .98, .19],...
                                'Enable', 'on');                    
        
        handles.fig1 = main_fig;                    
        handles.axes1 = axes(...
            'Position', [0.25 0 0.5 1]);
        
        handles.im1 = imagesc(-mammogram); axis image; hold on; colormap(gray(256));
        set(handles.im1, 'ButtonDownFcn', @plot_border_Callback);
        %handles.curr_point = breast_border(end,:);
        plot_points_m
        handles.not_saved = 0;
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
            axes(handles.axes1);
            zoom on
        else
            axes(handles.axes1);
            zoom off
        end
    end

    % --------------------------------------------------------------------
    function pan_main_Callback(hObject, eventdata) %#ok
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    % --------------------------------------------------------------------
    function continue_Callback(hObject, eventdata) %#ok
    % Callback to the "file->save" function used to save the annotations
    % as data in a ".mat" file
        delete(handles.fig1);
        return;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % --------------------------------------------------------------------
    function show_border_Callback(hObject, eventdata) %#ok
    % Callback to the "Show border" button that allows the user to show or
    % hide the border outline
        if get(handles.show_border, 'Value')
            set(plots.l_m1, 'Visible', 'on');
        else
            set(plots.l_m1, 'Visible', 'off');
        end
    end 

    % --------------------------------------------------------------------
    function delete_point_m_Callback(hObject, eventdata) %#ok
    % Callback to the "Delete point" button that delete's the current point
    % from the border outline
        breast_border(end, :) = [];
        if isempty(breast_border)
            delete(plots.p_m1); clear plots.p_m1;             

            delete(plots.l_m1); clear plots.l_m1;

            delete(plots.cp_m1); clear plots.cp_m1;

            set(handles.show_border, 'Enable', 'off');
            set(handles.delete_border, 'Enable', 'off');
            set(handles.delete_point_m, 'Enable', 'off');
            set(handles.show_border, 'Value', 0)
        else
            %handles.curr_point = breast_border(end,:);
            refreshdata(handles.axes1, 'caller');
        end
    end

    function delete_border_Callback(hObject, eventdata) %#ok
    % Callback to the "Delete border" button that deletes all the points from
    % the current border outline
        selection = questdlg('Delete all points?',...
                         'Delete all points',...
                         'Yes','No','No');
        switch selection,
            case 'Yes',
                breast_border = [];
                
                delete(plots.p_m1); clear plots.p_m1;
                delete(plots.l_m1); clear plots.l_m1;
                delete(plots.cp_m1); clear plots.cp_m1;
               
                set(handles.show_border, 'Enable', 'off');
                set(handles.delete_border, 'Enable', 'off');
                set(handles.delete_point_m, 'Enable', 'off');
                set(handles.show_border, 'Value', 0)
                
                refreshdata(handles.axes1, 'caller');
            case 'No',
                return
        end
    end 
    
    % --------------------------------------------------------------------
    function plot_border_Callback(hObject, eventdata) %#ok
    % Callback that handles a mouse-click in the mammogram when in
    % border mode (and not zooming/panning) 
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
    % Auxiliary function that adds a new point to the border outline - and if
    % it's the 1st point calls plot_points_m to initialise the plot handles
        handles.not_saved = 1;
        %handles.curr_point = [x y];
        if isempty(breast_border)
            breast_border(1,:) = [x,y];
            set(handles.show_border, 'Enable', 'on');
            set(handles.delete_border, 'Enable', 'on');
            set(handles.delete_point_m, 'Enable', 'on');
            plot_points_m
        else
            breast_border = [breast_border; [x,y]];
            refreshdata(handles.axes1, 'caller');
        end
        
    end

    function highlight_point_m(x, y)
    % Auxiliary function selects the closest point in the border outline to
    % the new right mouse-click
        d = (breast_border(:,1) - x).^2 + (breast_border(:,2) - y).^2;
        idx = find(d == min(d));
        
        if length(breast_border(:,1)) ~= idx
            breast_border =...
            [breast_border(idx+1:end,:); breast_border(1:idx,:)];
            %handles.curr_point = breast_border(end, :);
        end
        
        refreshdata(handles.axes1, 'caller');
    end
    
    % --------------------------------------------------------------------
    function move_point_m(x, y)
    % Auxiliary function that moves the current point in the border outline
    % to the new mouse-click
        breast_border(end, :) = [x y];
        %handles.curr_point = [x y];
        refreshdata(handles.axes1, 'caller');
    end

    % --------------------------------------------------------------------
    function plot_points_m
    % Auxiliary function that initialises the plot handles for 1) the points 
    % in the border outline; 2) the border outline; 3) the current point in the
    % border outline
        %axes1
        axes(handles.axes1);
        plots.p_m1 =...
            plot(breast_border(:,1),...
            breast_border(:,2), 'mx',...
            'XDataSource', 'breast_border(:,1)',...
            'YDataSource', 'breast_border(:,2)',...
            'HitTest', 'off');
        
        plots.l_m1 =...
        plot([breast_border(:,1); breast_border(1,1)],...
            [breast_border(:,2); breast_border(1,2)],...
            'g', 'LineWidth', 1.5,...
            'XDataSource', '[breast_border(:,1); breast_border(1,1)]',...
            'YDataSource', '[breast_border(:,2); breast_border(1,2)]',...
            'HitTest', 'off');
        plots.cp_m1 =...
        plot([breast_border(end,1); breast_border(1,1)],...
            [breast_border(end,2); breast_border(1,2)],...
            'r', 'LineWidth', 1.5,...
            'XDataSource', '[breast_border(end,1); breast_border(1,1)]',...
            'YDataSource', '[breast_border(end,2); breast_border(1,2)]',...
            'HitTest', 'off');
        
        if not(get(handles.show_border, 'Value'));
            set(plots.l_m1, 'Visible', 'off');
        end
    end


end


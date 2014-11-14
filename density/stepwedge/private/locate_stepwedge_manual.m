function [step_centres_xy ui errorcheck] = locate_stepwedge_manual(mammo, filmsize,left_breast, num_steps, notch_idx, debug_mode)
%LOCATE_STEPWEDGE_MANUAL method for locating the stepwedge within a
%mammogram given user input
%   [step_xy] = locate_stepwedge_manual(mammo,debug_mode)
%
% Inputs:
%	mammo- Mammogram containing stepwedge
%
%	filmsize - 1 = 18x24, 2 = 24x30 film
%
%	num_steps - number of steps on wedge
%
%	debug_mode- *Insert description of input variable here*
%
%
% Outputs:
%	step_centres_xy - Coordinates of the centre of each step
%
%	ui.step_pos_xy - Coordinates of the initial rectangle selected by the
%      user
%
%	ui.step_4_xy - Coordinates of the initial 4 points selected by the user
%
%   ui.step_nums - Numbers of the 2 steps selected by the user
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Sep-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if nargin < 6
    debug_mode = 0;
end

if left_breast
    mammo = rot90(mammo, 2);
end

%Create a figure to display stepwedge
step_fig = figure(...
    'Name','Select stepwedge region',...
    'Units', 'normalized',...
    'OuterPosition', [0 0 1 1]);
subplot(1,2,1);
imagesc(mammo); axis image; colormap(gray(256));
title('Select stepwedge region');
set(gca, 'xticklabel', [], 'yticklabel', []);

%Check size of film and choose default location for box
if filmsize == 1 %18x24 mammogram
    y_start = size(mammo,1) - 599;%150;
    y_end = 600;
else %filmsize == 2 %24x30 mammogram
    y_start = size(mammo,1) - 599;
    y_end = 600;
end
x_start = 1;%75;
x_end = 100;%100;

%Draw resizable rectangle in initial position, this cane resized by the
%user
rect_h = imrect(gca, [x_start y_start x_end y_end]);
api = iptgetapi(rect_h);
pos = [x_start y_start x_end y_end];
fcn = makeConstrainToRectFcn('imrect', get(gca,'XLim'), get(gca,'YLim'));

api.setPositionConstraintFcn(fcn);
api.addNewPositionCallback(@(p) select_stepwedge(p));

figure(step_fig); subplot(1,2,2);
imagesc(mammo(pos(2):pos(2)+pos(4)-1, pos(1):pos(1)+pos(3)-1));
axis image;

yes = 'Yes';
no = 'Stepwedge missing';
answer = questdlg('Is the stepwedge present?', 'Stepwedge', yes, no, yes);

errorcheck = 0;
if strcmpi(answer, no)
    step_centres_xy = []; 
    ui = [];
    errorcheck = 1;
    close(step_fig);
    return;
end

% Query the user if the current position is ok, when they click yes we
% continue execution
region_ok = helpdlg('Is the region OK?', 'Select stepwedge');
waitfor(region_ok);
                    
step_pos = round(api.getPosition());
clear api;

%Extract box defined by user into sw_image
sw_image = rot90(mammo(step_pos(2):step_pos(2)+step_pos(4)-1, step_pos(1):step_pos(1)+step_pos(3)-1));

%Now go through process with user of finding steps
ok = -3;
while (ok < 0)

    if(ok == -3)
        
        %Display the stepwedge to the user and ask them to select 4 points
        set(step_fig, 'Name',...
            'Select steps');
        subplot(1,1,1);
        imagesc(sw_image); axis image; colormap(gray(256));
        title('Select the left and right edges of two steps, then press return (4 points in all)');
        set(gca, 'xticklabel', [], 'yticklabel', []);
        
        %Check they've selected 4 points, otherwise ask them what's wrong
        problem = 1;
        while problem
            [user_pts_x, user_pts_y, p] = impixel; %#ok
            problem = length(user_pts_x) ~= 4;
            if problem
                new_roi = 'Select new stepwedge region';
                new_steps = 'Select new step points';
                answer = questdlg('Incorrect number of points selected. Please either:',...
                    'Stepwedge error', new_roi, new_steps, new_roi);

                if strcmpi(answer, new_roi);
                    close(step_fig);
                    [step_centres_xy] = locate_stepwedge_manual(mammo, filmsize,left_breast, num_steps, notch_idx, debug_mode);
                    return;
                end
            end
        end                
        
        %sort points by x-vals in case user clicked in the wrong order
        [user_pts_x idx] = sort(user_pts_x);
        user_pts_y = user_pts_y(idx);
        
        %plot the points so they stay on the mammo
        hold on;
        plot(user_pts_x, user_pts_y, 'r+', 'MarkerSize', 4);
        
        swx=zeros(2,1);
        swy=zeros(2,1);
    
        swx(1) = (user_pts_x(1) + user_pts_x(2))/2;
        swx(2) = (user_pts_x(3) + user_pts_x(4))/2;
        swy(1) = (user_pts_y(1) + user_pts_y(2))/2;
        swy(2) = (user_pts_y(3) + user_pts_y(4))/2;
        
        % get step numbers but try to guard against mistakes
        problem = 1;
        while problem
            
            prompt = {'What number is the first step?','What number is the second step?:'};
            dlg_title = 'Select step numbers';
            num_lines = 1;
            screen_size = get(0,'ScreenSize');
            options.Position = [screen_size(3)/4, 3*screen_size(4)/5];
            answer = inputdlg(prompt,dlg_title,num_lines,{'',''}, options);
            
            i_min = str2num(answer{1}); %#ok
            i_max = str2num(answer{2}); %#ok
            
            %if step numbers accidentally entered wrong way round, don't crash!
            if(i_min > i_max)
                i_temp = i_max;
                i_max = i_min;
                i_min = i_temp;
            end
            
            %check steps are numbers (if non-numeric input is given str2num
            %returns []), and are in a suitable range
            problem = isempty(i_min) || isempty(i_max) ||...
                i_min<1 || i_min>num_steps || i_max<1 || i_max>num_steps;
            
            if problem
                h = warndlg('Please choose again',...
                    'Invalid step numbers', 'modal');
                uiwait(h);
            end
            
        end
    
    elseif(ok==-1)
        %Nudge steps down
        i_min = i_min+5;
        i_max = i_max+5;
        
        figure(step_fig);
        imagesc(sw_image); axis image; colormap(gray(256)); axis off;
        
    elseif(ok==-2)
        %Nudge steps up
        i_min = i_min-5;
        i_max = i_max-5;
        
        figure(step_fig);
        imagesc(sw_image); axis image; colormap(gray(256)); axis off;    
    end    
    numsteps_selected = i_max - i_min + 1;

    %Work out the x and y offsets between the steps
    x_off = (swx(2)-swx(1))/(numsteps_selected-1);
    y_off = (swy(2)-swy(1))/(numsteps_selected-1);
    
    %Calculate the centre of each step on the wedge
    step_centres_x = swx(1) + ((1:num_steps) - i_min)*x_off;
    step_centres_y = swy(1) + ((1:num_steps) - i_min)*y_off;
    
    %Display the step centres to the user
    plot(step_centres_x, step_centres_y, 'bs', 'MarkerSize', 15);
    plot(step_centres_x(notch_idx(1:end-1)), step_centres_y(notch_idx(1:end-1)), 'rs', 'MarkerSize', 15);
    plot(step_centres_x(notch_idx(end)), step_centres_y(notch_idx(end)), 'gs', 'MarkerSize', 15);
    
    xlabel(['Red boxes should match widely spaced notches,',...
        'the green box should match the far right notch']);
    
    %Check whether steps need adjusting
    dlg = create_choice_dialog;
    uiwait(dlg);
    %will return ok = ...
    % 0 for yes points good,
    % -3 for no, choose points again
    % -2 for nudge steps up
    % -1 for nudge steps down
    % Note option of 2 no longer selectable

end
close(step_fig);
ui.step_pos_xy = step_pos;
ui.step_4_xy = [user_pts_x, user_pts_y];
ui.step_nums = [i_min i_max];

%Finally, we need to make the step centres relative to the original
%mammogram: x = step_pos(1)+step_pos(3)-yi; y = step_pos(2)+xi-1;
step_centres_xy = [step_pos(1) + step_pos(3) - step_centres_y'...
                   step_pos(2) + step_centres_x' - 1];
               
%For left breasts make the step_centres relative to the original mammo
if left_breast
    [r c] = size(mammo);
    step_centres_xy(:,1) = c - step_centres_xy(:,1) + 1;
    step_centres_xy(:,2) = r - step_centres_xy(:,2) + 1;
    if debug_mode
        mammo = rot90(mammo, 2);
    end
end
    
               
%if in debug mode, check the steps match up with the original mammo
if debug_mode
    figure('Name', 'Stepwedge located'); 
    imagesc(mammo); axis image; colormap(gray(256)); hold on;
    plot(step_centres_xy(:,1), step_centres_xy(:,2), 'bs', 'MarkerSize', 15);
    plot(step_centres_xy(notch_idx(1:end-1),1), step_centres_xy(notch_idx(1:end-1),2), 'rs', 'MarkerSize', 15);
    plot(step_centres_xy(notch_idx(end),1), step_centres_xy(notch_idx(end),2), 'gs', 'MarkerSize', 15);
    title('Mammogram with stepwedge located');
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    function dlg = create_choice_dialog

        qstring = 'Are these regions ok?';
        choose_again = 'No, choose new pts';
        nudge_u = 'Nudge 5 steps left';
        nudge_d = 'Nudge 5 steps right';
        yes = 'Yes';
        
        screen_size = get(0,'ScreenSize');
        
        dlg = figure(...
            'Position',[screen_size(3)/4 - 165,3*screen_size(4)/5 - 45,330,95],...
            'Visible','on',...
            'Name', 'Selected steps',...
            'WindowStyle', 'normal',...
            'NumberTitle', 'off'...
            );
        uicontrol(...
            'Style','text',...
            'String', qstring,...
            'Position',[15 70 300 25],...
            'BackgroundColor', get(dlg, 'Color'),...
            'Parent',dlg);
        % Create three radio buttons in the button group.
        uicontrol(...
            'Style','ToggleButton',...
            'String', yes,...
            'Position', [10 45 150 25],...
            'Parent', dlg,...
            'Callback', @yes_Callback,...
            'HandleVisibility','off');
        uicontrol(...
            'Style','ToggleButton',...
            'String', choose_again,...
            'Position',[170 45 150 25],...
            'parent', dlg,...
            'Callback', @no_Callback,...
            'HandleVisibility','off');
        uicontrol(...
            'Style','ToggleButton',...
            'String', nudge_u,...
            'Position',[10 10 150 25],...
            'Parent', dlg,...
            'Callback', @nudge_u_Callback,...
            'HandleVisibility','off');
        uicontrol(...
            'Style','ToggleButton',...
            'String', nudge_d,...
            'Position',[170 10 150 25],...
            'Parent', dlg,...
            'Callback', @nudge_d_Callback,...
            'HandleVisibility','off');
    end

    % --------------------------------------------------------------------
    function yes_Callback(hObject, eventdata) %#ok     
        ok = 0;
        delete(dlg);
    end
    % --------------------------------------------------------------------
    function no_Callback(hObject, eventdata) %#ok     
        ok = -3;
        delete(dlg);
    end
    % --------------------------------------------------------------------
    function nudge_u_Callback(hObject, eventdata) %#ok     
        ok = -1;
        delete(dlg);
    end
    % --------------------------------------------------------------------
    function nudge_d_Callback(hObject, eventdata) %#ok     
        ok = -2;
        delete(dlg);
    end
    % --------------------------------------------------------------------
    function select_stepwedge(pos)
        pos = round(pos);
        figure(step_fig); subplot(1,2,2);
        imagesc(mammo(pos(2):pos(2)+pos(4)-1, pos(1):pos(1)+pos(3)-1));
        axis image; axis off; colormap(gray(256));
        figure(region_ok);
    end
end

% function wedgevals = wedgecal_jen(filename)
% 
% % when calling the programme as part of 'calibrate.m'
% % don't forget to change (filename) to (IMAGE, map)
% % and comment out line below
% 
% [IMAGE,map] = imread(filename);

function [wedgevals swx swy errorcheck] = wedgecal_jen(IMAGE, debug_mode)

if nargin < 2
    debug_mode = 0;
end

numsteps = 39;

errorcheck = 0;

%Create a figure to display stepwedge
step_fig = figure('Name','Select stepwedge region');
subplot(1,2,1);
imagesc(IMAGE); axis image; colormap(gray(256));
title('Select stepwedge region');
set(gca, 'xticklabel', [], 'yticklabel', []);

%Need to check Matlab version here
v = version('-release');
if strcmpi(v, '2007b')
    
    rect_h = imrect(gca, [75 150 100 size(IMAGE,1)-200]);
    api = iptgetapi(rect_h);
    pos = [75 150 100 size(IMAGE,1)-200];
    fcn = makeConstrainToRectFcn('imrect',...
         get(gca,'XLim'),get(gca,'YLim'));
    
    api.setPositionConstraintFcn(fcn);
    api.addNewPositionCallback(@(p) select_stepwedge(p));
    
    figure(step_fig); subplot(1,2,2);
    imagesc(IMAGE(pos(2):pos(2)+pos(4)-1, pos(1):pos(1)+pos(3)-1));
    axis image;
    
    region_ok = helpdlg('Is the region OK?', 'Select stepwedge');
    waitfor(region_ok);
    pos = round(api.getPosition());
    clear api;
    
elseif strcmpi(v, '2007a')
    rect_h = imrect(gca, [50 20 100 size(IMAGE,1)-40]);
    api = iptgetapi(rect_h);
    pos = [40 20 100 size(IMAGE,1)-40];
    fcn = makeConstrainToRectFcn('imrect',...
         get(gca,'XLim'),get(gca,'YLim'));
    
    api.setPositionConstraintFcn(fcn);
    api.addNewPositionCallback(@(p) select_stepwedge(p));
    
    figure(step_fig); subplot(1,2,2);
    imagesc(IMAGE(pos(2):pos(2)+pos(4)-1, pos(1):pos(1)+pos(3)-1));
    axis image;
    
    region_ok = helpdlg('Is the region OK?', 'Select stepwedge');
    waitfor(region_ok);
    pos = round(api.getPosition());
    clear api;
else
    %Put old method back here
end



% find the image size
%disp('select opposite corners of the step wedge area, then press return');
%dimensions = size(IMAGE);
% pointsgiven = 1;
% 
% while pointsgiven ~= 2
% 
%     [swx, swy, p] = impixel(IMAGE, []); %#ok 
%     dimsx = size(swx);
%     pointsgiven = dimsx(1);
%     
% end
% 
% % if box goes off edge of image then shift it in a bit:
% for ipoint=1:2
%     if(swx(ipoint) < 1)
%         swx(ipoint)=1;
%     elseif(swx(ipoint) > dimensions(2))
%         swx(ipoint)=dimensions(2);
%     end    
%     if(swy(ipoint) < 1)
%         swy(ipoint)=1;
%     elseif(swy(ipoint) > dimensions(1))
%         swy(ipoint)=dimensions(1);
%     end
% end
%sw_image = IMAGE(swy(1):swy(2),swx(1):swx(2));
sw_image = rot90(IMAGE(pos(2):pos(2)+pos(4)-1, pos(1):pos(1)+pos(3)-1));
[r_sw c_sw] = size(sw_image);


ok = -4;

while (ok<0)

    if(ok<=-3)
        
        set(step_fig, 'Name',...
            'Select steps');
        subplot(1,1,1);
        imagesc(sw_image); axis image; colormap(gray(256));
        title('Select the top and bottom of two steps, then press return (4 points in all)');
        set(gca, 'xticklabel', [], 'yticklabel', []);
        
        problem = 1;
        while problem
            [swx4, swy4, p] = impixel; %#ok
            problem = length(swx4) ~= 4;
            if problem
                new_roi = 'Select new stepwedge region';
                new_steps = 'Select new step points';
                answer = questdlg('Incorrect number of points selected. Please either:',...
                    'Stepwedge error', new_roi, new_steps, new_roi);

                if strcmpi(answer, new_roi);
                    close(step_fig);
                    [wedgevals swx swy errorcheck] = wedgecal_jen(IMAGE, debug_mode);
                    return;
                end
            end
        end
                
        
        %sort points by x-vals in case user clicked in the wrong order
        [swx4 idx] = sort(swx4);
        swy4 = swy4(idx);
        
        %plot the points so they stay on the image
        hold on;
        plot(swx4, swy4, 'r+', 'MarkerSize', 4);
        
        swx=zeros(2,1);
        swy=zeros(2,1);
    
        swx(1) = (swx4(1)+swx4(2))/2;
        swx(2) = (swx4(3)+swx4(4))/2;
        swy(1) = (swy4(1)+swy4(2))/2;
        swy(2) = (swy4(3)+swy4(4))/2;
    
        i_min = numsteps+1;
        i_max = 0;
        
        % get step numbers but try to guard against mistakes
        problem = 1;
        while problem
            
            prompt = {'What number is the first step?','What number is the second step?:'};
            dlg_title = 'Select step numbers';
            num_lines = 1;
            screen_size = get(0,'ScreenSize');
            options.Position = [screen_size(3)/4, screen_size(4)/2];
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
                i_min<1 || i_min>numsteps || i_max<1 || i_max>numsteps;
            
            if problem
                h = warndlg('Please choose again',...
                    'Invalid step numbers', 'modal');
                uiwait(h);
            end
            %i_min = input('what number is the first step? ');
            %i_max = input('what number is the second step? ');
            
        end
    
    elseif(ok==-1)
        i_min = i_min+1;
        i_max = i_max+1;
        
        figure(step_fig);
        imagesc(sw_image); axis image; colormap(gray(256)); axis off;
        
    elseif(ok==-2)
        i_min = i_min-1;
        i_max = i_max-1;
        
        figure(step_fig);
        imagesc(sw_image); axis image; colormap(gray(256)); axis off;    
    end
    
    numsteps_selected = i_max - i_min + 1;

    % draw a line down centre of wedge
    separation = (swx(2)-swx(1))/(numsteps_selected-1);
    offset = (swy(2)-swy(1))/(numsteps_selected-1);
    
    wedgevals = zeros(numsteps,1);
    wedgestd = zeros(numsteps,1);
    
    % draw rectangles on each step, representing area to be sampled for
    % grey level
    dodgy_std = 0;
    for i=1-i_min:numsteps-i_min
        
        figure(step_fig);
        drawrect([swx(1)+i*separation-3 swx(1)+i*separation+3],...
            [swy(1)+i*offset+5 swy(1)+i*offset-5]);
        
        % cjb mark a couple of steps with step number to aid debug!
        if((i+i_min)==10)
            text(swx(1)+i*separation,swy(1)+i*offset+5,'10');
        elseif((i+i_min)==20)
            text(swx(1)+i*separation,swy(1)+i*offset+5,'20');
        end
        % define indicies of this step's region
        xmin = round(swx(1)+i*separation-3);
        xmax = round(swx(1)+i*separation+3);
        ymin = round(swy(1)+i*offset-5);
        ymax = round(swy(1)+i*offset+5);
        
        % don't go outside image - this assumes y_max and x_max can't be
        % less than 1, and equivalently, that y_min and x_min can't be greater
        % than the image dims. However this could occur if a whole step
        % lies off the image
        if (ymin<0.5)
            ymin = 1;
        elseif(ymax>r_sw)
            ymax = r_sw;
        end
        if (xmin<0.5)
            xmin = 1;
        elseif(xmax>c_sw)
            xmax = c_sw;
        end
        this_step = sw_image(ymin:ymax,xmin:xmax);
        
        wedgevals(i+i_min) = mean(double(this_step(:)));
        wedgestd(i+i_min) = std(double(this_step(:)));
        
%       if(wedgestd(i+i_min) > 7.5) use 7.5 for 256 grey levels

        if(wedgestd(i+i_min) > 1000)
            if debug_mode
                disp('step has large std - look at intensity histogram');
                disp('press enter to continue...');
                pause;
                figure;
                imhist(this_step);
            else
                dodgy_std = 1;
            end
        end
        
    end
    if dodgy_std
        dstd = warndlg('One or more steps has large std');
        uiwait(dstd);
    end
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

if any(isnan(wedgevals))
    repeat_step = 'Repeat stepwedge selection';
    skip = 'Skip this mammogram';
    answer = questdlg('An error has occurred selecting step values. Please either:',...
        'Stepwedge error', repeat_step, skip, repeat_step);
    
    if strcmpi(answer, skip);
        errorcheck = 1;
        return;
    else
        [wedgevals swx swy errorcheck] = wedgecal_jen(IMAGE, debug_mode);
        return
    end
end

while(ok==2)
    
    figure(8);
    plot(wedgevals);
    
    bad_first = input('what is first bad step? ');
    bad_last = input('what is last bad step? ');
    
    if(bad_first == 1)
        wedgevals(bad_first:bad_last) = wedgevals(bad_last+1);
        wedgestd(bad_first:bad_last) = wedgestd(bad_last+1);
    elseif(bad_last == numsteps)
        wedgevals(bad_first:bad_last) = wedgevals(bad_first-1);
        wedgestd(bad_first:bad_last) = wedgestd(bad_first-1);
    else
        temp1 = 1:bad_last-bad_first+1;
        temp2 = temp1.*(wedgevals(bad_last+1)-wedgevals(bad_first-1))/(bad_last-bad_first+2) + wedgevals(bad_first-1);
        wedgevals(bad_first:bad_last) = temp2;    
    end
    figure(8);
    plot(wedgevals);
    
    % check if this is OK
    disp('Is this OK?');
    disp('enter "1" for yes');
    disp('enter "2" for no, try again');
    ok = input(':');
   
end


% first check there isn't a problem with the first step -
% if this goes too close to edge of film, then it may not be as dark as it should

if(wedgevals(1) < wedgevals(2))
    if debug_mode
        disp('problem with first step - not as dark as it should be');
    end
    wedgevals(1) = wedgevals(2);
    wedgestd(1) = wedgestd(2);
end

if debug_mode
    figure; hold on;
    plot(wedgevals,'bo');
    xlabel('step number');
    ylabel('mean intensity');
    title('stepwedge intensity vs step number');
end


% correct steps which have unphysical pixel value (eg. due to scatter/edge effects or stepwedge obstruction)
% i.e. as thickness of wedge increases, pixel value can only decrease.
% cjb commented this out because it does really really stupid things.
% 28/6/04
halfwaystep = round((numsteps-1)/2);
for i=halfwaystep:numsteps-1
    if(wedgevals(i)<wedgevals(i+1))
        if debug_mode
            disp(['step above ', num2str(i), ' has intensity higher - setting same']);
        end
        wedgevals(i+1)=wedgevals(i);
        wedgestd(i+1) = wedgestd(i);
    end
end
for i=halfwaystep:-1:2
    if(wedgevals(i-1)<wedgevals(i))
        if debug_mode
            disp(['step below ', num2str(i), ' has intensity higher - setting same']);
        end
        wedgevals(i-1)=wedgevals(i);
        wedgestd(i-1) = wedgestd(i);
    end
end

if debug_mode
    plot(wedgevals,'r+');
    legend('before correction','after correction');
    disp(wedgevals);
    
    for i=1:numsteps
    %     if(wedgestd(i) > 7.5) use 7.5 for 256 grey levels
        if(wedgestd(i) > 1000)
            disp('possible large std dodgyness on step:');
            disp(i);
            disp('press enter to continue...');
            pause
        end
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    function dlg = create_choice_dialog

        qstring = 'Are these regions ok?';
        choose_again = 'No, choose new pts';
        nudge_u = 'Nudge up a step';
        nudge_d = 'Nudge down a step';
        yes = 'Yes';
        
        screen_size = get(0,'ScreenSize');
        
        dlg = figure(...
            'Position',[screen_size(3)/4 - 165,screen_size(4)/2 - 45,330,95],...
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
    function select_stepwedge(pos) %#ok
        pos = round(pos);
        figure(step_fig); subplot(1,2,2);
        imagesc(IMAGE(pos(2):pos(2)+pos(4)-1, pos(1):pos(1)+pos(3)-1));
        axis image; axis off; colormap(gray(256));
        figure(region_ok);
    end
end





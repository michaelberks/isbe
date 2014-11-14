function [] = chitra_single_gui2()
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
color2 = [212 208 200]/255;
color3 = [0 0 0];
valid_ids = 1:99;
buff = 10;
screen_size = get(0,'ScreenSize');
centre_x = screen_size(3)/2;
%centre_y = screen_size(4)/2;
blob_path = 'C:\isbe\ug_project_options\2011\chitra\blobs\';
user_path = 'C:\isbe\ug_project_options\2011\chitra\users\user';
blob_files = dir([blob_path '*.mat']);

number_of_blobs = [];

%Create empty variables that exist globally and will be filled auxilliary
%functions - this set do not need to be reset for each session
ui = [];
axes_pos = [];
curr_blob = [];
instructions = [];

%main program runs in while loop - each loop runs one study session.
%Sessions will continue until the main fig is closed - however the user is
%blocked from doing this unless they have saved their data
continue_study = 1;
while true

    display('Start of study loop');

    %Assume we're going to show instructions and show the welcome screen
    do_instructions = 1;
    display_welcome_screen;
    uiwait(ui.main_fig);
    
    %User allowed to quit at this point - so check if continue study has 
    %been negated
    if ~continue_study
        break;
    end
    
    %Depending on whether the user ID entered in the welcome screen exists or
    %not, we'll either show the instructions and start a new study
    if do_instructions
        display_instruction_screen;
        uiwait(ui.main_fig);
        curr_blob = 1;
    else
        %Or else we'll have found and loading the existing user data - we then
        %need to find which mass to start from
        curr_blob = find(~user_data.blob_timings, 1);
    end

    %If we've loaded in data it could be all the blob_order are complete, in which
    %case curr_blob will be empty
    if ~isempty(curr_blob)

        %If the curr_blob isn't empty, lets go to the study
        display_observer_screen;
        uiwait(ui.main_fig);
    else
        % if curr_blob is empty we must have finished
        h = msgbox({'You already have a full set desnity readings for this part of the study.';... 
            'Please check you have entered the correct user ID'}, 'Study completed','help');
        uiwait(h);
    end

    %We're finished with this user, reset the user_data structure and
    %return to the start, 
    user_data = [];
    display_end_screen;
    display('End of study loop');
    uiwait(ui.main_fig);
    
    %User allowed to quit at this point - so check if continue study has 
    %been negated
    if ~continue_study
        break;
    end
end

display('End of function');
set(0,'DefaultFigureWindowStyle',orig_window_style);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Auxilliary functions that control data
%
%--------------------------------------------------------------------------
    function initialise_data

        %reset the rand seed based on user id
        rand('twister', 100 + user_data.user_id);
        
        %be careful of random seed here - must load a new one in or
        %we'll keep generating the same numbers!
        number_of_blobs = 10;
        
        %workout the missing images that must be included
        fid = fopen('C:\isbe\ug_project_options\2011\chitra\software\missing_single.txt');
        missing_txt = textscan(fid, '%s %s', 'delimiter', '=');
        fclose(fid);
        missing_ims = ...
            str2num(missing_txt{2}{strncmpi(missing_txt{1},...
            ['reader_' zerostr(user_data.user_id,2)], 9)});

        %now get a random selection of the other images to make up 10 in
        %total
        num_missing = length(missing_ims);
        extra_idx = setdiff(1:100, missing_ims);
        extra_ims = extra_idx(randperm(100 - num_missing));
        
        %Combine the missing and extra images and randomise the order
        all_ims = [missing_ims extra_ims]; 
        user_data.blob_order = all_ims(randperm(number_of_blobs))';
        user_data.blob_timings = zeros(number_of_blobs,1);
    end 

%--------------------------------------------------------------------------
    function save_user_data
        %save the users blob_order to a file labelled with their unique file
        %ID
        filename = [user_path, zerostr(user_data.user_id,2), '_1a.mat'];
        save(filename, 'user_data');
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Auxilliary functions so set up UI stuff (axes, buttons etc)
%
%--------------------------------------------------------------------------
    function display_welcome_screen
        %Generate main figure if it doens't already exist
        if ~isfield(ui, 'main_fig');
            ui.main_fig = figure(...
                'Position', [0 100 screen_size(3), screen_size(4)-100],...
                'Visible','on',...
                'Name', 'Synthetic mammogram density study: Reading Model I: Individual Assessment',...
                'NumberTitle', 'off',...
                'MenuBar', 'none',...
                'WindowStyle', 'normal',...
                'Color', color1,...
                'CloseRequestFcn', @quit_Callback);
        end
        
        instructions = {...
            ['Please enter your participant number in the box below',...
            ' then press continue.'];... 
            [];...
            ['If you have previously started a study',...
            ' this will load your existing density estimations and you will continue',...
            ' from the point at which you quit.'];...
            [];...
            'Otherwise a new session will be created.'};
        
        ui.instructions = uicontrol(...
            'Style','text',...
            'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3)/2, screen_size(4)/2],...
            'FontSize', 14,...
            'BackgroundColor', color1,...
            'Enable', 'on');
        
        [instructions new_pos] = textwrap(ui.instructions,instructions);
        set(ui.instructions, 'String', instructions, 'Position', new_pos);
        
        ui.user_id = uicontrol(...
            'Style','text',...
            'String', 'Participant number:',...
            'FontSize', 14,...
            'Position', [screen_size(3)/2 - 100, 3*buff+80, 200, 40],...
            'BackgroundColor', color1,...
            'Enable', 'on');
        
        ui.user_id_input = uicontrol(...
            'Style','edit',...
            'String',[],...
            'Position', [screen_size(3)/2 - 50, 2*buff+40, 100, 40],...
            'BackgroundColor', color1,...
            'Enable', 'on');
        
        ui.continue_button = uicontrol(...
            'Style','pushbutton',...
            'String','Continue',...
            'Tag','continue',...
            'Callback', @user_id_Callback,...
            'Position', [screen_size(3)/2 - 50, buff, 100, 40],...
            'BackgroundColor', color2,...
            'Enable', 'on');
        
        movegui(ui.main_fig, 'northwest');
    end

%--------------------------------------------------------------------------
    function display_instruction_screen
        figure(ui.main_fig);
        set(ui.main_fig, 'Color', color1);
        ui.continue_button = uicontrol(...
            'Style','pushbutton',...
            'String','Continue',...
            'Tag','continue',...
            'Callback', @continue_Callback,...
            'Position', [screen_size(3)/2 - 50, buff, 100, 40],...
            'BackgroundColor', color2,...
            'Enable', 'on');
        
        instructions = observer_instructions;
        
        ui.instructions = uicontrol(...
            'Style','text',...
            'Tag','file_open',...
            'BackgroundColor', color1,...
            'FontSize', 14,...
            'FontName', 'Arial',...
            'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3)/2, screen_size(4)/2]);
        
        [instructions new_pos] = textwrap(ui.instructions, instructions);
        set(ui.instructions,...
            'String', instructions,...
            'Position', new_pos);
    end
%--------------------------------------------------------------------------
    function display_end_screen
        figure(ui.main_fig);
        set(ui.main_fig, ...
            'Color', color1,...
            'CloseRequestFcn', @quit_Callback);
        ui.goodbye = uicontrol(...
            'Style','text',...
            'String', {'You are finished!'; []; 'Thank you for your time'},...
            'Tag','file_open',...
            'BackgroundColor', color1,...
            'FontSize', 16,...
            'FontName', 'Arial',...
            'Position', [screen_size(3)/2-150, screen_size(4)/2-150, 300, 300]);
        ui.continue_button = uicontrol(...
            'Style','pushbutton',...
            'String','Continue',...
            'Tag','continue',...
            'Callback', @continue_Callback,...
            'Position', [screen_size(3)/2 - 50, buff, 100, 40],...
            'BackgroundColor', color2,...
            'Enable', 'on');
    end

%--------------------------------------------------------------------------
    function display_observer_screen
        
        
        button_h = 40;
        button_w = 100;
        text_w = 250;
 
        axes_y = screen_size(4) - 100 - 3*buff - button_h;
        axes_x = 0.8 * axes_y;
        axes_pos = [centre_x - axes_x/2, 2*buff + button_h, axes_x, axes_y];
        
        
        figure(ui.main_fig);
        set(ui.main_fig,...
            'Color', color3,...
            'CloseRequestFcn', @no_quit_Callback);
                                   
        ui.next_image = uicontrol(...
            'Style','pushbutton',...
            'String','Next',...
            'Tag','pan_main',...
            'Callback', @next_image_Callback,...
            'Position', [centre_x - (button_w+buff-text_w)/2, buff, button_w, button_h],...
            'BackgroundColor', color1,...
            'ForegroundColor', color3,...
            'FontSize', 14,...
            'FontName', 'Arial',...
            'Enable', 'on');
        
        ui.regions_text = uicontrol(...
            'Style','text',...
            'Tag','file_open',...
            'BackgroundColor', color3,...
            'ForegroundColor', color1,...
            'FontSize', 14,...
            'FontName', 'Arial',...
            'Position', [centre_x - (text_w+button_w+buff)/2, buff, text_w, 0.8*button_h]);
                            
        ui.axes_main = axes(...
            'Parent', ui.main_fig,...
            'Units', 'pixels',...
            'Visible', 'on');
         
        set(ui.main_fig,...
            'Colormap', gray(256));
        ui.region = image([]);
        %ui.region = imagesc([]);
        
        set(ui.axes_main,...
            'Position', axes_pos,...
            'Xtick', [],...
            'Ytick', [],...
            'YDir','reverse',...
            'Xlim', [0.5 axes_pos(3)+0.5],...
            'Ylim', [0.5 axes_pos(4)+0.5]);
        
        update_observer_screen;

    end

% --------------------------------------------------------------------
    function update_observer_screen
    % 
        if curr_blob <= number_of_blobs

            %load in the new blob
            blob_mammogram = u_load([blob_path, blob_files(user_data.blob_order(curr_blob,1)).name]);
            
            %Compute the apsect ratios for these images (they may vary from
            %pair to pair)
            aspect_ratio = size(blob_mammogram,2) / size(blob_mammogram,1);
            axes_x = axes_pos(4)*aspect_ratio;
        
            %Update the size and position of the axes
            axes_pos(1) = centre_x - axes_x/2;
            axes_pos(3) = axes_x;
        
            set(ui.axes_main,...
                'Position', axes_pos,...
                'Xlim', [0.5 axes_pos(3)+0.5],...
                'Ylim', [0.5 axes_pos(4)+0.5]);
            
            %Display the new blob
            set(ui.region,...
                'Visible', 'on',...
                'CData', blob_mammogram,...
                'XData', [1 axes_pos(3)],...
                'YData', [1 axes_pos(4)]);

            %Start timer for recording viewing time of each mass
            tic;
            
            %Update the number of blobs text
            regions_text = ['Mammogram: ', num2str(curr_blob), ' of ', num2str(number_of_blobs)];
            set(ui.regions_text, 'String', regions_text);
        else
            %Delete the buttons and stuff from the figure
            delete(get(ui.main_fig, 'Children'));

            %Save the user ratings
            save_user_data;

            %Allow main function to continue
            uiresume(ui.main_fig);
        end

        
    end

%--------------------------------------------------------------------------
    function ins = observer_instructions
        ins = {...
    'Hello, and thank you for agreeing to take part in our study.';...
    [];...
    ['The aims of the project are (1) to assess the accuracy of '...
    'estimation of synthetic mammographic breast density and (2) to determine whether the accuracy of '...
    'detecting change in density is improved by comparing two synthetic mammograms side by side. '];...
    [];...
    ['You will be viewing cranio-caudal views of the breast. Light grey areas are to be interpreted as dense tissue, '...
    'and dark grey areas as fatty tissue.  Estimates will be based on percentage (%) density – that is '...
    'the percentage of the breast area that is occupied by dense tissue.'];...
    [];...
    [];...
    'Reading Model I : Individual Assessment';...
    [];...
    ['In this reading model you will view 100 synthetic mammograms one at a time.  You will be required to '...
    'estimate the density (%) of each mammogram. Please mark your answers using a biro on the visual '...
    'analogue scales that are provided. On the scale the left end point represents 0% (density) and the '...
    'right end point represents 100% (density). Please mark the scale as appropriate. The time taken for '...
    'you to view and assess every image will be recorded.']...
            };


    end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% UI Callbacks
%
%--------------------------------------------------------------------------
% --------------------------------------------------------------------
    function user_id_Callback(hObject, eventdata) %#ok
    % Callback to allow us to continue
        user_id = str2double(get(ui.user_id_input, 'String'));
        if ismember(user_id, valid_ids)
            filename = [user_path, zerostr(user_id,2), '_1a.mat'];
            if exist(filename, 'file')
                %Set do_instructions to 0 to load existing data and resume study
                do_instructions = 0;
                user_data = u_load(filename);
                number_of_blobs = size(user_data.blob_order, 1);
            else
                %Set do_instructions to 1 to display instructions, generate
                % a new random ordering of masses then start new study
                user_data.user_id = user_id;
                do_instructions = 1;
                initialise_data;
            end
            %delete all the objects from the figure
            delete(get(ui.main_fig, 'Children'));
            
            %Pass control back to the main function
            uiresume(ui.main_fig);
        else
            warndlg('Invalid participant ID, please enter again', 'Invalid ID', 'modal');
            set(ui.user_id_input, 'String', []);
        end   
            
    end

% --------------------------------------------------------------------
    function continue_Callback(hObject, eventdata) %#ok
    % Callback to allow us to continue
        %delete all the objects from the figure
        delete(get(ui.main_fig, 'Children'));
        uiresume(ui.main_fig);
    end
    
% --------------------------------------------------------------------
    function quit_Callback(hObject, eventdata) %#ok
    % Callback executed if the user tries to quit
        %Check if they're ok to quit
        continue_study = 0;
        delete(ui.main_fig);
    end

% --------------------------------------------------------------------
    function no_quit_Callback(hObject, eventdata) %#ok
    % Callback executed if the user tries to quit
        
        % Callback executed if the user tries to quit
        if curr_blob <= number_of_blobs

            %Study not completed - query whether they want to continue
            queststr = {'You have not estimated density for all the mammograms in the study.';...
                'Are you sure want to quit?'};
            continue_study = 'Continue the study';
            save_quit = 'Quit';

            answer = questdlg(queststr, 'Study not finished', continue_study, save_quit, continue_study);

            %if continue just exit this function
            if strcmpi(answer, continue_study);
                return;
            end
        end
        
        %If we're still quitting:
        
        %1 save the user data
        save_user_data
        
        %Clear the main figure
        delete(get(ui.main_fig, 'Children'));
        
        %3 Allow the main fig to be exited
        set(ui.main_fig, 'CloseRequestFcn', @quit_Callback);
        
        %4. Pass control back to the main function
        uiresume(ui.main_fig);
    end

% --------------------------------------------------------------------
    function previous_image_Callback(hObject, eventdata) %#ok
        curr_blob = curr_blob - 1;
        update_feedback_screen;
    end

% --------------------------------------------------------------------
    function next_image_Callback(hObject, eventdata) %#ok
        
        %record time taken to read image if less than a second assume the
        %user hit next by accident
        tt = toc;
        if tt < 1
            return;
        end
        
        %Otherwise save the time and continue to next image
        user_data.blob_timings(curr_blob) = tt;
        
        %Update blob count and screen
        curr_blob = curr_blob + 1;
        update_observer_screen;
    end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---------------------- END OF FUNCTION -----------------------------------
%--------------------------------------------------------------------------
end
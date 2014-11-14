function [] = chitra_study_gui()
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
blob_path = 'C:\isbe\ug_project_options\2011\chitra\blobs\';
user_path = 'C:\isbe\ug_project_options\2011\chitra\user';
mass_files = dir([blob_path '*.mat']);

number_of_masses = [];

%Create empty variables that exist globally and will be filled auxilliary
%functions - this set do not need to be reset for each session
ui = [];
axes_pos = [];
panel_pos = [];
curr_mass = [];
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
        display_instruction_screen('observer');
        uiwait(ui.main_fig);
        curr_mass = 1;
    else
        %Or else we'll have found and loading the existing user data - we then
        %need to find which mass to start from
        curr_mass = find(~user_data.ratings(:,2), 1);
    end

    %If we've loaded in data it could be all the ratings are complete, in which
    %case curr_mass will be empty
    if ~isempty(curr_mass)

        %If the curr_mass isn't empty, lets go to the study
        do_feedback = 0;
        display_observer_screen;
        uiwait(ui.main_fig);
    else
        % if curr_mass is empty ask the user if they want to review 
        % their ratings, but they don't have the option of repeating the
        % study
        answer = questdlg(...
            {'You already have a full set of ratings.';... 
            'Would you like to review your feedback of synthetic masses?'},...
            'Study completed', 'Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            do_feedback = 1;
        else
            do_feedback = 0;
        end
    end

    % If we're continuing to the feedback section (we won't do this if the user
    % has chosen to quit before completing the rating of each mass), display
    % the feedback instructions then start the feedback study
    if do_feedback;
        display_instruction_screen('feedback');
        uiwait(ui.main_fig);
        start_feedback_study;
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

        %be careful of random seed here - must load a new one in or
        %we'll keep generating the same numbers!
        number_of_masses = 60;
        
        %Set the random seed to be the same for consecutive pairs of user
        %IDs (i.e. 1,2 -> 1000, 3,4 ->2000,...). This keeps the permutation
        %for choosing synthetic masses of type A or B the same for each
        %pair. We then reverse whether we choose from A or B for each pair,
        %ensuring that each normal background appears equally as a type A
        %or B throughout the whole study
        rand('twister', ceil(user_data.user_id/2)*1000);
        
        %Select which of the 30 synthetic masses will come from each
        %method (ensuring 15 from each) 
        
        %Create random permutation of 1-30, for odd user_ids the first 15 
        %will be type A, 15-30 mass B. For even sessions, the first 15 will
        %be type B, 16-30 will be type A
        rp = randperm(30);
        
        %Get the set of masses to be used for this session
        %real masses have number 1-30
        %syn mass type A have number 31-60
        %syn mass type B have number 61-90
        if rem(user_data.user_id, 2)
            mass_id = [1:30, 30+rp(1:15), 60+rp(16:30)];
        else
            mass_id = [1:30, 30+rp(16:30), 60+rp(1:15)];
        end
        
        %now lets randomly order these masses, setting the random seed
        %based on the clock to generate different orders for consecutive
        %pairs of sessions
        rand('twister', sum(100*clock));
        mass_order = mass_id(randperm(number_of_masses))';
        
        %Generate random ordering of mass regions and pre-allocate ratings
        %array to store user responses
        user_data.ratings = zeros(number_of_masses, 3);
        user_data.ratings(:,1) = mass_order;
    end

%--------------------------------------------------------------------------
    function start_feedback_study
        
        %Get the masses identified as def or prob synthetic
        user_synthetic = find(user_data.ratings(:,2) <= 2);
        number_of_masses = length(user_synthetic);
        
        %If no masses, then inform user and return
        if ~number_of_masses
            h = warndlg([...
                'You rated no masses as definitely or probably synthetic '...
                'therefore we do not require any feedback'],...
                'No feedback required', 'modal');
            display('dialog box open');
            waitfor(h);
            return;
            
        end
        
        %Otherwise continue with feedback...
        
        %set the current mass back to 1
        curr_mass = 1;
        
        % Create the feedback structures if the don't aldreay exist (i.e.
        % we might be loading in a saved session)
        if ~isfield(user_data, 'feedback')
            user_data.feedback = zeros(number_of_masses, num_feedback_features+1);
            user_data.feedback(:,1) = user_data.ratings(user_synthetic',1);
            user_data.feedback_other = cell(number_of_masses,1);
        end
        
        %Enter feedback GUI screen
        display_feedback_screen;
        
        %Stop main function until main_fig resumes
        uiwait(ui.main_fig);
    end

%--------------------------------------------------------------------------
    function save_user_data
        %save the users ratings to a file labelled with their unique file
        %ID
        filename = [user_path, zerostr(user_data.user_id,2), '.mat'];
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
                'Name', 'Synthetic mass observer study',...
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
            ' this will load your existing answers and you will continue',...
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
    function display_instruction_screen(study_type)
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
        
        if strcmp(study_type, 'observer');
            instructions = observer_instructions;
        else
            instructions = feedback_instructions;
        end
        
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
    function display_feedback_screen
        
        x_max = screen_size(3);
        y_max = screen_size(4)-100;
        button_h = 40;
        button_w = 100;
        text_w = 250;
        y_bar = button_h + buff;
        y_menu = 20;
        axes_dim = y_max - y_bar - 2*buff - y_menu;
        centre = buff + axes_dim/2;
        
        axes_pos = [buff, y_bar+buff, axes_dim, axes_dim];
        panel_pos = [2*buff+axes_dim, y_bar+buff, x_max-axes_dim-3*buff, axes_dim];
        
        figure(ui.main_fig);
        set(ui.main_fig,...
            'Color', color3,...
            'CloseRequestFcn', @no_quit_Callback);
        
        %Re-create all the buttons and stuff (note old method of hiding
        %them then making them visible has been scrapped for the case that
        %we've loaded ina session straight into feedback mode and therefore
        %haven't already created them in study mode)
        ui.previous_image = uicontrol(...
            'Style','pushbutton',...
            'String','Previous',...
            'Tag','zoom_main',...
            'Callback', @previous_image_Callback,...
            'Position', [centre - button_w - text_w/2 - buff, buff, button_w, button_h],...
            'BackgroundColor', color2,...
            'Enable', 'off');
        ui.next_image = uicontrol(...
            'Style','pushbutton',...
            'String','Next',...
            'Tag','pan_main',...
            'Callback', @next_image_Callback,...
            'Position', [centre + text_w/2 + buff, buff, button_w, button_h],...
            'BackgroundColor', color2,...
            'Enable', 'on');
        ui.submit_answers = uicontrol(...
            'Style','pushbutton',...
            'String','Submit feedback',...
            'Tag','add_mass',...
            'Callback', @submit_feedback_Callback,...
            'Position', [panel_pos(1), buff, text_w, button_h],...
            'BackgroundColor', color2,...
            'Enable', 'on');
        
        ui.regions_text = uicontrol(...
            'Style','text',...
            'Tag','file_open',...
            'BackgroundColor', color3,...
            'ForegroundColor', color1,...
            'FontSize', 16,...
            'FontName', 'Arial',...
            'Position', [centre - text_w/2, buff, text_w, button_h]);
                            
        ui.axes_main = axes(...
            'Parent', ui.main_fig,...
            'Units', 'pixels',...
            'Visible', 'on');
        
        ui.panel = uipanel(...
            'Parent', ui.main_fig,...
            'Units', 'pixels',...
            'Position', panel_pos,...
            'Visible', 'on',...
            'BackgroundColor', color1);
        
        f_size = 0.8 / num_feedback_features;
        
        for ii = 1:num_feedback_features
            ui.feature(ii) = uicontrol(...
                'Style', 'Radio',...
                'String', feedback_features{ii},...
                'Units', 'normalized',...
                'Position',[0.05 0.2+(ii-1)*f_size 0.95 f_size],...
                'FontSize', 12,...
                'FontName', 'Arial',...
                'Parent',ui.panel,...
                'BackgroundColor', color1,...
                'Callback',@feedback_Callback);
        end
        ui.feature_other = uicontrol(...
            'Style', 'Edit',...
            'String', [],...
            'Units', 'normalized',...
            'Position',[0.025 0.025 0.95 0.15],...
            'FontSize', 9,...
            'FontName', 'Arial',...
            'Parent',ui.panel,...
            'BackgroundColor', color1,...
            'Min', 0,...
            'Max', 2,...
            'HorizontalAlignment', 'left',...
            'Callback',@feedback_Callback);
        ui.additional_comments = uicontrol(...
            'Style', 'Text',...
            'String', 'Any other comments:',...
            'Units', 'normalized',...
            'Position',[0.025 0.176 0.95 0.022],...
            'FontSize', 9,...
            'FontName', 'Arial',...
            'Parent',ui.panel,...
            'BackgroundColor', color1,...
            'HorizontalAlignment', 'left',...
            'Callback',@feedback_Callback);
        
        ui.region = image([]);
        %ui.region = imagesc([]);
        
        set(ui.axes_main,...
            'Position', axes_pos,...
            'Xtick', [],...
            'Ytick', [],...
            'YDir','reverse',...
            'Xlim', [0.5 axes_pos(3)+0.5],...
            'Ylim', [0.5 axes_pos(4)+0.5]);
        
        set(ui.main_fig,...
            'Colormap', gray(256));
        
        update_feedback_screen;
        
    end

%--------------------------------------------------------------------------
    function display_observer_screen
        
        x_max = screen_size(3);
        y_max = screen_size(4)-100;
        button_h = 40;
        button_w = 100;
        text_w = 250;
        y_bar = button_h + buff;
        y_menu = 20;
        axes_dim = y_max - y_bar - 2*buff - y_menu;
        centre = buff + axes_dim/2;
        
        figure(ui.main_fig);
        set(ui.main_fig,...
            'Color', color3,...
            'CloseRequestFcn', @no_quit_Callback);
        
        %bar_pos = [0 0 x_max, y_bar];
        axes_pos = [buff, y_bar+buff, axes_dim, axes_dim];
        panel_pos = [2*buff+axes_dim, y_bar+buff, x_max-axes_dim-3*buff, axes_dim];
                                   
        ui.next_image = uicontrol(...
            'Style','pushbutton',...
            'String','Next',...
            'Tag','pan_main',...
            'Callback', @next_image_Callback,...
            'Position', [panel_pos(1), buff, button_w, button_h],...
            'BackgroundColor', color2,...
            'Enable', 'on');
        ui.submit_answers = uicontrol(...
            'Style','pushbutton',...
            'String','Submit ratings',...
            'Tag','add_mass',...
            'Callback', @submit_ratings_Callback,...
            'Position', [panel_pos(1)+buff+button_w, buff, button_w, button_h],...
            'BackgroundColor', color2,...
            'Enable', 'on');
        
        ui.regions_text = uicontrol(...
            'Style','text',...
            'Tag','file_open',...
            'BackgroundColor', color3,...
            'ForegroundColor', color1,...
            'FontSize', 16,...
            'FontName', 'Arial',...
            'Position', [centre - text_w/2, buff, text_w, button_h]);
                            
        ui.axes_main = axes(...
            'Parent', ui.main_fig,...
            'Units', 'pixels',...
            'Visible', 'on');
        
        ui.panel = uibuttongroup(...
            'Parent', ui.main_fig,...
            'Units', 'pixels',...
            'Position', panel_pos,...
            'Visible', 'on',...
            'BackgroundColor', color1,...
            'SelectionChangeFcn', @rating_Callback);
        
        ui.def_s = uicontrol(...
            'Style', 'Radio',...
            'String', 'Definitely synthetic',...
            'Units', 'normalized',...
            'Position',[0.05 0.8 0.9 0.2],...
            'FontSize', 16,...
            'FontName', 'Arial',...
            'Parent',ui.panel,...
            'BackgroundColor', color1,...
            'HandleVisibility','off');
        ui.pos_s = uicontrol(...
            'Style', 'Radio',...
            'String', 'Probably synthetic',...
            'Units', 'normalized',...
            'Position',[0.05 0.6 0.9 0.2],...
            'FontSize', 16,...
            'FontName', 'Arial',...
            'Parent',ui.panel,...
            'BackgroundColor', color1,...
            'HandleVisibility','off');
        ui.unsure = uicontrol(...
            'Style', 'Radio',...
            'String', 'Possibly real',...
            'Units', 'normalized',...
            'Position',[0.05 0.4 0.9 0.2],...
            'FontSize', 16,...
            'FontName', 'Arial',...
            'Parent',ui.panel,...
            'BackgroundColor', color1,...
            'HandleVisibility','off');
        ui.pos_r = uicontrol(...
            'Style', 'Radio',...
            'String', 'Probably real',...
            'Units', 'normalized',...
            'Position',[0.05 0.2 0.9 0.2],...
            'FontSize', 16,...
            'FontName', 'Arial',...
            'Parent',ui.panel,...
            'BackgroundColor', color1,...
            'HandleVisibility','off');
        ui.def_r = uicontrol(...
            'Style', 'Radio',...
            'String', 'Definitely real',...
            'Units', 'normalized',...
            'Position',[0.05 0.0 0.9 0.2],...
            'FontSize', 16,...
            'FontName', 'Arial',...
            'Parent',ui.panel,...
            'BackgroundColor', color1,...
            'HandleVisibility','off');
        
        set(ui.panel, 'SelectedObject', []);
         
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

        %load in the new mass and display
        mass_ROI = u_load([blob_path, mass_files(user_data.ratings(curr_mass,1)).name]);
        set(ui.region,...
            'Visible', 'on',...
            'CData', mass_ROI,...
            'XData', [1 axes_pos(3)],...
            'YData', [1 axes_pos(4)]);
        
        %Main study:
        % Always the next mass to never have a selected object
        set(ui.panel, 'SelectedObject', []);
        set(ui.next_image, 'Enable', 'off');

        %Start timer for recording viewing time of each mass
        tic;               

        %Update the number of regions text
        regions_text = ['Mass region: ', num2str(curr_mass), ' of ', num2str(number_of_masses)];
        set(ui.regions_text, 'String', regions_text);
    end

% --------------------------------------------------------------------
    function update_feedback_screen
    %           
        %load in the new mass and display
        mass_ROI = u_load([blob_path, mass_files(user_data.feedback(curr_mass,1)).name]);
        set(ui.region,...
            'Visible', 'on',...
            'CData', mass_ROI,...
            'XData', [1 axes_pos(3)],...
            'YData', [1 axes_pos(4)]);
        
        %Update the number of regions text
        regions_text = ['Mass region: ', num2str(curr_mass), ' of ', num2str(number_of_masses)];
        set(ui.regions_text, 'String', regions_text);
        
        %Enable/disable next/previous as appropriate
        if curr_mass == number_of_masses
            set(ui.next_image, 'Enable', 'off');
        else
            set(ui.next_image, 'Enable', 'on');
        end
        if curr_mass == 1
            set(ui.previous_image, 'Enable', 'off');
        else
            set(ui.previous_image, 'Enable', 'on');
        end

        %Feedback study:
        % Check if a feedback exists for this mass a check radio buttons
        % accordingly
        for ii = 1:num_feedback_features;
            set(ui.feature(ii), 'Value', user_data.feedback(curr_mass, ii+1));
        end
        set(ui.feature_other, 'String', user_data.feedback_other{curr_mass});
        
    end

%--------------------------------------------------------------------------
    function ins = observer_instructions
        ins = {...
    'Hello, and thank you for agreeing to take part in our study.';...
    [];...
    ['You will now be shown a randomly ordered sequence of mammogram regions,'...
    ' each of which contains a mass lesion. Approximately half of the masses',...
    ' are from real digitised mammograms. The other half have been synthetically generated.'];...
    [];...
    ['Please rate each mass as one of ''Definitely synthetic'', ''Probably synthetic'',',...
    ' ''Possibly real'', ''Probably real'' or ''Definitely real'' using the',...
    ' buttons displayed to the right of each mass.'];...
    [];...
    ['When you are happy with your rating, click ''Next'' to move to the',...
    ' next mass in the sequence. After rating the last mass, please click ''Submit ratings''.'];...
    [];...
    ['If necessary, you may click ''Submit ratings'' at any point during the study.',...
    ' This will save your current ratings and allow you to complete the study',...
    ' at a later time. However, we would much prefer you to complete the study in one session.'];...
            };

    end

%--------------------------------------------------------------------------
    function ins = feedback_instructions
        ins = {...
    ['You have finished rating the masses, thank you!!',...
    ' We would now like you to provide further information on how you identified synthetic masses.'];...
    [];...
    ['You will be shown each of the masses that you rated as either definitely',...
    ' or probably synthetic. In each case, use the buttons on the right of the',...
    ' mass to select any features of the mass that caused you to rate it as synthetic.',...
    ' You may select as many of the features as you like. Additionally, if you have any',...
    ' other comments on the appearance of the mass, please type these in the box provided.'];...
    [];...
    ['If, on review, you think that the mass is not synthetic, there is an option',...
    ' to select this. Note however, this will not change your original rating.'];...
    [];...
    ['You may use the ''Previous'' and ‘Next'' buttons to move between the masses as you wish.',...
    ' When you have finished, please click ''Submit feedback''. You will then have completed the study.'];...
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
            filename = [user_path, zerostr(user_id,2), '.mat'];
            if exist(filename, 'file')
                %Set do_instructions to 0 to load existing data and resume study
                do_instructions = 0;
                l = load(filename);
                %this will load a structure called user_data - avoid u_load
                %here as minimising non-standard functions to be compiled
                user_data = l.user_data; clear l;
                number_of_masses = size(user_data.ratings, 1);
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
        
        %Check if they're ok to quit
        h = warndlg(...
            'Quit has been disabled. To exit the study please subtmit your answers',...
            'Study not finished', 'modal');
        waitfor(h);
    end

% --------------------------------------------------------------------
    function previous_image_Callback(hObject, eventdata) %#ok
    % Callback to the "zoom" button allowing the user to zoom in/out the 
    % mammogram
        curr_mass = curr_mass - 1;
        update_feedback_screen;
    end

% --------------------------------------------------------------------
    function next_image_Callback(hObject, eventdata) %#ok
    % Callback to the "pan" button allowing the user to pan around the 
    % mammogram
        curr_mass = curr_mass + 1;
        if do_feedback
            update_feedback_screen;
        else
            update_observer_screen
        end
    end

% --------------------------------------------------------------------
    function submit_ratings_Callback(hObject, eventdata) %#ok
    % Callback to the "select nipple" allowing the user to mark the 
    % position of the nipple on the mammogram
    
        if any(~user_data.ratings(:,2))
            %Study not completed - query whether they want to continue
            queststr = 'You have not rated all of the masses in the study. Please either: ';
            continue_study = 'Continue the study';
            save_quit = 'Save answers and quit';

            answer = questdlg(queststr, 'Study not finished', continue_study, save_quit, continue_study);
            
            %if continue just exit this function
            if strcmpi(answer, continue_study);
                return;
            end
            %if they choose to quit now don't show the feedback instructions
            % do_feedback already set to zero 
        else
            %Study completed, show feedback instructions
            do_feedback = 1;
        end

        %Delete the buttons and stuff from the figure
        delete(get(ui.main_fig, 'Children'));
        
        %Save the user ratings
        save_user_data;
        
        %Allow main function to continue
        uiresume(ui.main_fig);
        
    end

% --------------------------------------------------------------------
    function submit_feedback_Callback(hObject, eventdata) %#ok
    % Callback to the "select nipple" allowing the user to mark the 
    % position of the nipple on the mammogram
        if any(all(~user_data.feedback(:,2:end),2) & cellfun(@isempty, user_data.feedback_other))
            %Study not completed - query whether they want to continue
            queststr = 'You have not given feedback for all of the masses rated as synthetic. Please either: ';
            continue_study = 'Continue the study';
            save_quit = 'Save answers and quit';

            answer = questdlg(queststr, 'Feedback not finished', continue_study, save_quit, continue_study);
            
            %if continue just exit this function
            if strcmpi(answer, continue_study);
                return;
            end
        end
        
        %Save the feedback
        save_user_data;
        
        %Delete the buttons, text etc. from figure
        delete(get(ui.main_fig, 'Children'));        
        
        %Save the user ratings
        save_user_data;
        
        %Allow main function to continue
        uiresume(ui.main_fig);
    end
% --------------------------------------------------------------------
    function rating_Callback(hObject, eventdata) %#ok
        
        %record time taken to make choice (from when mass first displayed)
        user_data.ratings(curr_mass, 3) = round(toc);
        
        %Record rating
        rating = get(ui.panel, 'SelectedObject');
        switch rating
            case ui.def_s
                user_data.ratings(curr_mass, 2) = 1;
            case ui.pos_s
                user_data.ratings(curr_mass, 2) = 2;
            case ui.unsure
                user_data.ratings(curr_mass, 2) = 3;
            case ui.pos_r
                user_data.ratings(curr_mass, 2) = 4;
            case ui.def_r
                user_data.ratings(curr_mass, 2) = 5;
        end
        
        %Enable the next image button
        if curr_mass < number_of_masses
            set(ui.next_image, 'Enable', 'on');
        end
    end
% --------------------------------------------------------------------
    function feedback_Callback(hObject, eventdata) %#ok
    %   
    %   If the callback due to a comment being entered
        if hObject == ui.feature_other;
            %Save comment in feedback other array
            user_data.feedback_other{curr_mass} = get(ui.feature_other, 'String');
        
        %If the callback was from the no longer synthetic button
        elseif hObject == ui.feature(1);
            not_syn = get(ui.feature(1), 'Value');
            user_data.feedback(curr_mass, 2) = not_syn;
            
            %If not synthetic
            if not_syn
                %Set all other buttons and associated feedback to 0
                for ii = 2:num_feedback_features;
                    user_data.feedback(curr_mass, ii+1) = 0;
                    set(ui.feature(ii), 'Value', 0);
                end
            end
        else
            %Callback for a synthetic feature - work out which feature    
            for ii = 2:num_feedback_features;
                if hObject == ui.feature(ii);
                    user_data.feedback(curr_mass, ii+1) = get(ui.feature(ii), 'Value');
                end
            end
            %Set the not synthetic button and feedback to 0
            user_data.feedback(curr_mass, 2) = 0;
            set(ui.feature(1), 'Value', 0);
        end
        
    end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---------------------- END OF FUNCTION -----------------------------------
%--------------------------------------------------------------------------
end
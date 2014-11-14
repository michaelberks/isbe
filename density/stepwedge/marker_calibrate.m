function marker_calibrate(varargin)

%Now using u_parkargs interface
args = u_packargs(varargin, 0, ...
            'DebugMode', 0, ...
            'User', [] ...
			);

%Swicth on/off figures to help check method is working
debug_mode = args.DebugMode;

format bank; % don't scale by 1000
iptsetpref('ImshowBorder','tight'); % show less border round image
warning off Images:initSize:adjustingMag; % don't clutter up output with image size warnings
warning off Images:imshow:magnificationMustBeFitForDockedFigure
warning off MATLAB:Figure:SetPosition

orig_window_style = get(0,'DefaultFigureWindowStyle');
if ~strcmp(orig_window_style, 'docked')
    display('Warning: changing window style to normal for function');
    set(0,'DefaultFigureWindowStyle','docked');
end

%select files to process
%disp('Select the files you wish to process');
[filenames, image_path] = uigetfile('E:\*.tif','Select the mammograms to process','Multiselect','on');

if isnumeric(filenames) && ~filenames
    h = warndlg('No files selected', 'You must select at least one calibration image to process');
    uiwait(h);
    marker_calibrate(varargin);
    return;
end

%If user selects only 1 file need to make sure name is saved in a cell
%array
if ~iscell(filenames)
        temp = filenames;
        filenames = cell(1);
        filenames{1} = temp; clear temp;
end

%select folder to save user data files to
%disp('Select the directory you wish to save the user input to');
data_path = uigetdir(image_path, 'Select the directory to save the user input to');

% loop over the selected mammograms



van_name = filenames{1}(1:4);

calibration_data_filename = [data_path, filesep, van_name, '_calibration_data.mat'];
%--------------------------------------------------------------------------
% Set up data structures for storing calibration data
%--------------------------------------------------------------------------
if exist(calibration_data_filename, 'file');
    calibration_data = u_load(calibration_data_filename); 
    if ~isfield(calibration_data, 'marker_dists');
        calibration_data(2).marker_dists = [];
    end
else
    calibration_data(2).marker_dists = [];
end

%--------------------------------------------------------------------------
%Pre-allocate a structure to save failed cases to
failures = [];
n_process = length(filenames);
max_pairs = [3 4];
%--------------------------------------------------------------------------
%Process each image
%--------------------------------------------------------------------------
for i_file=1:n_process
    
    %Generate a filename to save the user input to
    user_data_filename = [data_path, filesep, filenames{i_file}(1:end-4), '_data.mat'];
    
    %Check whether data already exists for this file
    if exist(user_data_filename, 'file');
        
        answer = questdlg(...
                {'Data already exists for this calibration image, continue to next image?';...
                 'Warning, selecting ''No'' will overwrite the old data'},'Mammogram already processed','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        end
    end
    disp(['processing: ' image_path, filenames{i_file}]);
    
    h = waitbar(0,'Loading and processing image. Please wait...');
    
    %----------------------------------------------------------------------
    % workout film size - 18x24 and 24x30 are digitised at different orientations
    % 18x24 films will need to be rotated
    if ~isempty(strfind(filenames{i_file}, '1824'))
        filmsizes = 1;
    elseif ~isempty(strfind(filenames{i_file}, '2430'))
        filmsizes = 2;
    else
        %define error action
    end
    
    % read in the image
    [IMAGE] = imread([image_path, filenames{i_file}],'tif');
    
    if filmsizes == 1
        % this is only needed for the 1824 film sizes
        IMAGE = rot90(IMAGE,1);
    end

    % Filter image down size and reduce bit-depth
    % median filter
    IMAGE = medfilt2(IMAGE);    
    close(h);
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % call markerdetect to get thickness image
    % use 'filmsizes' to choose magnification calibration constants
    [x_marker,y_marker,r_marker,selected_markers,marker_ui_xy,errorcheck] = markerdetect(IMAGE, max_pairs(filmsizes), 0, debug_mode);
    user_data.errorcheck = errorcheck;
    user_data.user = args.User;
    user_data.user_auto = getenv('UserName');
    save(user_data_filename, 'user_data');
    
    %Check we have marker pair data to continue with, otherwise skip this
    %mammogram
    if errorcheck
        disp('No marker pair data. Skipping this mammogram');
        failures{end+1,1} = filenames{i_file}(1:end-4); %#ok
        answer = questdlg(...
            'Continue to next image','Mammogram Complete','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        else
            break;
        end
    end
    
    %----------------------------------------------------------------------
    %If no error, compute distances between each marker pair
    %----------------------------------------------------------------------
    distances = sqrt(...
        (x_marker(1,:)-x_marker(2,:)).^2 +...
        (y_marker(1,:)-y_marker(2,:)).^2);  
    
    %----------------------------------------------------------------------
    % Get thickness from the filename
    t_idx = strfind(filenames{i_file},'cm');
    t_error = true;
    if ~isempty(t_idx)
        thickness = str2double(filenames{i_file}(t_idx-1));
        if ~isnan(thickness)
            t_error = false;
        end
    end
    
    if ~t_error
        answer = questdlg(...
                {['Thickness for file ', filenames{i_file}];...
                 '';...
                 ['Thickness = ', num2str(thickness) 'cm'];...
                 '';...
                 'Is this correct?'},'Thickness data','Yes', 'No', 'Yes');
        if strcmpi(answer, 'no')
            t_error = true;
        end
    end
    if t_error
        %manually enter thickness data
        prompt = {'What is the thickness for this image?'};
        dlg_title = 'Select thicknesses';
        num_lines = 1;
        screen_size = get(0,'ScreenSize');
        options.Position = [screen_size(3)/4, 3*screen_size(4)/5];
        answer = inputdlg(prompt,dlg_title,num_lines,{''}, options);
        thickness = str2num(answer{1}); %#ok
            
    end
    %Convert thickness from cm to mm
    thickness = 10*thickness;
    
    %----------------------------------------------------------------------
    % Save calibration data - do this after each image so we don't lose
    % existing data in a crash
         
    %Check whether this thicknesses already exists in the calibration
    %data
    if ~isempty(calibration_data(filmsizes).marker_dists)
        data_exists = calibration_data(filmsizes).marker_dists(:,1) == thickness;
    else
        data_exists = false;
    end
    if any(data_exists)
        %If so, overwrite the marker distance data
        calibration_data(filmsizes).marker_dists(data_exists, selected_markers+1) = distances;
    else
        %Otherwise append to end
        calibration_data(filmsizes).marker_dists(end+1, 1) = thickness; %#ok
        calibration_data(filmsizes).marker_dists(end, selected_markers+1) = distances;
    end
    save(calibration_data_filename, 'calibration_data');

    answer = questdlg(...
            'Continue to next image','Mammogram Complete','Yes', 'No', 'Yes');
    
    if strcmpi(answer, 'no')
        break;
    end
     
end
%--------------------------------------------------------------------------
% Now we've processed all the calibration images get marker pair distance
% relationships for each filmsize

marker_fig = figure;
size_labels = {'18x24','24x30'};
colours = 'rgby';

for filmsizes = 1:2
    subplot(1,2,filmsizes); hold on;
    if ~isempty(calibration_data(filmsizes).marker_dists)
        calibration_data(filmsizes).M_mag = zeros(1,max_pairs(filmsizes));
        calibration_data(filmsizes).C_mag = zeros(1,max_pairs(filmsizes));

        %For each marker pair get relationship between distance and thickness
        for ii = 1:max_pairs(filmsizes)
            %For 18x24 films
            %Only using markers we have distance data for
            idx = calibration_data(filmsizes).marker_dists(:,ii+1) > 0;

            %Fit a polynomial of thickness against distance
            p = polyfit(...
                calibration_data(filmsizes).marker_dists(idx,1),...
                calibration_data(filmsizes).marker_dists(idx,ii+1),1);

            %Record the results in calibration data
            calibration_data(filmsizes).M_mag(ii) = p(1);
            calibration_data(filmsizes).C_mag(ii) = p(2);
            
            plots(ii) = plot([20 100], p(1)*[20 100]+p(2), [colours(ii), ':'], 'LineWidth', 2); %#ok
            labels{ii} = ['Pair ', num2str(ii)]; %#ok
            plot(...
                calibration_data(filmsizes).marker_dists(idx,1),...
                calibration_data(filmsizes).marker_dists(idx,ii+1), [colours(ii), 'x'], 'MarkerSize', 10);
            hold on;
            
            xlabel('Breast thickness (mm)');
            ylabel('Marker distance (pixels)');
        end
        title([van_name, ' marker calibration, ', size_labels{filmsizes}, ' films']);
        legend(plots, labels, 'location', 'northwest');
        clear plots labels
    end
end

%Save final calibration data
save(calibration_data_filename, 'calibration_data');

h = helpdlg('Finished processing all selected mammograms','Finished!');
disp('Finished processing all selected mammograms');
uiwait(h);
if exist('marker_fig', 'var')
    close(marker_fig);
end

%turn warnings back on
warning on Images:initSize:adjustingMag;
warning on Images:imshow:magnificationMustBeFitForDockedFigure
warning on MATLAB:Figure:SetPosition
set(0,'DefaultFigureWindowStyle',orig_window_style);

%--------------------------------------------------------------------------
%------------------------ END OF FUNCTION ---------------------------------
%--------------------------------------------------------------------------

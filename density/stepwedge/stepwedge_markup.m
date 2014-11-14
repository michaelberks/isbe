function stepwedge_markup(varargin)

%Now using u_parkargs interface
args = u_packargs(varargin, 0, ...
			'M_mag', [], ...
			'C_mag', [], ...
            'ResizeFactor', 0.176, ...
            'DebugMode', 0, ...
            'MaxPairs', [],...
            'User', [] ...
			);

% this is the amount by which resolution is reduced after markers have been found (44um to 250um)        
resize_factor = args.ResizeFactor;

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

%  CALIBRATION DATA  
disp('***************** CALIBRATION DATA **************************');

%  Marker calibration  
%  x_pixels = M_mag * x_b + C_mag  

if isempty(args.M_mag)
    % use values from 2008 digitisation of 2008 calibration images
    % [18x24 pair 1 - 4; 24x30 pair 1 - 4]
    % format is M_mag(i_sw,i_point)
    M_mag = [8.10, 8.35, 8.35, 8.36; 10.32, 10.33, 10.41, 10.42]; % units = pixels per mm
    disp('Using marker calibration from 2008 digitisation');
else
    M_mag = args.M_mag;
end
if isempty(args.C_mag)
    C_mag = [4475.2, 4472.2, 4477.7, 4479.5; 5515.3, 5517.5, 5522.6, 5527.5]; % units = pixels
else
    C_mag = args.C_mag;
end

% Values for Bury daylight processor films digitised using CadPro
% The relationship is not straightforward so a Matlab interpolation
% function has been used, rather than using an equation to define the
% relationship (as was done previously)
disp('Using stepwedge calibration from 2008 digitisation of Bury films');


% kV constants were previously applied but it was found that calibration
% was independent of kV (Tom Marchant MSc)

disp('***************** CALIBRATION DATA **************************');

%select files to process
%disp('Select the files you wish to process');
[filenames, image_path] = uigetfile('E:\*.tif','Select the mammograms to process','Multiselect','on');

%select folder to save thicknesses/sw locations to
%disp('Select the directory you wish to save the user input to');
data_path = uigetdir([], 'Select the directory to save the user input to');
  

% loop over the selected mammograms

%If user selects only 1 file need to make sure name is saved in a cell
%array
if ~iscell(filenames)
    temp = filenames;
    filenames = cell(1);
    filenames{1} = temp; clear temp;
end

n_process = length(filenames);
%%
for i_file=1:n_process
    
    %Generate a filename to save the user input to
    data_filename = [data_path, filesep, filenames{i_file}(1:end-4), '_data.mat'];
    
    %Check whether data already exists for this file
    if exist(data_filename, 'file');
        
        answer = questdlg(...
                {'Data already exists for this mammogram, continue to next image?';...
                 'Warning, selecting ''No'' will overwrite the old data'},'Mammogram already processed','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        end
    end
    disp(['processing: ' image_path, filenames{i_file}]);
    
    %ikv=1; %ikv = which_kv(i_file); % no longer using kv dependent constants    
    
    % define film size
    % 18x24 and 24x30 are digitised at different orientations
    % 18x24 films will need to be flipped
    if ~isempty(strfind(filenames{i_file}, '1824'))
        filmsizes = 1;
        max_pairs = 3;
    elseif ~isempty(strfind(filenames{i_file}, '2430'))
        filmsizes = 2;
        max_pairs = 4;
    else
        %define error action
    end
    %if max_pairs has been set manually, override the automatic choice
    if ~isempty(args.MaxPairs)
        max_pairs = args.MaxPairs;
    end
    
    %is this a left or right breast (search for R not L as L in MLO!)
    left_breast = isempty(strfind(filenames{i_file}, 'R'));
    
    %is this an MLO or a CC
    mlo = ~isempty(strfind(filenames{i_file}, 'ML'));
    
    % read in the image
    [IMAGE] = imread([image_path, filenames{i_file}],'tif');
    
    if filmsizes == 1
        % this is only needed for the 1824 film sizes
        IMAGE = imrotate(IMAGE,90);
    end

    % What is the justification for median filtering at this point?!
    % median filter
    IMAGE = medfilt2(IMAGE);

    % call markerdetect to get thickness image
    % use 'filmsizes' to choose magnification calibration constants
    [x_marker,y_marker,r_marker,selected_markers,errorcheck] = markerdetect(IMAGE, max_pairs, left_breast, debug_mode);
    density_data.errorcheck = errorcheck;
    density_data.user = args.User;
    density_data.user_auto = getenv('UserName');
    save(data_filename, 'density_data');
    
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
    
    [x_b_info] = thickness_from_markers(...
        x_marker, y_marker, ...
        M_mag(filmsizes,selected_markers), C_mag(filmsizes,selected_markers));
        
    density_data.x_b_info = x_b_info;
    density_data.x_marker = x_marker;
    density_data.y_marker = y_marker;
    density_data.r_marker = r_marker;
    save(data_filename, 'density_data');

    % reduce image sizes here (don't need such high resolution now magnification markers have been located)
    IMAGE = imresize(IMAGE,resize_factor);

    %mark and save the position of the nipple in case this is used later
    [nipple_position] = nipple_select(IMAGE);
    density_data.nipple_position = nipple_position;
    save(data_filename, 'density_data');   

    % Now reduce the image to 4096 grey levels and invert to make it compatible
    % with sw_lookup_table.m
    IMAGE = IMAGE./16;
    IMAGE = 4095-IMAGE;

    % next need to locate stepwedge (do wedgevals calculations here for back
    % compatibility, though strictly it would be nicer not to - it's a
    % calculation not user input)
    [step_centres_xy] = locate_stepwedge_manual(IMAGE, filmsizes, debug_mode);
    [wedgevals, errorcheck] = calculate_step_values(IMAGE, step_centres_xy, debug_mode);
    
    %save outputs
    density_data.errorcheck = errorcheck;
    density_data.wedgevals = wedgevals;
    density_data.step_centres_xy = step_centres_xy;
    save(data_filename, 'density_data');
    
    %If MLO we need to mark the pectoral muscle
    if mlo
        %if left breast rotate image to be upright
        if left_breast
            IMAGE = rot90(IMAGE, 2);
        end
        [pectoral_mask pectoral_x pectoral_y] = pectoral_user(IMAGE);
        
        %Save the pectoral masks (less memory to save xy coord of region)
        density_data.pectoral_x = pectoral_x;
        density_data.pectoral_y = pectoral_y;
        save(data_filename, 'density_data');
    end
    
    %Clear density_data and ready for next mammogram
    clear density_data;
    
    answer = questdlg(...
            'Continue to next image','Mammogram Complete','Yes', 'No', 'Yes');
    
    if strcmpi(answer, 'no')
        break;
    end
     
end
helpdlg('Finished processing all selected mammograms','Finished!')
disp('Finished processing all selected mammograms');

%turn warnings back on
warning on Images:initSize:adjustingMag;
warning on Images:imshow:magnificationMustBeFitForDockedFigure
warning on MATLAB:Figure:SetPosition
set(0,'DefaultFigureWindowStyle',orig_window_style);

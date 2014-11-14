function stepwedge(varargin)
%STEPWEDGE Main user interface for marking up images to compute breast
%density using the Manchetser stepwedge method
%
% STEPWEDGE uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Required Arguments: none
%
% Optional Arguments:
%
% 'MarkerCalibration ([])'
%   - Can be used to input specific calibration data, otherwise (i.e. []) the user
%   will be asked to select file containing the calibration data
%
% 'ResizeFactor (0.176)'
%   - Factor by which image will be reduced once the breast thickness
%   markers have been detected
%
% 'DebugMode (0)'
%   - Controls what images are displayed during execution. Can be set to 1
%   or 2.
%
% 'ThicknessMethod (1)' 
%   - Determine which method for computing breast thickness across the breast
%     once markers have been determined. Has the following options
%   (1) Fit straight line through thickness at each marker
%   (2) Linearly interpolate thicknesses between markers
%
% 'TaperMethod (1)' 
%   - Determine which method to use to estimate the tapering in breast
%   thickness at edge of breast. Currently only one method implemented, but
%   option is here in case changes are made (it may be main source of error
%   in the method)
%
%  'MaxPairs ([])' 
%   - Set how many pairs of points a user is expected to mark in an image.
%   If not set (i.e. []), defaults to 3 pairs for 18x24 films and 4 pairs
%   for 24x30 (4 pairs are present in both film sizes, but the chest egde
%   pair is nearly always obscured in the smaller films)
%
%  'DoNipple (false)'
%   - Choose whether or not to mark the nipple
%
%  'GenerateTextResults (false)'
%   - Choose whether to generate text files with results of process images
%   (note that full results should always be generated in the form of
%   matlab results files)
%
%  'StepHeights []'
%  - Height (in mm) of each step on the stepwedge, the length of this
%  vector fixes the number of steps
%
%  'NotchIdx []'
%  - index of steps with notches in them (default 4:5:29 39)
%
% 'NumBits (12)'
% - Num of bits in grey scale, default is 12 bits = 0:4905 grey-levels
%
% Outputs: none
%
% Notes: Full documentation for this method is provided in the file...?
%
% See also:...
%
% Created: 09-May-2006
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5125 
% Copyright: (C) University of Manchester

args = u_packargs(varargin, 0, ...
            'MammoNames', [],...
            'MammoDir', [],...
            'SegmentationDir', [],...
            'MarkupDir', [],...
            'ResultsDir', [],...
            'CalibrationName', [],...
            'CalibrationDir', [],...
			'MarkerCalibration', [], ...
            'ResizeFactor', 0.176, ...
            'DebugMode', 0, ...
            'DoAutoMarker', 1,...
            'MarkerTemplate', 'marker_template.mat',...
            'MarkerReductionFactor', 0.25,...
            'MarkerRadius',    50,...
            'ThicknessMethod', 1, ...
            'TaperMethod', 1, ...
            'MaxPairs', [],...
            'DoNipple', 0,...
            'DoAnodeHeel', 0,...
            'GenerateTextResults', 0,...
            'StepHeights', [],...
            'NotchIdx',[],...
            'NumBits', 12,...
            'ImageFormat', 'tif',...
            'User', [] ...
			);
clear varargin;

% this is the amount by which resolution is reduced after markers have been found (44um to 250um)        
resize_factor = args.ResizeFactor;

%Swicth on/off figures to help check method is working
debug_mode = args.DebugMode;

%Select the method for fitting thickness profiles having found the markers
thickness_method = args.ThicknessMethod;

%Select the method for determing the taper at the breast edge
taper_method = args.TaperMethod;

%Set user
if isempty(args.User)
    density_data.user = getenv('username');
else
    density_data.user = args.User;
end

%Assign standard step_heights if none passed by user
if isempty(args.StepHeights)
    step_heights = zeros(39,1);
    step_heights(1:29) = 0.4 + 0.2*((1:29)-1);
    step_heights(30:33) = step_heights(29) + 0.5*(1:4)';
    step_heights(34:39) = step_heights(33) + (1:6)';
else
    step_heights = args.StepHeights;
end
num_steps = length(step_heights);
if isempty(args.NotchIdx)
    notch_idx = [4:5:29 39];
else
    notch_idx = args.NotchIdx;
end

%Load in the template to detect width markers
marker_template = imresize(u_load(args.MarkerTemplate), args.MarkerReductionFactor);

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

%--------------------------------------------------------------------------
% Get user to select mammograms to process and directory locations for
% loading/saving associated data
%--------------------------------------------------------------------------

%select files to process
if isempty(args.MammoNames)
    [mammo_names mammo_dir] = ...
        uigetfile([args.MammoDir '\*.' args.ImageFormat], 'Select the mammograms to process','Multiselect','on');
    if ~isnumeric(mammo_names)
        args.MammoNames = mammo_names;
        args.MammoDir = mammo_dir;
    end
end
%select folder containing breast segmentation
if isempty(args.SegmentationDir)
    segmentation_dir = ...
    uigetdir([], 'Select the directory containing breast segmentations');
    if segmentation_dir
        args.SegmentationDir = segmentation_dir;
    end
end
%select folder to save user data files to
if isempty(args.MarkupDir)
    markup_dir = ...
        uigetdir([], 'Select the directory to save markup data to');
    if markup_dir
        args.MarkupDir = markup_dir;
    end
end
%select folder to density results to
if isempty(args.ResultsDir)
    results_dir = ...
        uigetdir([], 'Select the directory to save results data to');
    if results_dir
        args.ResultsDir = results_dir;
    end
end
%Get calibration data
if isempty(args.CalibrationName)
    [calibration_name, calibration_dir] = ...
        uigetfile([args.CalibrationDir '*.mat'],...
            'Select the calibration data for thickness marker and stepwedge');
    if calibration_name
        args.CalibrationDir = calibration_dir;
        args.CalibrationName = calibration_name;
    end
end
calibration_data = u_load([args.CalibrationDir filesep args.CalibrationName]);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%If we generating text results, open file IDs for the various files
if args.GenerateTextResults
    fid1 = fopen([args.ResultsDir, '/area.txt'],'at');
    fprintf(fid1, '%s \n', ['New session: Stepwedge density with user input:  ' datestr(now)]);
    fprintf(fid1, '%s \n','filename         area_b     area_g     adens(%)');

    fid2 = fopen([args.ResultsDir,'/volume.txt'],'at');
    fprintf(fid2, '%s \n', ['New session: Stepwedge density with user input:  ' datestr(now)]);
    fprintf(fid2, '%s \n','filename         vol_b      vol_g      vdens(%)');

    fid3 = fopen([args.ResultsDir,'/area definition.txt'],'at');
    fprintf(fid3, '%s \n', ['New session: Stepwedge density with user input:  ' datestr(now)]);
    fprintf(fid3, '%s \n','filename         area_g0    area_g5    area_g10  area_g15  area_g20  area_g25  adensg0(%) adensg5(%) adensg10(%) adensg15(%) adensg20(%) adensg25(%)');

    fid4 = fopen([args.ResultsDir, '/breast thickness.txt'],'at');
    fprintf(fid4, '%s \n', ['New session: Stepwedge density with user input:  ' datestr(now)]);
    fprintf(fid4, '%s \n','filename         max_bt     min_bt');
end
%--------------------------------------------------------------------------

%If user selects only 1 file need to make sure name is saved in a cell
%array
if ~iscell(args.MammoNames)
    args.MammoNames = cellstr(args.MammoNames);
end

%Pre-allocate a structure to save failed cases to
failures = [];
n_process = length(args.MammoNames);

%--------------------------------------------------------------------------
%% Main Program Loop
%--------------------------------------------------------------------------
for i_file=1:n_process
    
%--------------------------------------------------------------------------
% Load in the mammogram and work out view, film size etc.
%--------------------------------------------------------------------------

    density_data.time_started = datestr(now);    
    tic;
    
    %Generate a filename to save the user input to
    data_filename = [args.MarkupDir, filesep, args.MammoNames{i_file}(1:end-4), '_data.mat'];
    
    %Check whether data already exists for this file
    if exist(data_filename, 'file');
        
        answer = questdlg(...
                {'Data already exists for this mammogram, continue to next image?';...
                 'Warning, selecting ''No'' will overwrite the old data'},'Mammogram already processed','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        end
    end
    disp(['processing: ' args.MammoDir, args.MammoNames{i_file}]);
    
    %Check whether mammogram is MLO or CC and right or left
    mlo = ~isempty(strfind(args.MammoNames{i_file}, 'LML')) || ~isempty(strfind(args.MammoNames{i_file}, 'RML'));
    left_breast = isempty(strfind(args.MammoNames{i_file}, 'RML')) && isempty(strfind(args.MammoNames{i_file}, 'RCC'));
    
    %Read in the image
    [mammo_large] = imread([args.MammoDir, args.MammoNames{i_file}], args.ImageFormat);
    
    %Reduce mammo size (only need high resolution for locating markers)
    mammo_small = imresize(mammo_large,resize_factor);
    dimensions = size(mammo_small);
    
    %Make sure mammo is upright
    if dimensions(1) < dimensions(2)
        mammo_small = rot90(mammo_small);
        mammo_large = rot90(mammo_large);
        dimensions = [dimensions(2) dimensions(1)];
    end
    
    % workout film size - 18x24 and 24x30 are digitised at different orientations
    % 18x24 films will need to be rotated. We used to do this by checking
    % the name, but better just to check the dimensions as we can't
    % guarantee naming conventions. 18/24 = 0.75, 24/30 = 0.8, so thresh at
    % 0.775
    if dimensions(2)/dimensions(1) < 0.775 %~isempty(strfind(args.MammoNames{i_file}, '1824'))
        filmsizes = 1;
        max_pairs = 3;
    else %~isempty(strfind(args.MammoNames{i_file}, '2430'))
        filmsizes = 2;
        max_pairs = 4;
    end    
    
    %if max_pairs has been set manually, override the automatic choice
    if ~isempty(args.MaxPairs)
        max_pairs = args.MaxPairs;
    end
    
    %Make sure the mammo is oriented correctly
    r1 = round(dimensions(1) / 3);
    r2 = 2*r1;
    c_half = round(dimensions(2) / 2);
    
    mam_mean = mean(mammo_small(:));
    centre_mean = mean(mammo_small(r1:r2,:));
    
    sum1 = sum(centre_mean(1:c_half) > mam_mean);
    sum2 = sum(centre_mean(c_half+1:end) > mam_mean);
    
    %if sum1 > sum2 then the breast is on the left of the image
    if (sum2 > sum1) == left_breast
        mammo_small = rot90(mammo_small, 2);
        mammo_large = rot90(mammo_large, 2);
    end
    
    %Check the image grey-level bit depth (only need to do this on small
    %copy)
    bit_depth = ceil(log2(double(max(mammo_small(:)))));
    if bit_depth < args.NumBits
        disp('Bit depth of mammogram less than target bit depth. Skipping this mammogram');
        failures{end+1,1} = args.MammoNames{i_file}(1:end-4); %#ok
        answer = questdlg(...
            'Continue to next image','Bit depth of mammogram less than target bit depth','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        else
            break;
        end
    end
    mammo_small = floor(mammo_small / 2^(bit_depth - args.NumBits));
    mammo_small = 2^args.NumBits - mammo_small - 1; %12 bits, 0:4095
    
    if args.DoAnodeHeel
        %   Anode heel correction currently not applied but may be in future
        disp('Applying geometric / anode heel correction...');
        
        %See old verison in repository
    end

    density_data.load_time = toc;
    
%--------------------------------------------------------------------------
% Locate the stepwedge
%--------------------------------------------------------------------------
%We do this first now, because if the stepwedge isn't there it is better to
%bail out now before wasting time marking up the rest of the image

    tic;
    % Locate the position of the steps in the step wedge
    [step_centres_xy ui errorcheck] = locate_stepwedge_manual(mammo_small, filmsizes, left_breast, num_steps, notch_idx, debug_mode);
    density_data.stepwedge_markup_time = toc;
    
    %save outputs
    density_data.errorcheck = errorcheck;
    density_data.step_ui = ui;
    density_data.step_centres_xy = step_centres_xy;
    save(data_filename, 'density_data');
    
    if errorcheck
        disp('Error locating stepwedge. Skipping this mammogram');
        failures{end+1,1} = args.MammoNames{i_file}(1:end-4); %#ok
        answer = questdlg(...
            'Continue to next image','Error locating stepwedge','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        else
            break;
        end
    end
    
%--------------------------------------------------------------------------
% Load in breast border and mark pectoral muscle
%--------------------------------------------------------------------------
    tic;
    try
        %load in the existing data
        segmentation_filename = [args.SegmentationDir, filesep, args.MammoNames{i_file}(1:end-4), '_segmentation'];
        load(segmentation_filename);
        
        %resize the breast border
        [breast_border breast_air] = segment_breast_resize(size(mammo_small), segmentation);
        
    catch err
        disp(err.message);
        %disp(['Error loading segmentation file: ', segmentation_filename, '. Check this file exists']);
        %failures{end+1,1} = args.MammoNames{i_file}(1:end-4); %#ok
        %answer = questdlg(...
        %    'Continue to next image','Error loading segmentation','Yes', 'No', 'Yes');
        %if strcmpi(answer, 'yes')
        %    continue;
        %else
        %    break;
        %end;   
        breast_border = [];
        breast_air = [];
    end
    
    %get the user to check the segmented breast is ok - if not skip
    [breast_border, breast_air, errorcheck] = accept_breast_border([], mammo_small, breast_border, breast_air, debug_mode);
    density_data.errorcheck = errorcheck;
    
    if errorcheck
        disp('Error with breast border. Skipping this mammogram');
        failures{end+1,1} = args.MammoNames{i_file}(1:end-4); %#ok
        answer = questdlg(...
            'Continue to next image','Error with breast border','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        else
            break;
        end
    end
    density_data.segmentation_time = toc;
    tic;
    if mlo
        %If MLO mammogram we must also mark the pectoral muscle
        [pectoral_mask pectoral_x pectoral_y] = pectoral_user(mammo_small);
        
        %Save the pectoral masks (less memory to save xy coord of region)
        density_data.pectoral_x = pectoral_x;
        density_data.pectoral_y = pectoral_y;
        save(data_filename, 'density_data');
    end
    density_data.pectoral_time = toc;
    if args.DoNipple
        %mark and save the position of the nipple in case this is used later
        [nipple_position] = nipple_select(mammo_small);
        density_data.nipple_position = nipple_position;
        save(data_filename, 'density_data');
    end

%--------------------------------------------------------------------------
% Locate thickness markers and compute map of breast thickness
%--------------------------------------------------------------------------

    tic;
    % Call markerdetect to locate the thickness markers in the mammo
    if args.DoAutoMarker
        [x_marker, y_marker, selected_markers, marker_ui_xy] = ...
            marker_auto_detect(mammo_large, max_pairs, marker_template, filmsizes==2, left_breast, args.MarkerReductionFactor, debug_mode);   
        r_marker = args.MarkerRadius*ones(size(x_marker));
    else
        [x_marker,y_marker,r_marker,selected_markers,marker_ui_xy,errorcheck] = markerdetect(mammo_large, max_pairs, left_breast, debug_mode);
    end
    density_data.marker_markup_time = toc;
    density_data.marker_ui_xy = marker_ui_xy;
    density_data.errorcheck = errorcheck;
    clear mammo_large;
    
    tic;
    %Check we have marker pair data to continue with, otherwise skip this
    %mammogram
    if errorcheck
        disp('No marker pair data. Skipping this mammogram');
        failures{end+1,1} = args.MammoNames{i_file}(1:end-4); %#ok
        save(data_filename, 'density_data');
        answer = questdlg(...
            'Continue to next image','Error locating markers','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        else
            break;
        end
    end
    
    %Get marker calibration data - either use loaded calibration data or
    %overwrite with specified input arguments
    if isempty(args.MarkerCalibration)        
        M_mag = calibration_data(filmsizes).M_mag; % units = pixels per mm
        C_mag = calibration_data(filmsizes).C_mag; % units = pixels
    else
        M_mag = args.MarkerCalibration(filmsizes).M_mag;
        C_mag = args.MarkerCalibration(filmsizes).C_mag;
    end
    
    %For left breasts, the order of the markers will be reversed so flip
    %calibration coefficients
    if left_breast
        M_mag = M_mag(end:-1:1);
        C_mag = C_mag(end:-1:1);
    end
    
    %Use calibration data to compute thickness across mammo (not specific
    %to breast shape yet)
    [thickness_profile] = thickness_from_markers(x_marker, y_marker, M_mag(selected_markers), C_mag(selected_markers));
    
    density_data.thickness_profile = thickness_profile;
    density_data.x_marker = x_marker;
    density_data.y_marker = y_marker;
    density_data.r_marker = r_marker;
    save(data_filename, 'density_data');

    % make a mask with ones inside marker and zeros outside - we'll use
    % this later to mask the breast so the bright white markers don't
    % influence the density calculations
    [marker_mask] = make_marker_mask(mammo_small, x_marker, y_marker, r_marker, resize_factor);

    %MB at Jenny's request this is now calculated using a straight line
    %fitted to the thickness, rather than interpolating between the
    %thicknesses - however the choice is there to switch between the
    %two
    switch thickness_method

        case 1
             %default method fit straight line
             breast_thickness = thickness_from_polyfit(thickness_profile,dimensions,resize_factor, debug_mode);
        case 2
             %linearly interpolate between each point
             breast_thickness = thickness(thickness_profile,dimensions,resize_factor);

        otherwise
            display('Not a valid thickness profile method, fitting straight line');
            breast_thickness = thickness_from_polyfit(thickness_profile,dimensions,resize_factor, debug_mode);
    end
    
    %Now call edge_profile to fit a tapered thickness at the breast
    %edge (designed to allow us to choose different methods in future)
    switch taper_method

        case 1
            %default method fit a semicircular profile
            [breast_thickness taper_region_idx] =...
                compute_edge_profile(breast_thickness, breast_border, breast_air, left_breast, resize_factor, mlo, debug_mode);          

        otherwise
            display('Not a valid taper method, fitting semicircular profiles');
            [breast_thickness taper_region_idx] =...
                compute_edge_profile(breast_thickness, breast_border, breast_air, left_breast, resize_factor, mlo, debug_mode);
    end
    
    %Throw away pixels that belong to image markers from the breast
    %thickness
    breast_thickness(marker_mask) = 0;
    
    %If MLO, throw away pectoral muscle
    if mlo
        breast_thickness(pectoral_mask) = 0;
    end
    density_data.marker_processing_time = toc;
%--------------------------------------------------------------------------
% Compute grey-levels associated with each stepwedge step
%--------------------------------------------------------------------------   
    tic;
    % Compute the mean grey-level associated with each step
    [step_g, errorcheck] = calculate_step_values(mammo_small, step_centres_xy, debug_mode);
    
    %save outputs
    density_data.errorcheck = errorcheck;
    density_data.step_ui = ui;
    density_data.step_g = step_g;
    density_data.step_centres_xy = step_centres_xy;
    save(data_filename, 'density_data');
    
    %Check there were no problems in calculating wedge values
    if errorcheck
        disp('Error in stepwedge data. Skipping this mammogram');
        failures{end+1,1} = args.MammoNames{i_file}(1:end-4); %#ok
        answer = questdlg(...
            'Continue to next image','Error in stepwedge data','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        else
            break;
        end
    end

    % create lookup table to convert pixel grey-level to step height
    sw_lookup = compute_sw_lookup(step_g, step_heights, args.NumBits);
    density_data.stepwedge_processing_time = toc;
    save(data_filename, 'density_data');
    
%--------------------------------------------------------------------------
% Use lookup table and breast thickness map to compute gland thickness
%--------------------------------------------------------------------------

    tic;
    %Can now produce gland_thickness
    gland_thickness = ...
        sw_vs_gland_thick(sw_lookup(mammo_small+1), breast_thickness,...
        calibration_data(filmsizes).stepwedge, debug_mode);

    % gland_thickness array contains NaN - so convert NaN to 0 or volume calculation will give NaN
    gland_thickness(isnan(gland_thickness)) = 0;
    density_data.gland_processing_time = toc;
    
    %Display the final gland map alongside the original mammogram
    gland_fig = figure(...
        'Name', 'Gland thickness map',...
        'Units', 'normalized',...
        'OuterPosition', [0 0 1 1]); 
    subplot(1,2,1); 
    imagesc(-double(mammo_small)); axis image; colormap(gray(256)); 
    title('Original mammogram');
    subplot(1,2,2);        
    imagesc(gland_thickness); axis image; colormap(gray(256));
    if left_breast
        colorbar('east', 'xcolor', 'r', 'ycolor', 'r');
    else
        colorbar('west', 'xcolor', 'r', 'ycolor', 'r');
    end
    title('Gland thickness map');       

    %Generate results, this also print the results to the correct files
    if args.GenerateTextResults
        [density_results] = ...
            compute_density_results(gland_thickness, breast_thickness, taper_region_idx,...
                fid1, fid2, fid3, fid4, args.MammoNames{i_file});
    else
        [density_results] = ...
            compute_density_results(gland_thickness, breast_thickness, taper_region_idx);
    end

    %Also store the thickness profile and taper methods used
    density_results.thickness_method = thickness_method;
    density_results.taper_method = taper_method;

    %Generate a filename to save the density results to
    results_filename = [args.ResultsDir, filesep, args.MammoNames{i_file}(1:end-4), '_results'];
    save(results_filename, 'density_results');

    %Total time for this mammo
    density_data.time_finished = datestr(now);
    save(data_filename, 'density_data');
    
    %Clear density_data and density_resuts ready for next mammogram
    clear density_data density_results;
    
    answer = questdlg(...
            'Continue to next image','Mammogram Complete','Yes', 'No', 'Yes');
    if ~debug_mode
        close(gland_fig);
    end
    
    if strcmpi(answer, 'no')
        break;
    end
     
end
%--------------------------------------------------------------------------
% End of main loop - tidy up below on exiting program
%--------------------------------------------------------------------------
%%
helpdlg('Finished processing all selected mammograms','Finished!')
disp('Finished processing all selected mammograms');

%turn warnings back on
warning on Images:initSize:adjustingMag;
warning on Images:imshow:magnificationMustBeFitForDockedFigure
warning on MATLAB:Figure:SetPosition
set(0,'DefaultFigureWindowStyle',orig_window_style);

if args.GenerateTextResults
    %Close the file IDs for the results text files
    fclose(fid1);
    fclose(fid2);
    fclose(fid3);
    fclose(fid4);
end

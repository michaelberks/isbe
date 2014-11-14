 function stepwedge_calibrate(varargin)

%Now using u_parkargs interface
args = u_packargs(varargin, 0, ...
            'ImageFormat', 'tif',...
            'NumBits', 12,...
            'ResizeFactor', 0.176, ...
            'StepHeights', [],...
            'NotchIdx', [],...
            'DebugMode', 0, ...
            'User', [] ...
			);

% this is the amount by which resolution is reduced after markers have been found (44um to 250um)        
resize_factor = args.ResizeFactor;

%Swicth on/off figures to help check method is working
debug_mode = args.DebugMode;

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

try
    fid = fopen('C:\program files\density_software\default_locations.txt');
    def_txt = textscan(fid, '%s %s', 'delimiter', '=');
    fclose(fid);
    default_dir = def_txt{2}{strncmpi(def_txt{1}, 'stepwedge_calibrate', length('stepwedge_calibrate'))};
catch
    display('Problem loading default locations');
    default_dir = [];
end

%select files to process
%disp('Select the files you wish to process');
[mammo_names, mammo_dir] = uigetfile([default_dir '\*' args.ImageFormat],'Select the mammograms to process','Multiselect','on');

%select folder to save user data files to
%disp('Select the directory you wish to save the user input to');
data_dir = uigetdir(default_dir, 'Select the directory to save the user input to');

% loop over the selected mammograms

%If user selects only 1 file need to make sure name is saved in a cell
%array
if ~iscell(mammo_names)
    temp = mammo_names;
    mammo_names = cell(1);
    mammo_names{1} = temp; clear temp;
end

van_name = mammo_names{1}(1:4);

calibration_data_filename = [data_dir, filesep, van_name, '_calibration_data.mat'];
%--------------------------------------------------------------------------
% Set up data structures for storing calibration data
%--------------------------------------------------------------------------
if exist(calibration_data_filename, 'file');
    calibration_data = u_load(calibration_data_filename); 
    if ~isfield(calibration_data, 'stepwedge');
        calibration_data(2).stepwedge = [];
    end
else
    calibration_data(2).stepwedge = [];
end

fid_1824 = fopen([data_dir, filesep, van_name, '_calibration_data_1824.txt'],'at');
fid_2430 = fopen([data_dir, filesep, van_name, '_calibration_data_2430.txt'],'at');
%--------------------------------------------------------------------------

%Pre-allocate a structure to save failed cases to
failures = [];
n_process = length(mammo_names);
%
for i_file=1:n_process
    
    %Generate a filename to save the user input to
    data_filename = [data_dir, filesep, mammo_names{i_file}(1:end-4), '_data.mat'];
    
    %Check whether data already exists for this file
    if exist(data_filename, 'file');
        
        answer = questdlg(...
                {'Data already exists for this calibration image, continue to next image?';...
                 'Warning, selecting ''No'' will overwrite the old data'},'Mammogram already processed','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        end
    end
    disp(['processing: ' mammo_dir, mammo_names{i_file}]);
    
    h = waitbar(0,'Loading and processing image. Please wait...');
    
    %----------------------------------------------------------------------
    % workout film size - 18x24 and 24x30 are digitised at different orientations
    % 18x24 films will need to be rotated
    if ~isempty(strfind(mammo_names{i_file}, '1824'))
        filmsizes = 1;
    elseif ~isempty(strfind(mammo_names{i_file}, '2430'))
        filmsizes = 2;
    else
        %define error action
    end
    
    % read in the image
    [mammo] = imread([mammo_dir, mammo_names{i_file}],'tif');
    
    if filmsizes == 1
        % this is only needed for the 1824 film sizes
        mammo = imrotate(mammo,90);
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------

    % reduce image sizes here (don't need such high resolution now magnification markers have been located)
    mammo = imresize(mammo,resize_factor); 

    % Now reduce the image to 4096 grey levels and invert to make it compatible
    % with sw_lookup_table.m
    bit_depth = ceil(log2(double(max(mammo(:)))));
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
    mammo = floor(mammo / 2^(bit_depth - args.NumBits));
    mammo = 2^args.NumBits - mammo - 1; %12 bits, 0:4095
    
    close(h);
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Locate stepwedge and compute average intensity of each step
    [step_centres_xy] = locate_stepwedge_manual(mammo, filmsizes, 0, num_steps, notch_idx, debug_mode);
    [step_g, errorcheck] = calculate_step_values(mammo, step_centres_xy, debug_mode);
    
    %save outputs
    film_calibration_data.errorcheck = errorcheck;
    film_calibration_data.step_g = step_g;
    film_calibration_data.step_centres_xy = step_centres_xy;
    save(data_filename, 'film_calibration_data');
    
    %Check there were no problems in calculating wedge values
    if errorcheck
        disp('Error in stepwedge data. Skipping this mammogram');
        failures{end+1,1} = mammo_names{i_file}(1:end-4); %#ok
        answer = questdlg(...
            'Continue to next image','Error in stepwedge data','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            continue;
        else
            break;
        end
    end
    %----------------------------------------------------------------------
    % Compute step height for each grey-level intensity and match to
    % intensity of claibration image
    
    % create lookup table to convert pixel value to x_sw
    sw_lookup = compute_sw_lookup(step_g, step_heights, args.NumBits);
    
    if filmsizes == 1
        r1 = 448;
        r2 = 488;
        c1 = 524;
        c2 = 563;
    else
        r1 = 569;
        r2 = 609;
        c1 = 772;
        c2 = 811;
    end
    
    phantom_area = mammo(r1:r2, c1:c2);
    phantom_intensity = round(mean(phantom_area(:)));   
    step_height = sw_lookup(round(phantom_intensity));
    %step_height = phantom_intensity / 100;
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    %Display outputs to the user
    calib_fig = figure; 
    imagesc(mammo); axis image; colormap(gray(256)); hold on;
    plot([c1 c2 c2 c1 c1], [r1 r1 r2 r2 r1], 'r', 'LineWidth', 2);
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Get fat, gland and total thickness from the filename
    fat_idx = find(mammo_names{i_file} == 'f', 1);
    f_error = true;
    if ~isempty(fat_idx)
        fat_thickness = str2double(mammo_names{i_file}(fat_idx-2:fat_idx-1));
        if ~isnan(fat_thickness)
            f_error = false;
        end
    end

    gland_idx = find(mammo_names{i_file} == 'g', 1);
    g_error = true;
    if ~isempty(gland_idx)
        gland_thickness = str2double(mammo_names{i_file}(gland_idx-2:gland_idx-1));
        if ~isnan(fat_thickness)
            g_error = false;
        end
    end
    
    if ~f_error && ~g_error
        answer = questdlg(...
                {['Fat and gland thicknesses for file ', mammo_names{i_file}];...
                 '';...
                 ['Fat = ', num2str(fat_thickness)];...
                 ['Gland = ', num2str(gland_thickness)];...
                 ['Breast thickness = ', num2str(fat_thickness+gland_thickness)];...
                 '';...
                 'Is this correct?'},'Thickness data','Yes', 'No', 'Yes');
        if strcmpi(answer, 'no')
            f_error = true;
        end
    end
    if f_error
        %manually enter thickness data
        prompt = {'What is the fat thickness?','What is the gland thickness?:'};
        dlg_title = 'Select fat and gland thicknesses';
        num_lines = 1;
        screen_size = get(0,'ScreenSize');
        options.Position = [screen_size(3)/4, 3*screen_size(4)/5];
        answer = inputdlg(prompt,dlg_title,num_lines,{'',''}, options);
        fat_thickness = str2num(answer{1}); %#ok
        gland_thickness = str2num(answer{2}); %#ok
            
    end
    breast_thickness = fat_thickness + gland_thickness;
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Save calibration data - do this after each image so we don't lose
    % existing data in a crash

    %Check whether this pair of breast/gland thicknesses already exists
    if ~isempty(calibration_data(filmsizes).stepwedge)
        pair_exists = ...
            calibration_data(filmsizes).stepwedge(:,1) == breast_thickness &...
            calibration_data(filmsizes).stepwedge(:,2) == gland_thickness;
    else
        pair_exists = false;
    end
    
    if any(pair_exists)
        %If so, overwrite the step height data
        calibration_data(filmsizes).stepwedge(pair_exists,3) = step_height;
        
    else
        %Otherwise append to end
        calibration_data(filmsizes).stepwedge = [...
            calibration_data(filmsizes).stepwedge;...
            breast_thickness, gland_thickness, step_height];
    end
    save(calibration_data_filename, 'calibration_data');
    
    if filmsizes == 1
        fprintf(fid_1824,'%s %7.4f %7.4f \n',mammo_names{i_file}(1:end-4),phantom_intensity,step_height);
    else
        fprintf(fid_2430,'%s %7.4f %7.4f \n',mammo_names{i_file}(1:end-4),phantom_intensity,step_height);
    end

    answer = questdlg(...
            'Continue to next image','Mammogram Complete','Yes', 'No', 'Yes');
    
    if strcmpi(answer, 'no')
        break;
    end
    close(calib_fig);
     
end

%--------------------------------------------------------------------------
% Plot the relationship between step heights, breast thickness and gland
% thickness
figure;
size_labels = {'18x24','24x30'};
for filmsizes = 1:2
    [sh bt] = meshgrid(linspace(0,14,100), linspace(0,max(calibration_data(filmsizes).stepwedge(:,1)),100));
    try
        [gt] = sw_vs_gland_thick(sh, bt, calibration_data(filmsizes).stepwedge, 0);
    catch %#ok
        continue;
    end
    subplot(1,2,filmsizes);
    mesh(sh,bt,gt);
    xlabel('Step heights (mm)');
    ylabel('Breast thickness (mm)');
    zlabel('Gland thickness (mm)');
    title([van_name, ' stepwedge calibration, ', size_labels{filmsizes}, ' films']);
end   
    
h = helpdlg('Finished processing all selected mammograms','Finished!');
disp('Finished processing all selected mammograms');
uiwait(h);

%turn warnings back on
warning on Images:initSize:adjustingMag;
warning on Images:imshow:magnificationMustBeFitForDockedFigure
warning on MATLAB:Figure:SetPosition
set(0,'DefaultFigureWindowStyle',orig_window_style);

%Close the file IDs for the results text files
fclose(fid_1824);
fclose(fid_2430);

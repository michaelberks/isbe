function failures = stepwedge_from_saved_old(varargin)

%Now using u_parkargs interface
args = u_packargs(varargin, 0, ...
			'M_mag', [], ...
			'C_mag', [], ...
            'ResizeFactor', 0.176, ...
            'DebugMode', 0, ...
            'ThicknessMethod', 1, ...
            'TaperMethod', 1, ...
            'DisplayGlandFig', 0, ...
            'SelectionMode', 1, ...
            'Save', 1 ...
			);

orig_window_style = get(0,'DefaultFigureWindowStyle');
if ~strcmp(orig_window_style, 'docked')
    display('Warning: changing window style to normal for function');
    set(0,'DefaultFigureWindowStyle','docked');
end

% this is the amount by which resolution is reduced after markers have been found (44um to 250um)        
resize_factor = args.ResizeFactor;

%Swicth on/off figures to help check method is working
debug_mode = args.DebugMode;

%Select the method for fitting thickness profiles having found the markers
thickness_method = args.ThicknessMethod;

%Select the method for determing the taper at the breast edge
taper_method = args.TaperMethod;

format bank; % don't scale by 1000
iptsetpref('ImshowBorder','tight'); % show less border round image

% all films are digitised as 'right' images -maintain for
% back-compatibility
leftright = 0;

switch args.SelectionMode
    case 1
        %select files to process
        %disp('Select the images you wish to process');
        [filenames, image_path] = uigetfile('G:\*.tif','select mammograms to process','Multiselect','on');

        %select folder to save thicknesses/sw locations to
        %disp('Select the directory containing the saved user data');
        data_path = uigetdir([], 'Select the directory containing the saved user data');
        
        %select folder to save thicknesses/sw locations to
        %disp('Select the directory you wish to the save the results to');
        if args.Save
            results_path = uigetdir(data_path, 'Select the directory to the save the results to');
        end
    case 2
        %select files to process
        %disp('Select the user data files you wish to process');
        [filenames, data_path] = uigetfile('J:\*.mat','select user data files to process','Multiselect','on');

        %select folder to save thicknesses/sw locations to
        %disp('Select the directory containing the saved user data');
        image_path = uigetdir('G:\', 'Select the directory containing the mammograms');

        %select folder to save thicknesses/sw locations to
        %disp('Select the directory you wish to the save the results to');
        if args.Save
            results_path = uigetdir(data_path, 'Select the directory to the save the results to');
        end
    case 3
        %select files to process
        %disp('Select the images you wish to process');
        [filenames, results_path] = uigetfile('J:\*.mat','select the results files to process','Multiselect','on');
        
        %select folder to save thicknesses/sw locations to
        %disp('Select the directory you wish to the save the results to');
        image_path = uigetdir('G:\', 'Select the directory containing the mammograms');
        
        %select folder to save thicknesses/sw locations to
        %disp('Select the directory containing the saved user data');
        data_path = uigetdir(results_path, 'Select the directory containing the saved user data');
end

if args.Save
    %Open file ids to write results
    fid1 = fopen([results_path, '/area.txt'],'at');
    fprintf(fid1, '%s \n', ['New session: Stepwedge density from saved data:  ' datestr(now)]);
    fprintf(fid1, '%s \n','filename         area_b     area_g     adens(%)');

    fid2 = fopen([results_path,'/volume.txt'],'at');
    fprintf(fid2, '%s \n', ['New session: Stepwedge density from saved data:  ' datestr(now)]);
    fprintf(fid2, '%s \n','filename         vol_b      vol_g      vdens(%)');

    fid3 = fopen([results_path,'/area definition.txt'],'at');
    fprintf(fid3, '%s \n', ['New session: Stepwedge density from saved data:  ' datestr(now)]);
    fprintf(fid3, '%s \n','filename         area_g0    area_g5    area_g10  area_g15  area_g20  area_g25  adensg0(%) adensg5(%) adensg10(%) adensg15(%) adensg20(%) adensg25(%)');

    fid4 = fopen([results_path, '/breast thickness.txt'],'at');
    fprintf(fid4, '%s \n', ['New session: Stepwedge density from saved data:  ' datestr(now)]);
    fprintf(fid4, '%s \n','filename         max_bt     min_bt');
end 
% loop over the files

%If user selects only 1 file need to make sure name is saved in a cell
%array
if ~iscell(filenames)
    temp = filenames;
    filenames = cell(1);
    filenames{1} = temp; clear temp;
end

failures = [];

n_process = length(filenames);
%%
for i_file=1:n_process
    
    switch args.SelectionMode
        case 1
            image_filename = [image_path, filesep, filenames{i_file}];
            data_filename = [data_path, filesep, filenames{i_file}(1:end-4), '_data'];
            if args.Save
                results_filename = [results_path, filesep, filenames{i_file}(1:end-4), '_results'];
            end
        case 2
            data_filename = [data_path, filesep, filenames{i_file}];
            image_filename = [image_path, filesep, filenames{i_file}(1:end-9), '.tif'];
            if args.Save
                results_filename = [results_path, filesep, filenames{i_file}(1:end-9), '_results'];
            end
        case 3
            if args.Save
                results_filename = [results_path, filesep, filenames{i_file}];
            end
            image_filename = [image_path, filesep, filenames{i_file}(1:end-12), '.tif'];
            data_filename = [data_path, filesep, filenames{i_file}(1:end-12), '_data'];
    end
    
    
    disp(['processing: ' image_filename]);
    
    %Don't crash if data file doesn't exist
    try
        %load in the existing data
        load(data_filename);
    catch
        disp(['Error loading data file: ', data_filename, '. Check this file exists']);
        failures{end+1,1} = data_filename; %#ok
        continue;       
    end
    
    %Check marker pair data exists (and hence errorcheck is neg)
    if do_errorcheck(density_data)
        display(['Incorrect user_data file: ', data_filename, '. Skipping this mammogram']);
        failures{end+1,1} = data_filename; %#ok
    else

        % define film size
        % 18x24 and 24x30 are digitised at different orientations
        % 18x24 films will need to be flipped
        if ~isempty(strfind(filenames{i_file}, '1824'))
            filmsizes = 1;
        elseif ~isempty(strfind(filenames{i_file}, '2430'))
            filmsizes = 2;
        else
            %define error action
        end
        
        try
            % read in the image
            [IMAGE] = imread(image_filename);
        catch
            disp(['Error loading mammogram: ', image_filename, '. Check this file exists']);
            failures{end+1,1} = data_filename; %#ok
            continue;
        end
        if filmsizes == 1
            % this is only needed for the 1824 film sizes
            IMAGE = imrotate(IMAGE,90);
        end

        % reduce image sizes here (don't need such high resolution now magnification markers have been located)
        IMAGE = imresize(IMAGE,resize_factor);
        dimensions = size(IMAGE);

        %MB at Jenny's request this is now calculated using a straight line
        %fitted to the thickness, rather than interpolating between the
        %thicknesses - however the choice is there to switch between the
        %two
        try
            switch thickness_method

                case 1
                     %default method fit straight line
                     x_b = thickness_from_polyfit(density_data.x_b_info,dimensions,resize_factor, debug_mode);
                case 2
                     %linearly interpolate between each point
                     x_b = thickness(density_data.x_b_info,dimensions,resize_factor);

                otherwise
                    display('Not a valid thickness profile method, fitting straight line');
                    x_b = thickness_from_polyfit(density_data.x_b_info,dimensions,resize_factor, debug_mode);
            end

            % Now reduce the image to 4096 grey levels and invert to make it compatible
            % with sw_lookup_table.m
            IMAGE = IMAGE./16;
            IMAGE = 4095-IMAGE;

            % create lookup table to convert pixel value to x_sw
            sw_lookup = sw_lookup_table_12bit(density_data.wedgevals);

            % get mask of breast area
            [background_mask, edgex, edgey] = breast_edge(IMAGE,...
                density_data.coarse_edgex, density_data.coarse_edgey, leftright, debug_mode);
            background_mask = double(background_mask > 0);
        catch
            disp(['Error processing: ', data_filename, '. Check user data file']);
            failures{end+1,1} = data_filename; %#ok
            continue;
        end
        
        %now call edge_profile to fit a tapered thickness at the breast edge
        %Designed to allow us to choose different methods in future
        switch taper_method

            case 1
                %default method fit a semicircular profile
                [breast_thickness taper_region_idx] =...
                    mb_edge_profile_old(x_b,background_mask,edgex,edgey,leftright,resize_factor, debug_mode);

            otherwise
                display('Not a valid taper method, fitting semicircular profiles');
                [breast_thickness taper_region_idx] =...
                    mb_edge_profile_old(x_b,background_mask,edgex,edgey,leftright,resize_factor, debug_mode);
        end

        % then produce x_gland image
        gland_thickness = sw_vs_gland_thick(sw_lookup(IMAGE+1), breast_thickness);

        % gland_thickness array contains NaN - so convert NaN to 0 or 
        % volume calculation will give NaN
        gland_thickness(isnan(gland_thickness)) = 0;

        if args.DisplayGlandFig
            figure('Name', data_filename); 
            subplot(1,2,1); 
            imagesc(-double(IMAGE)); axis image; colormap(gray(256)); 
            title('Original mammogram');

            subplot(1,2,2);        
            imagesc(gland_thickness); axis image; colormap(gray(256)); 
            colorbar('west', 'xcolor', 'r', 'ycolor', 'r');
            title('Gland thickness map');
        end
        
        if args.Save
            %Generate results, this also print the results to the correct files
            [density_results] = ...
                compute_density_results(gland_thickness, breast_thickness, taper_region_idx,...
                    fid1, fid2, fid3, fid4, filenames{i_file}); %#ok

            %Also store the thickness profile and taper methods used
            density_results.thickness_method = thickness_method;
            density_results.taper_method = taper_method;
        
        
            %save the density results
            save(results_filename, 'density_results');
            
            %Clear density_data and density_resuts ready for next mammogram
            clear density_results;
        end
            
    end
    
end
disp('Finished processing all selected mammograms');

%Close the file IDs for the results text files
if args.Save
    fclose(fid1);
    fclose(fid2);
    fclose(fid3);
    fclose(fid4);
end

set(0,'DefaultFigureWindowStyle',orig_window_style);

end
%--------------------------------------------------------------------------
function data_error = do_errorcheck(density_data)

    data_error = ...
        ~isfield(density_data, 'errorcheck') || ...
        density_data.errorcheck || ...
        ~isfield(density_data, 'user') || ...
        ~isfield(density_data, 'x_b_info') || ...
        ~isfield(density_data, 'nipple_position') || ...
        ~isfield(density_data, 'wedgevals') || ...
        ~isfield(density_data, 'swx') || ...
        ~isfield(density_data, 'swy') || ...
        ~isfield(density_data, 'coarse_edgex') || ...
        ~isfield(density_data, 'coarse_edgey');
end

    


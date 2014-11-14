function [] = convert_image_format()
%DENSITY_WRITE_VAS_FORMS implements a GUI to create VAS density forms
%   [] = density_write_VAS_forms()
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
% Created: 17-Jul-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%create args to be passed to form reader function
args.input_type = 1;
args.input_bit = 1;
args.output_type = 1;
args.output_bit = 1;
args.num_images = 0;
args.bit_lookup = [8; 12; 16];
args.type_lookup = {'bmp', 'tif', 'dcm'};

%Create the main figure;
main_fig = figure(...
    'Name', 'Convert images',...
    'WindowStyle', 'normal',...
    'MenuBar', 'none',...
    'NumberTitle', 'off',...
    'Position', [100 100 840 160]);

%Create the UI controls for the main fig
input_images_box = uicontrol(...
    'Style', 'text',...
    'Position', [100 110 400 40],...
    'BackgroundColor', [1 1 1],...
    'Parent', main_fig,...
    'String', 'No images selected');

input_images_select = uicontrol(... 
    'Style', 'pushbutton',...
    'Position', [510 110 100 40],...
    'String', 'Select images',...
    'Parent', main_fig,...
    'Callback', @input_images_select_Callback);%#ok

input_images_bit = uicontrol(... 
    'Style', 'listbox',...
    'Position', [620 110 100 40],...
    'String', {'8', '12', '16'},...
    'Parent', main_fig,...
    'Callback', @input_images_bit_Callback);%#ok

input_images_type = uicontrol(... 
    'Style', 'listbox',...
    'Position', [730 110 100 40],...
    'String', args.type_lookup,...
    'Parent', main_fig,...
    'Callback', @input_images_type_Callback);%#ok

input_images_text = uicontrol(... 
    'Style','text',...
    'BackgroundColor', get(main_fig, 'Color'),...
    'FontName', 'Arial',...
    'String', 'Input images:',...
    'HorizontalAlignment', 'right',...
    'Position', [10 110 80 25]);%#ok

output_dir_box = uicontrol(...
    'Style', 'edit',...
    'Position', [100 60 400 40],...
    'BackgroundColor', [1 1 1],...
    'Parent', main_fig,...
    'String', []);

output_dir_select = uicontrol(... 
    'Style', 'pushbutton',...
    'Position', [510 60 100 40],...
    'String', 'Select folder',...
    'Parent', main_fig,...
    'Callback', @output_dir_select_Callback);%#ok

output_dir_text = uicontrol(... 
    'Style','text',...
    'BackgroundColor', get(main_fig, 'Color'),...
    'FontName', 'Arial',...
    'String', 'Output images:',...
    'HorizontalAlignment', 'right',...
    'Position', [10 60 80 25]);%#ok

output_dir_box = uicontrol(...
    'Style', 'edit',...
    'Position', [100 60 400 40],...
    'BackgroundColor', [1 1 1],...
    'Parent', main_fig,...
    'String', []);

output_images_bit = uicontrol(... 
    'Style', 'listbox',...
    'Position', [620 60 100 40],...
    'String', {'8', '12', '16'},...
    'Parent', main_fig,...
    'Callback', @output_images_bit_Callback);%#ok

output_images_type = uicontrol(... 
    'Style', 'listbox',...
    'Position', [730 60 100 40],...
    'String', args.type_lookup,...
    'Parent', main_fig,...
    'Callback', @output_images_type_Callback);%#ok

convert_images = uicontrol(... 
    'Style', 'pushbutton',...
    'Position', [100 10 200 40],...
    'String', 'Convert images',...
    'Parent', main_fig,...
    'Callback', @convert_images_Callback);%#ok

% advanced_options = uicontrol(...
%     'Style', 'pushbutton',...
%     'Position', [310 10 100 40],...
%     'Parent', main_fig,...
%     'String', 'Advanced options');


%specify default locations for browser to open
try
    fid = fopen('C:\program files\density_software\default_locations.txt');
    def_txt = textscan(fid, '%s %s', 'delimiter', '=');
    fclose(fid);

    default_output_dir = ...
        def_txt{2}{strncmpi(def_txt{1}, 'output_images_dir', length('output_images_dir'))};
    default_input_dir = ...
        def_txt{2}{strncmpi(def_txt{1}, 'input_images_dir', length('input_images_dir'))};

catch
    display('Problem opening default file locations');
    default_output_dir = '';
    default_input_dir = '';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks for UI controls
% --------------------------------------------------------------------
    function output_dir_select_Callback(hObject, eventdata) %#ok
    % Callback to...
        args.output_dir = ...
            uigetdir(default_output_dir, 'Select the directory to the save the forms to');
        if args.output_dir
            set(output_dir_box, 'string', args.output_dir);
        end
    end
%--------------------------------------------------------------------------
    function input_images_select_Callback(hObject, eventdata) %#ok
    % Callback to...
        [args.input_names, args.input_dir] = uigetfile(...
            ['*.' args.type_lookup{args.input_type}],...
            'Select input images to convert',default_input_dir,...
            'MultiSelect', 'on');
        
        if ischar(args.input_names)
            args.input_names = cellstr(args.input_names);
        end
        if iscell(args.input_names)
            args.num_images = length(args.input_names);
            set(input_images_box, 'string', [num2str(args.num_images) ' images selected']);
        end
    end

%--------------------------------------------------------------------------
    function input_images_bit_Callback(hObject, eventdata) %#ok
    % Callback to...
        args.input_bit = get(input_images_bit, 'value');
        %display(['Input bit depth is ' num2str(args.bit_lookup(args.input_bit))]);
        
    end

%--------------------------------------------------------------------------
    function input_images_type_Callback(hObject, eventdata) %#ok
    % Callback to...
        args.input_type = get(input_images_type, 'value');
        %display(['Input image type is ' args.type_lookup{args.input_type}]);
    end

%--------------------------------------------------------------------------
    function output_images_bit_Callback(hObject, eventdata) %#ok
    % Callback to...
        temp_bit = get(output_images_bit, 'value');
        if (temp_bit ~= 1) && (args.output_type ~= 3)
            warndlg({...
                'Images can only be saved in 12 or 16 bit format as DICOM images.';
                'If you wish to select this option, please first select ''dcm'' from the output file types.'},...
                'Bit depth and file type incompatible');
            set(output_images_bit, 'value', 1);
        else
            args.output_bit = temp_bit;
        end
        %display(['Output bit depth is ' num2str(args.bit_lookup(args.output_bit))]);
        
    end

%--------------------------------------------------------------------------
    function output_images_type_Callback(hObject, eventdata) %#ok
    % Callback to...
        temp_type = get(output_images_type, 'value');
        if (temp_type ~= 3) && (args.output_bit ~= 1)
            warndlg({...
                'Images can only be saved in 12 or 16 bit format as DICOM images.';
                'If you wish to select this option, please first select 8 from the output bit depths.'},...
                'Bit depth and file type incompatible');
            set(output_images_type, 'value', 3);
        else
            args.output_type = temp_type;
        end 
        %display(['Output image type is ' args.type_lookup{args.output_type}]);
    end


%--------------------------------------------------------------------------
    function convert_images_Callback(hObject, eventdata) %#ok
    % Callback to...
        output_folder = [get(output_dir_box, 'string'), filesep];
        
        if isempty(output_folder)
            warndlg('Please select an excel file containing participant data',...
                'Missing data');
        else
            set(get(main_fig, 'Children'), 'Enable', 'off');
            h = waitbar(0,'Creating forms. Please wait...');
            
            warnings = 0;
            
            %Do actual conversions
            for ii = 1:args.num_images
                
                try

                    %load in image
                    [dummy, image_name ext] = fileparts(args.input_names{ii});
                    if strcmp(ext, '.dcm')
                        image_in = dicomread([args.input_dir args.input_names{ii}]);
                    else
                        image_in = imread([args.input_dir args.input_names{ii}]);
                    end

                    %Check image isn't in RGB form
                    if size(image_in, 3) == 3;
                        image_in = rgb2gray(image_in);
                    end
                    
                    if max(image_in(:)) > 2^args.bit_lookup(args.input_bit)
                        display(['Warning: image ', args.input_names{ii}, 'appears to have a larger ',...
                            'input bit depth than specified.']);
                        warnings = 1;
                    end       

                    %Convert to double so we can easily scale image
                    image_in = double(image_in);

                    %check image bit depth and convert if necessary
                    conversion_factor = 2^(args.bit_lookup(args.output_bit) - args.bit_lookup(args.input_bit));
                    image_in = image_in * conversion_factor;

                    if args.bit_lookup(args.output_bit) == 8
                        image_in = uint8(image_in);
                    else
                        image_in = image_in ./ ((2^args.bit_lookup(args.output_bit))-1);
                    end

                    %form new image name depending on output file type
                    new_image_name = [args.output_dir filesep image_name '.' args.type_lookup{args.output_type}];

                    %save new images
                    if strcmp(args.type_lookup{args.output_type}, 'dcm')
                         dicomwrite(image_in, new_image_name);
                    else
                        imwrite(image_in, new_image_name);
                    end
                    
                    display(['Image ', num2str(ii), ' saved as: ', new_image_name]);
                catch
                    display(lasterr);
                    display(['Skipping image: ', args.input_names{ii}, '. Please check file type and bit depth.']);
                end
            end
            
            if warnings
                warndlg('One or more images had inconsistent input bit depths. Please check ouput images.',...
                    'Problems converting images');
            end
            
            close(h);
            answer = questdlg(...
                'Do you have more images to convert or would you like to exit?',...
                'Images successfully converted!','Continue', 'Exit', 'Exit');
            if strcmpi(answer, 'exit');
                close(main_fig);
            else
                set(get(main_fig, 'Children'), 'Enable', 'on');
            end
        end
        
            
        
        
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---------------------- END OF FUNCTION -----------------------------------
%--------------------------------------------------------------------------
end

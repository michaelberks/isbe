function stepwedge_pectoral_markup(varargin)

%Now using u_parkargs interface
args = u_packargs(varargin, 0, ...
            'ResizeFactor', 0.176 ...
			);

% this is the amount by which resolution is reduced after markers have been found (44um to 250um)        
resize_factor = args.ResizeFactor;

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
disp('***************** PECTORAL MUSCLE MARKUP **************************');

%select files to process
%disp('Select the user data files you wish to process');
[filenames, data_path] = uigetfile('J:\*.mat','select user data files to process','Multiselect','on');

%select folder that stores the mammograms
%disp('Select the directory containing the saved user data');
image_path = uigetdir('G:\', 'Select the directory containing the mammograms');
  

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
    
    
    %Generate data and image filenames
    data_filename = [data_path, filesep, filenames{i_file}];
    image_filename = [image_path, filesep, filenames{i_file}(1:end-9), '.tif'];
    
    if ~isempty(strfind(filenames{i_file}(1:end-9), 'ML'));        
        % workout film size - 18x24 and 24x30 are digitised at different orientations
        % 18x24 films will need to be rotated
        if ~isempty(strfind(filenames{i_file}(1:end-9), '1824'))
            filmsizes = 1;
        elseif ~isempty(strfind(filenames{i_file}(1:end-9), '2430'))
            filmsizes = 2;
        else
            %define error action
        end
    
        %is this a left or right breast (search for R not L as L in MLO!)
        left_breast = isempty(strfind(filenames{i_file}(1:end-9), 'R'));
    
        % read in the image
        [IMAGE] = imread(image_filename,'tif');
    
        %load in the existing data
        load(data_filename);
    
        if filmsizes == 1
            % this is only needed for the 1824 film sizes
            IMAGE = imrotate(IMAGE,90);
        end

        % What is the justification for median filtering at this point?!
        % median filter
        IMAGE = medfilt2(IMAGE);

        % reduce image sizes here (don't need such high resolution now magnification markers have been located)
        IMAGE = imresize(IMAGE,resize_factor);
        IMAGE = IMAGE./16;
        IMAGE = 4095-IMAGE;
        if left_breast
            IMAGE = rot90(IMAGE, 2);
        end
        [pectoral_mask pectoral_x pectoral_y] = pectoral_user(IMAGE);
        
        %Save the pectoral masks (less memory to save xy coord of region)
        density_data.pectoral_x = pectoral_x;
        density_data.pectoral_y = pectoral_y;
        save(data_filename, 'density_data');
        
        %Clear density_data and ready for next mammogram
        clear density_data;

        answer = questdlg(...
                'Continue to next image','Mammogram Complete','Yes', 'No', 'Yes');

        if strcmpi(answer, 'no')
            break;
        end
        
    end  
end
helpdlg('Finished processing all selected mammograms','Finished!')
disp('Finished processing all selected mammograms');

%turn warnings back on
warning on Images:initSize:adjustingMag;
warning on Images:imshow:magnificationMustBeFitForDockedFigure
warning on MATLAB:Figure:SetPosition
set(0,'DefaultFigureWindowStyle',orig_window_style);

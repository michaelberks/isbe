function image_id_data = get_image_id_data(varargin)
%GET_IMAGE_ID_DATA *Insert a one line summary here*
%   [] = get_image_id_data(varargin)
%
% GET_IMAGE_ID_DATA uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 05-Dec-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, 0, ...
    'main_xls_list',        'full_image_list_19-12-2013.xlsx',...
    'alternative_xls_lists',{'RSA_new_allocation_observers.xlsx', 'RSA_new_allocation_tonia.xlsx'},...
    'derm_xls_list',        'RSA_dermoscopy_17-10-2014.xlsx',...
    'main_column_order',    [2 3 4 5 6 7],... %Column order = 'id', 'visit', 'hand', 'digit', 'im_name', 'category' 
    'alt_column_orders',    [2 4 5 6 9 3; 1 3 4 5 6 2],... 
    'data_dir',             [nailfoldroot 'data/rsa_study/'],...
    'data_lists',           'data_lists',...
    'image_dir',            'images',...
    'vessel_dir',           'apexes',...
    'contour_dir',          'vessel_contours',...
    'validation_dir',       'test_half/images',...
    'test_dir',             'final_test/images',...
    'markup_dir',           'markup',...
    'marker_list',          [],...
    'save_path',            'image_id_data.mat'); 
clear varargin;

data_list_dir = [args.data_dir args.data_lists '\'];
%%
%--------------------------------------------------------------------------
%Get the main list of images from the full XLS spreadsheet
[~, ~, raw] = ...
    xlsread([data_list_dir args.main_xls_list]);
%num(1,:) = [];
raw(1,:) = [];

num_images = size(raw,1);

image_id_data.im_names = cell(num_images,1);
for i_im = 1:num_images
    image_id_data.im_names{i_im} = raw{i_im, args.main_column_order(5)}(1:6);
end
image_id_data.people_id =  cell2mat( raw(:, args.main_column_order(1)) );
image_id_data.visit =  cell2mat( raw(:, args.main_column_order(2)) );
image_id_data.hand = raw(:, args.main_column_order(3));
image_id_data.digit =  cell2mat( raw(:, args.main_column_order(4)) );
image_id_data.category = raw(:, args.main_column_order(6));
image_id_data.alternative_names = cell(num_images, 1);
image_id_data.derm_names = cell(num_images, 1);
%%
%--------------------------------------------------------------------------
%Now sort out repeats
discard_rows = false(num_images,1);
for id = 1:max(image_id_data.people_id)
    for visit = 1:3
        for hand = {'L', 'R'}
            for digit = 1:5
                im_idx = image_id_data.people_id == id & ...
                    image_id_data.visit == visit & ...
                    strcmpi(image_id_data.hand, hand{1}) & ...
                    image_id_data.digit == digit;
                
                if sum(im_idx) > 1
                    im_idx = find(im_idx);
                    
                    for i_im = 2:length(im_idx)
                        im_name = image_id_data.im_names{im_idx(i_im)};
                        if isempty(image_id_data.alternative_names{im_idx(1)})
                            image_id_data.alternative_names{im_idx(1)} = {im_name};
                        else
                            image_id_data.alternative_names{im_idx(1)}{end+1} = im_name;
                        end
                        discard_rows(im_idx(i_im)) = 1;
                    end
                end
            end
        end
    end
end
image_id_data.im_names(discard_rows,:) = [];
image_id_data.people_id(discard_rows,:) = [];
image_id_data.visit(discard_rows,:) = [];
image_id_data.hand(discard_rows,:) = [];
image_id_data.digit(discard_rows,:) = [];
image_id_data.category(discard_rows,:) = [];
image_id_data.alternative_names(discard_rows,:) = [];        
image_id_data.derm_names(discard_rows,:) = []; 
num_images = length(image_id_data.im_names);
%%
%--------------------------------------------------------------------------
% Now get the alternative names for each image   
for i_alt = 1:length(args.alternative_xls_lists)
    [~, ~, raw] = ...
        xlsread([data_list_dir args.alternative_xls_lists{i_alt}]);
    
    column_order = args.alt_column_orders(i_alt,:);
    for i_im = 2:size(raw,1)

        %Get details for this new image
        id = raw{i_im, column_order(1)};
        visit = raw{i_im, column_order(2)};
        hand = raw{i_im, column_order(3)};
        digit = raw{i_im, column_order(4)};
        im_name = raw{i_im, column_order(5)}(1:end-4);

        %Work out which image it is in th eoriginal list
        im_idx = image_id_data.people_id == id & ...
            image_id_data.visit == visit & ...
            strcmpi(image_id_data.hand, hand) & ...
            image_id_data.digit == digit;

        if ~any(im_idx)
            display(['Warning: no original image found for ' im_name]);
            display(['ID: ' num2str(id) ', visit: ' num2str(visit) ', hand: ' hand ', digit: ' num2str(digit)]);
        elseif sum(im_idx) > 1
            display(['Warning: more than 1 original image found for ' im_name]);
            display(['ID: ' num2str(id) ', visit: ' num2str(visit) ', hand: ' hand ', digit: ' num2str(digit) ' matches images']);
            display(num2str(find(im_idx)));
        elseif isempty(image_id_data.alternative_names{im_idx})
            image_id_data.alternative_names{im_idx} = {im_name};
        else
            image_id_data.alternative_names{im_idx}{end+1} = im_name;
        end
    end
end 

%--------------------------------------------------------------------------
% Now get the dermatoscope image names for each image
[~, ~, raw] = ...
    xlsread([data_list_dir args.derm_xls_list]);  

for i_im = 2:size(raw,1)
    %cols: ID	Visit	Hand	Digit	Filename
    id = raw{i_im,1};
    visit = raw{i_im,2};
    hand = raw{i_im,3};
    digit = raw{i_im,4};
    filename = raw{i_im,5}(2:end);
    
    match = ...
        image_id_data.people_id == id &...
        image_id_data.visit == visit &...
        strcmpi(image_id_data.hand, hand) &...
        image_id_data.digit == digit;
    
    if any(match)
    
        matched_names = image_id_data.derm_names{match};
        
        if ~ismember(filename, matched_names)
            matched_names = [matched_names, {filename}]; %#ok
            image_id_data.derm_names{match} = matched_names;
        end
    else
        display(['No match for ' ...
            'ID ' num2str(id) ...
            ' Visit	' num2str(visit)  ...
            ' Hand	' hand  ...
            ' Digit	' num2str(digit)]);
    end
end    


%%
%--------------------------------------------------------------------------
%Work out who has marked each image...
markup_dir = [args.data_dir args.markup_dir '/'];

%Create a cell array of marker names (that match the folder names in the
%markup dir - it may be easier just to explicitly list these)
if isempty(args.marker_list)
    markers = dir(markup_dir);
    markers = markers(3:end);
    markers_cell = struct2cell(markers);
    markers = markers_cell(1,:)';
    markers(~cell2mat(markers_cell(4,:))) = []; 
else
    markers = args.marker_list;
end
num_markers = length(markers);
image_id_data.marker_list = markers;

image_id_data.markers = cell(num_images,1);
image_id_data.marker_idx = cell(num_images,1);
image_id_data.marker_files = cell(num_images,1);

for i_im = 1:num_images
    for i_ma = 1:num_markers

        im_names = [image_id_data.im_names(i_im) image_id_data.alternative_names{i_im}];
        
        for i_n = 1:length(im_names)
            %Check if marker has marked this image
            markup_list = dir([markup_dir markers{i_ma} '\*#' im_names{i_n} '*_markup.txt']); 
            if isempty(markup_list); 
                continue; 
            end 
        
            %Get the last file in the list (should always be the newest because
            %of the naiming convention in the markup software but we'll sort just in case)
            num_files = length(markup_list);
            marker_files = cell(num_files,1);       
            for i_f = 1:num_files
                marker_files{i_f} = markup_list(i_f).name;
            end
            marker_files = sort(marker_files);
            marker_file = [markup_dir markers{i_ma} '\' marker_files{end}];

            if isempty(image_id_data.markers{i_im})
                image_id_data.markers{i_im} = markers(i_ma);
                image_id_data.marker_idx{i_im} = i_ma;
                image_id_data.marker_files{i_im} = {marker_file};
            else
                image_id_data.markers{i_im}{end+1} = markers{i_ma};
                image_id_data.marker_idx{i_im}(end+1) = i_ma;
                image_id_data.marker_files{i_im}{end+1} = marker_file;
            end
        end
    end
    %if length(image_id_data.markers{i_im}) > length(unique(image_id_data.markers{i_im}))
    %    display(['Image ' image_id_data.im_names{i_im} ' marked multiple by times by same marker']);
    %    display(image_id_data.markers{i_im});
    %end
end

%%
%--------------------------------------------------------------------------
%Set up lists of indices to each subject group (SSc, Primary Raynaud's and
%controls)
image_id_data.category_idx.ss = strcmpi(image_id_data.category, 'S');
image_id_data.category_idx.hc = strcmpi(image_id_data.category, 'HC');
image_id_data.category_idx.pr = strcmpi(image_id_data.category, 'P');
image_id_data.category_idx.u = strcmpi(image_id_data.category, 'U');

image_id_data.category_image_counts.ss = sum(image_id_data.category_idx.ss);
image_id_data.category_image_counts.hc = sum(image_id_data.category_idx.hc);
image_id_data.category_image_counts.pr = sum(image_id_data.category_idx.pr);
image_id_data.category_image_counts.u = sum(image_id_data.category_idx.u);

image_id_data.category_subject_counts.ss = length(unique(image_id_data.people_id(image_id_data.category_idx.ss)));
image_id_data.category_subject_counts.hc = length(unique(image_id_data.people_id(image_id_data.category_idx.hc)));
image_id_data.category_subject_counts.pr = length(unique(image_id_data.people_id(image_id_data.category_idx.pr)));
image_id_data.category_subject_counts.u = length(unique(image_id_data.people_id(image_id_data.category_idx.u)));
%%
%--------------------------------------------------------------------------
%Set up lists of indices to the training, validation and test datasets

%Work out which images were used to train the data
vessel_sizes = {'normal', 'enlarged', 'giant'};
vessel_names = cell(1,3);

for i_sz = 1:3
    vessel_dir = [args.data_dir args.vessel_dir '\' vessel_sizes{i_sz} '\'];
    contour_dir = [args.data_dir args.contour_dir '\' vessel_sizes{i_sz} '\'];
    v_files = dir([contour_dir '*contour.mat']);
    
    num_vessels = length(v_files);
    vessel_names{i_sz} = cell(num_vessels,1); %e   
    
    for i_v = 1:length(v_files)
        %load data
        apex_struc = load([vessel_dir v_files(i_v).name(1:8) '.mat']);       
        [~, im_name] = fileparts(apex_struc.apex_properties.nailfold_name);
        vessel_names{i_sz}{i_v} = im_name;            
    end
end
training_names = unique([vessel_names{1}; vessel_names{2}; vessel_names{3}]);
image_id_data.dataset_idx.training = ismember(image_id_data.im_names, training_names);

%Work out which images were used in the validation set
vali_dir = [args.data_dir args.validation_dir '\'];
vali_list = dir([vali_dir '*.mat']);
validation_names = cell(length(vali_list),1);
for i_im = 1:length(vali_list)
    validation_names{i_im} = vali_list(i_im).name(1:6);
end
image_id_data.dataset_idx.validation = ismember(image_id_data.im_names, validation_names);

%Assign all images that have been markerd and aren't in the other two sets
%to the test set
image_dir = [args.data_dir args.image_dir '\'];
image_id_data.dataset_idx.test = false(num_images, 1);
for i_im = 1:num_images
    if ~image_id_data.dataset_idx.training(i_im) && ... Not in the training set
            ~image_id_data.dataset_idx.validation(i_im) &&... Not in the validation set
            ~isempty(image_id_data.markers{i_im}) &&... Has been marked
            exist([image_dir image_id_data.im_names{i_im} '.png'], 'file') %We've got the actual image
        
        image_id_data.dataset_idx.test(i_im) = 1;
    end
end
test_names = image_id_data.im_names(image_id_data.dataset_idx.test);
%%
%--------------------------------------------------------------------------
%Sanity check, make sure the training, validation and test names all in the
%complete image list
if all(ismember(training_names, image_id_data.im_names))
    display('Training names present and correct!');
else
    display('Hmmm, that''s not right! Missing training names: ');
    display(setdiff(training_names, image_id_data.im_names));
end
if all(ismember(validation_names, image_id_data.im_names))
    display('Validation names present and correct!');
else
    display('Hmmm, that''s not right! Missing training names: ');
    display(setdiff(validation_names, image_id_data.im_names));
end
if all(ismember(test_names, image_id_data.im_names))
    display('Test names present and correct!');
else
    display('Hmmm, that''s not right! Missing test names: ');
    display(setdiff(test_names, image_id_data.im_names));
end
%%
%--------------------------------------------------------------------------
%Display output for how many subjects are in each set
display('*******');
display(['Number of SSc subjects in validation set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.ss & image_id_data.dataset_idx.validation))))]);
display(['Number of Primary Raynaud subjects in validation set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.pr & image_id_data.dataset_idx.validation))))]);
display(['Number of control subjects in validation set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.hc & image_id_data.dataset_idx.validation))))]);
display(['Number of undefined subjects in validation set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.u & image_id_data.dataset_idx.validation))))]);

display('*******');
display(['Number of SSc subjects in training set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.ss & image_id_data.dataset_idx.training))))]);
display(['Number of Primary Raynaud subjects in training set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.pr & image_id_data.dataset_idx.training))))]);
display(['Number of control subjects in training set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.hc & image_id_data.dataset_idx.training))))]);
display(['Number of undefined subjects in training set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.u & image_id_data.dataset_idx.training))))]);

display('*******');
display(['Number of SSc subjects in test set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.ss & image_id_data.dataset_idx.test))))]);
display(['Number of Primary Raynaud subjects in test set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.pr & image_id_data.dataset_idx.test))))]);
display(['Number of control subjects in test set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.hc & image_id_data.dataset_idx.test))))]);
display(['Number of undefined subjects in test set ' ...
    num2str(length(unique(image_id_data.people_id(image_id_data.category_idx.u & image_id_data.dataset_idx.test))))]);

%%
%--------------------------------------------------------------------------
%Display output for how many images are in each set
display('*******');
display(['Number of SSc images in validation set ' ...
    num2str(sum(image_id_data.category_idx.ss & image_id_data.dataset_idx.validation))]);
display(['Number of Primary Raynaud images in validation set ' ...
    num2str(sum(image_id_data.category_idx.pr & image_id_data.dataset_idx.validation))]);
display(['Number of control images in validation set ' ...
    num2str(sum(image_id_data.category_idx.hc & image_id_data.dataset_idx.validation))]);
display(['Number of undefined images in validation set ' ...
    num2str(sum(image_id_data.category_idx.u & image_id_data.dataset_idx.validation))]);

display('*******');
display(['Number of SSc images in training set ' ...
    num2str(sum(image_id_data.category_idx.ss & image_id_data.dataset_idx.training))]);
display(['Number of Primary Raynaud images in training set ' ...
    num2str(sum(image_id_data.category_idx.pr & image_id_data.dataset_idx.training))]);
display(['Number of control images in training set ' ...
    num2str(sum(image_id_data.category_idx.hc & image_id_data.dataset_idx.training))]);
display(['Number of undefined images in training set ' ...
    num2str(sum(image_id_data.category_idx.u & image_id_data.dataset_idx.training))]);

display('*******');
display(['Number of SSc images in test set ' ...
    num2str(sum(image_id_data.category_idx.ss & image_id_data.dataset_idx.test))]);
display(['Number of Primary Raynaud images in test set ' ...
    num2str(sum(image_id_data.category_idx.pr & image_id_data.dataset_idx.test))]);
display(['Number of control images in test set ' ...
    num2str(sum(image_id_data.category_idx.hc & image_id_data.dataset_idx.test))]);
display(['Number of undefined images in test set ' ...
    num2str(sum(image_id_data.category_idx.u & image_id_data.dataset_idx.test))]);

save C:\isbe\nailfold\data\rsa_study\image_id_data.mat image_id_data           

%%
%--------------------------------------------------------------------------
%Finally, save the output
if ~isempty(args.save_path)
    save([data_list_dir args.save_path], 'image_id_data');
end







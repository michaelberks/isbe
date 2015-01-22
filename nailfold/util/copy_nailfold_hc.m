function copy_nailfold_hc(xls_file_in, new_dir, copy_file, xls_lookup_old, xls_file_lookup_new)

%Load in the spreadsheets from the Excel file
[~,~,sheet1]=xlsread(xls_file_in,2);
[~,~,sheet2]=xlsread(xls_file_in,3);

%Workout what rows have been selected in sheet 2 and get these names
selected_idx = strncmpi(sheet2(:,11), 'z', 1);
selected_idx(1) = 0; %Just in case!
selected_names = sheet2(selected_idx, 3:4);

%Pre-allocate cells for the list of names, digits and visit numbers
im_filenames = [];
im_dirs = [];

%Loop through the names
num_names = length(selected_names);

for i_n = 1:num_names
    
    %Now work out what rows match this name in the first sheet and get
    %these images
    selected_images_idx = ...
        strcmpi(sheet1(:,5), selected_names{i_n,1}) & ...
        strcmpi(sheet1(:,6), selected_names{i_n,2});
    
    selected_images = sheet1(selected_images_idx,1);
    selected_dirs = sheet1(selected_images_idx,4);
    
    %Now loop through these images separating the study name, visit number
    %and digit
    num_ims = length(selected_images);
    
    study_names = cell(num_ims, 1);
    visits = cell(num_ims, 1);
    digits = cell(num_ims, 1);
    
    for i_im = 1:num_ims
        idx = regexp(selected_images{i_im}, 'V\d[LR]D\dX\d?', 'once');
        if ~isempty(idx)
            study_names{i_im} = selected_images{i_im}(1:idx-1);
            visits{i_im} = selected_images{i_im}(idx+1);
            digits{i_im} = selected_images{i_im}(idx+(2:4));
        else
            digits{i_im} = num2str(rand);
        end
    end
    
    %Workout the main digit that has been imaged (usually L or R 4)
    [u, ~, ib] = unique(digits);
    main_digit = u{mode(ib)};
    
    %Finally select only the images of the main digit (we'll assume there's
    %a unique one fo reach visit although we could enforce this)
    final_selected_idx = strcmpi(digits, main_digit);
    final_selected_images = selected_images(final_selected_idx);
    final_selected_dirs = selected_dirs(final_selected_idx);
    final_selected_names = study_names(final_selected_idx);
    final_selected_visits = visits(final_selected_idx);
    final_selected_digits = digits(final_selected_idx);
    
    %Loop through the final list of images and write a copy command
    for i_im = 1:length(final_selected_images)
        sub_dir = ['\' final_selected_names{i_im} '\Visit' num2str(final_selected_visits{i_im}) '\' final_selected_digits{i_im}(1)...
            'hand\Digit' final_selected_digits{i_im}(3) '\X300\FullResMosaic'];
        
        final_selected_dirs{i_im} = [final_selected_dirs{i_im} sub_dir];
    end
    
    %Add the filenames and dirs to the master list
    im_filenames = [im_filenames; final_selected_images]; %#ok
    im_dirs = [im_dirs; final_selected_dirs]; %#ok
end

%Generate set of numbers for these images
num_ims = length(im_filenames);

[~,~,prev_names]=xlsread(xls_lookup_old, 5);
num_names = length(prev_names);
used_numbers = zeros(num_names,1);
for i_n = 1:num_names 
    used_numbers(i_n) = str2double(prev_names{i_n,2}(1:5));
end
num_list = setdiff(1:1e5, used_numbers);
num_list = num_list(randperm(length(num_list)));

%Open up file stream
fid = fopen(copy_file, 'wt');
new_names = cell(num_ims,1);
for i_im = 1:num_ims
        
        new_names{i_im} = [zerostr(num_list(i_im),5) '.png'];
        copy_line = [
            'copy "' ...
            im_dirs{i_im} '\' im_filenames{i_im} ...
            '" ' new_dir new_names{i_im}];
        
        
        fprintf(fid, '%s \n', copy_line);
end
fclose(fid);

%Write out the matching list to an xls file
xlswrite(xls_file_lookup_new, {'Original image name', 'New random image name'}, 1, 'A1');
xlswrite(xls_file_lookup_new, [im_filenames new_names], 1, 'A2');


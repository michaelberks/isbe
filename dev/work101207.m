%Get path of all sub-directories
dir_path = genpath('C:\isbe\asymmetry_project');

%Workout where each sub-directory starts and ends
dir_ends = strfind(dir_path, pathsep);
num_dirs = length(dir_ends);

% Loop through the sub-directories looking for dicom files
id2 = -1;

for ii = 1:num_dirs

    id1 = id2 + 2;
    id2 = dir_ends(ii) - 1;
    dir_name = [dir_path(id1:id2) filesep];
    o04_list = dir([dir_name 'o04_*.mat']);

    %For each dicom file found, try converting it's name and moving to the
    %output directory
    for jj = 1:length(o04_list)
        movefile([dir_name o04_list(jj).name], [dir_name o04_list(jj).name(5:end)]);
        display(['Moved ' dir_name o04_list(jj).name ' to ' o04_list(jj).name(5:end)]);
    end
end
%%
data_type = {'abnormals', 'normals'};
    
for ii = 1:2
    mammo_dir = ['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\' data_type{ii} '\'];

    mammo_list = dir([mammo_dir '*ML.mat']);
    mammo_names = get_mammo_info(mammo_list);
    
    for jj = 1:2;
        
        %load segmentation
        segmentation = u_load(['C:\isbe\asymmetry_project\data\segmentations\2004_screening\' data_type{ii} '\'...
            mammo_names{jj} '_segmentation.mat']);
        
        %load K-line map
        line_map = load(['C:\isbe\asymmetry_project\data\k_line_maps\2004_screening_processed\' data_type{ii} '\'...
            mammo_names{jj} '_data.mat'], 'line_map');
        [r c] = size(line_map);
        
        %resize the segmentation to the map
        [breast_border breast_air pectoral_edge] = segment_breast_resize([r c], segmentation);

        %Make a mask of the pectoral triangle
        if strcmpi(mammo_names{jj}(4), 'R')
            sx = c;
            sy = 1;
        else
            sx = 1;
            sy = 1;
        end
        mask = poly2mask([sx; pectoral_edge(:,1)], [sy; pectoral_edge(:,2)], r, c);

        %Set the line map to zero inside the mask
        line_map(mask) = 0;
        
        figure('Name', mammo_names{jj});
        imagesc(line_map); axis image; colormap(gray(256));
        clear line_map;
        
        %Repeat for the RF maps
        %load K-line map
        line_map = load_uint8(['C:\isbe\asymmetry_project\data\k_line_maps\2004_screening_processed\' data_type{ii} '\'...
            mammo_names{jj} '_class.mat']);
        [r c] = size(line_map);
        
        %resize the segmentation to the map
        [breast_border breast_air pectoral_edge] = segment_breast_resize([r c], segmentation);

        %Make a mask of the pectoral triangle
        if strcmpi(mammo_names{jj}(4), 'R')
            sx = c;
            sy = 1;
        else
            sx = 1;
            sy = 1;
        end
        mask = poly2mask([sx; pectoral_edge(:,1)], [sy; pectoral_edge(:,2)], r, c);

        %Set the line map to zero inside the mask
        line_map(mask) = 0;
        
        figure('Name', mammo_names{jj});
        imagesc(line_map); axis image; colormap(gray(256));
        clear line_map;
    end
end
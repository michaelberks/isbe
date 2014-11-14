function mb_quick_manual_segmentation(mammo_dir, start_num, file_type)

if nargin < 2
    start_num = 1;
end
if nargin < 3
    file_type = '.mat';
end

mammo_list = dir([mammo_dir, '*', file_type]);

seg_fig = figure;
for ii = start_num:length(mammo_list)
    image_name = [mammo_dir, mammo_list(ii).name];
    
    if strcmp(file_type, '.mat')
        i1 = u_load(image_name);
    else
        i1 = imread(image_name);
    end
    figure(seg_fig); 
    breast_region = roipoly(i1); %#ok
    set(seg_fig, 'Name', ['Mass ', num2str(ii), ' of ', num2str(length(mammo_list)),...
        ': ', mammo_list(ii).name]); 
    save_name = [image_name(1:end-4), '_mask.mat'];
    save(save_name, 'breast_region');
    clear i1 breast_region image_name save_name
    
    
    set(0,'DefaultFigureWindowStyle','normal');
    selection = questdlg('Continue to next image?',...
                     'Warning',...
                     'Yes','No','Yes');
    set(0,'DefaultFigureWindowStyle','docked');             
    if strcmpi(selection, 'No')
        break
    end
    
end

close(seg_fig);
display(['Manual segmentation terminated at image ', num2str(ii)]);
            
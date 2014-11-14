function [mammo_names stepwedge_timings] = stepwedge_timings_to_excel(markup_dir, excel_file)

if nargin < 2 || isempty(excel_file);
    excel_file = 'stepwedge_timings.xls';
end
if nargin < 1 || isempty(markup_dir);
    try
        fid = fopen('C:\program files\density_software\default_locations.txt');
        def_txt = textscan(fid, '%s %s', 'delimiter', '=');
        fclose(fid);

        markup_dir = ...
            def_txt{2}{strncmpi(def_txt{1}, 'markup_dir', length('markup_dir'))};

    catch %#ok
        display('Problem opening default file locations');
        markup_dir = ...
            uigetdir([], 'Select the directory to save markup data to');
    end
end

if ~strcmp(markup_dir(end), filesep)
    markup_dir = [markup_dir, filesep];
end

markup_list = dir([markup_dir, '*.mat']);
num_markups = length(markup_list);

stepwedge_timings = zeros(num_markups, 7);
start_and_end_times = cell(num_markups, 2);
mammo_names = cell(num_markups, 1);

for ii = 1:num_markups
    
    density_data = u_load([markup_dir, markup_list(ii).name]);
    
    mammo_names{ii} = correct_name(markup_list(ii).name(1:end-12), 4);
    
    if isfield(density_data, 'load_time'); stepwedge_timings(ii, 1)                 = density_data.load_time; end
    if isfield(density_data, 'segmentation_time'); stepwedge_timings(ii, 2)         = density_data.segmentation_time; end
    if isfield(density_data, 'pectoral_time'); stepwedge_timings(ii, 3)             = density_data.pectoral_time; end
    if isfield(density_data, 'marker_markup_time'); stepwedge_timings(ii, 4)        = density_data.marker_markup_time; end
    if isfield(density_data, 'marker_processing_time'); stepwedge_timings(ii, 5)    = density_data.marker_processing_time; end
    if isfield(density_data, 'stepwedge_markup_time'); stepwedge_timings(ii, 6)     = density_data.stepwedge_markup_time; end
    if isfield(density_data, 'stepwedge_processing_time'); stepwedge_timings(ii, 7) = density_data.stepwedge_processing_time; end
    if isfield(density_data, 'gland_processing_time'); stepwedge_timings(ii, 8)     = density_data.gland_processing_time; end

    if isfield(density_data, 'time_started'); start_and_end_times{ii, 1}            = density_data.time_started; end
    if isfield(density_data, 'time_finished'); start_and_end_times{ii, 2}           = density_data.time_finished; end
   
end

   
%Check whether the file exists, and if so, whether to add to the
%existing data or overwrite it
if exist([markup_dir excel_file], 'file')
    num = xlsread([markup_dir excel_file]);
    if ~isempty(num)        
        headings = false;
        start_row = num2str(size(num,1)+2);
    else
        start_row = '2';
        headings = true;
    end
else
    start_row = '2';
    headings = true;
end
if headings
    xlswrite([markup_dir excel_file], {'Name'}, 1, 'A1');
    xlswrite([markup_dir excel_file], {'Load time'}, 1, 'B1');
    xlswrite([markup_dir excel_file], {'Segmentation'}, 1, 'C1');
    xlswrite([markup_dir excel_file], {'Pectoral'}, 1, 'D1');
    xlswrite([markup_dir excel_file], {'Marker markup'}, 1, 'E1');
    xlswrite([markup_dir excel_file], {'Marker processing'}, 1, 'F1');
    xlswrite([markup_dir excel_file], {'Stepwedge markup'}, 1, 'G1');
    xlswrite([markup_dir excel_file], {'Stepwedge processing'}, 1, 'H1');
    xlswrite([markup_dir excel_file], {'Gland processing'}, 1, 'I1');
    xlswrite([markup_dir excel_file], {'Start time'}, 1, 'J1');
    xlswrite([markup_dir excel_file], {'End time'}, 1, 'K1');
end

xlswrite([markup_dir excel_file], mammo_names, 1, ['A', start_row]);
xlswrite([markup_dir excel_file], stepwedge_timings, 1, ['B', start_row]);
xlswrite([markup_dir excel_file], start_and_end_times, 1, ['J', start_row]);

function corrected_name = correct_name(name, len)
    
    views = {'LCC', 'RCC', 'LML', 'RML'}; 
    pos = []; ii = 1; 
    while isempty(pos); pos = strfind(name,views{ii}); ii = ii+1; end
    corrected_name = name;
    for ii = 1:len - pos + 1
        corrected_name = ['0', corrected_name]; %#ok
    end
    
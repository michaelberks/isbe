function [failures error_codes unknown] = stepwedge_failures(varargin)

%Now using u_parkargs interface
args = u_packargs(varargin, 0, ...
            'DataPath', [], ...
            'ExcelFile', [] ...
			);          

if isempty(args.DataPath)
        %select files to process
        %disp('Select the user data files you wish to process');
        [filenames, data_path] = uigetfile('J:\*.mat','select user data files to check','Multiselect','on');
else
    data_path = args.DataPath;
    temp = dir([data_path, filesep, '*data*.mat*']);
        filenames = cell(length(temp),1);
        for ii = 1:length(temp); filenames{ii} = temp(ii).name; end
        clear temp;   

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
error_codes = [];
users = [];
unknown = [];

n_process = length(filenames);
%%
for i_file=1:n_process  

    data_filename = [data_path, filesep, filenames{i_file}];
    
    %Don't crash if data file doesn't exist
    try
        %load in the existing data
        density_data = u_load(data_filename);
    catch
        disp(['Error loading data file: ', data_filename, '. Check file is of correct format']);
        unknown{end+1,1} = data_filename; %#ok
        continue;       
    end
    
    %Check marker pair data exists (and hence errorcheck is neg)
    error_code = do_errorcheck(density_data);
    if ~isempty(error_code);
        display(['Incorrect user_data file: ', data_filename, '. Skipping this mammogram']);
        failures{end+1,1} = correct_name(filenames{i_file}(1:end-9), 4); %#ok
        error_codes{end+1,1} = error_code; %#ok
        users{end+1,1} = density_data.user; %#ok
    end
    
end

if ~isempty(args.ExcelFile)
    xlswrite([data_path, filesep, args.ExcelFile], {'Name'}, 1, 'A1');
    xlswrite([data_path, filesep, args.ExcelFile], {'Error Code'}, 1, 'B1');
    xlswrite([data_path, filesep, args.ExcelFile], {'User'}, 1, 'C1');
    if ~isempty(failures);
        xlswrite([data_path, filesep, args.ExcelFile], failures, 1, 'A2');
        xlswrite([data_path, filesep, args.ExcelFile], error_codes, 1, 'B2');
        xlswrite([data_path, filesep, args.ExcelFile], users, 1, 'C2');
    end

end

end
%--------------------------------------------------------------------------
function data_error = do_errorcheck(density_data)
    
    data_error = [];
    if ~isfield(density_data, 'x_b_info')
        data_error = 'No marker points';
    elseif ~isfield(density_data, 'nipple_position')
        data_error = 'No nipple selected';
    elseif ~isfield(density_data, 'wedgevals') || ...
            ~isfield(density_data, 'swx') || ...
            ~isfield(density_data, 'swy')
        data_error = 'Error selecting stepwedge';
    elseif ~isfield(density_data, 'coarse_edgex') || ...
            ~isfield(density_data, 'coarse_edgey')
        data_error = 'Error processing stepwedge or marking breast boundary';
    end
    
end
%--------------------------------------------------------------------------
function corrected_name = correct_name(name, len)
    
    views = {'LCC', 'RCC', 'LML', 'RML'}; 
    pos = []; ii = 1; 
    while isempty(pos); pos = strfind(name,views{ii}); ii = ii+1; end
    corrected_name = name;
    for ii = 1:len - pos + 1
        corrected_name = ['0', corrected_name]; %#ok
    end
end

    


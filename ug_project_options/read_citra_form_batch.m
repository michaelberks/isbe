function [density_percents] = ...
    read_chitra_form_batch(form_name, num_forms, excel_file, reader_initials, debug_flag)
%READ_DENSITY_FORM_BATCH Load in the set of density forms for a given
%reader and call read_density_form to compute density percentages
%
% Inputs:
%      reader_dir - scanned image of density form
%
%      num_forms - the number of forms expected in the folder
%
%      debug_flag- flag {0,1} to turn debug_flag output on or off
%
%
% Outputs:
%
%      density_percents - Percentages read from each of the 4 line-scales
%
%
%
% Example:
%
% Notes:
%
% See also: READ_DENSITY_FORM
%
% Created: 11-Jun-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% Set the default for the debug parameter to 0
if nargin < 1
    form_name = [];
end

%Open a UI box to select the form pdf if not already specified
if isempty(form_name)
    [form_name, form_dir] = uigetfile('*.pdf','Select the form to process','Multiselect','off');
    form_name = form_name(1:end-4);
else
    [form_dir, form_name] = fileparts(form_name);
end
if ~strcmp(form_dir(end), filesep)
    form_dir = [form_dir, filesep];
end

%Now make all the tiff files for the form using ghostscript
cmd = find_ghostscript;
options = [' -dBATCH -q -dNOPAUSE -sDEVICE=tiffg4 -r300 '...
    '-sOutputFile="' form_dir form_name, '_%03d.tif" "' form_dir form_name '.pdf" -c quit'];
[status result] = system([cmd options]);

if status
    %something has gone wrong!
    display(result);
    error('Error generating tiff files');
end

%Get the list of tiff files associated with the form PDF file
form_list = dir([form_dir, form_name, '*.tif']);

% Set the default number of forms to the number of matching tif files found
if nargin < 2
    num_forms = [];
end
if isempty(num_forms)
    num_forms = length(form_list);
end

% Set the default for the execl file
if nargin < 3
    excel_file = [];
end

% Set the default for the reader initials
if nargin < 4
    reader_initials = [];
end

% Set the default for the debug parameter to 0
if nargin < 5
    debug_flag = 0;
end

%pre-allocate space for the density %s of each form
density_percents = zeros(num_forms, 4);
nhs_ids = zeros(num_forms, 1);
dobs = cell(num_forms, 1);
reader_initials_all = cell(num_forms, 1);
warnings = cell(num_forms, 1);
timings = cell(num_forms, 1);

% Go through each form...
for ii = 1:num_forms
    
    %Load form
    form = imread([form_dir, form_list(ii).name]);
    
    %Read percentages, nhs number and DOB from form
    [dp nhs_id dob debug_out] = read_density_form(form, debug_flag);
    
    %Save data into structures
    density_percents(ii, :) = dp;
    nhs_ids(ii) = nhs_id;
    dobs{ii} = dob;
    reader_initials_all{ii} = reader_initials;
    
    %Check if we received any warnings reading the form
    warnings{ii} = [];
    
    %Were the page markers correctly found
    if debug_out.warning
        display(['Warning in form ', num2str(ii)]);
        warnings{ii} = [warnings{ii} 'Form may not have been aligned correctly; '];
    end
    
    %Were any density readings missing
    if any(debug_out.no_mark)
        display(['Warning in form ', num2str(ii)]);
        warnings{ii} = [warnings{ii} 'No density reading detected for:'];
        if debug_out.no_mark(1)
            warnings{ii} = [warnings{ii}, ' RCC;'];
        end
        if debug_out.no_mark(2)
            warnings{ii} = [warnings{ii}, ' RML;'];
        end
        if debug_out.no_mark(3)
            warnings{ii} = [warnings{ii}, ' LCC;'];
        end
        if debug_out.no_mark(4)
            warnings{ii} = [warnings{ii}, ' LML;'];
        end
    end
    timings{ii} = datestr(now);
end

%If we've been given an excel file, write data to it
if ~isempty(excel_file);
    
    %Check whether the file exists, and if so, whether to add to the
    %existing data or overwrite it
    if exist(excel_file, 'file')
        num = xlsread(excel_file);
        
        if ~isempty(num)        
%             answer = questdlg(...
%                 'Add density readings to existing data or overwrite?',...
%                 'Excel file contains data','Add', 'Overwrite', 'Add');
%             if strcmpi(answer, 'add')
%                 start_row = num2str(size(num,1)+2);
%                 headings = false;
%             else
%                 answer = questdlg(...
%                 'Are you sure you want to overwrite the existing data?',...
%                 'Excel file contains data','Yes', 'No', 'No');
%                 if strcmpi(answer, 'yes')
%                     start_row = '2';
%                     headings = true;
%                 else
%                     headings = false;
%                     start_row = num2str(size(num,1)+2);
%                 end
%             end
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
        xlswrite(excel_file, {'NHS Number'}, 1, 'A1');
        xlswrite(excel_file, {'D.O.B.'}, 1, 'B1');
        xlswrite(excel_file, {'RCC'}, 1, 'C1');
        xlswrite(excel_file, {'RMLO'}, 1, 'D1');
        xlswrite(excel_file, {'LCC'}, 1, 'E1');
        xlswrite(excel_file, {'LMLO'}, 1, 'F1');
        xlswrite(excel_file, {'Reader initials'}, 1, 'G1');
        xlswrite(excel_file, {'Warnings'}, 1, 'H1');
        xlswrite(excel_file, {'Time read'}, 1, 'I1');
    end
    
    xlswrite(excel_file, nhs_ids, 1, ['A', start_row]);
    xlswrite(excel_file, dobs, 1, ['B', start_row]);
    xlswrite(excel_file, density_percents, 1, ['C', start_row]);
    xlswrite(excel_file, reader_initials_all, 1, ['G', start_row]);
    xlswrite(excel_file, warnings, 1, ['H', start_row]);
    xlswrite(excel_file, timings, 1, ['I', start_row]);
end

%Now ask the user if they want to delete the form .tiff files
answer = questdlg(...
            ['Do you want to delete the .tiff files associated with' ...
            ' the forms (the PDF file will still be saved)'],...
            'Density reading complete','Yes', 'No', 'Yes');
if strcmpi(answer, 'yes')
    for ii = 1:num_forms
        delete([form_dir, form_list(ii).name]);
    end
end

function cmd = find_ghostscript
% Find the full path to a ghostscript executable
cmd = '';
if ispc
    % For Windows, look in the default location
    default_location = 'C:\Program Files\gs\';
    executable = '\bin\gswin32c.exe';
else
    % This case isn't supported. Contact me if you have a fix.
    return
end
dir_list = dir(default_location);
ver_num = 0;
for a = 1:numel(dir_list)
    % If there are multiple versions, use the newest
    ver_num2 = sscanf(dir_list(a).name, 'gs%g');
    if ~isempty(ver_num2) && ver_num2 > ver_num
        cmd2 = [default_location dir_list(a).name executable];
        if exist(cmd2, 'file') == 2
            cmd = ['"' cmd2 '"'];
            ver_num = ver_num2;
        end
    end
end
return
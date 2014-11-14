function [density_percents_1 density_percents_2] = ...
    read_chitra_form_batch(form_name, num_forms, excel_file, debug_flag)
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

% Set the default for the debug parameter to 0
if nargin < 4
    debug_flag = 0;
end

%pre-allocate space for the density %s of each form
density_percents_1 = zeros(100, 1, 15);
density_percents_2 = zeros(50, 2, 15);
readers_1 = false(15,1);
readers_2 = false(15,1);

% Go through each form...
for ii = 1:num_forms
    
    %Load form
    form = imread([form_dir, form_list(ii).name]);
    
    %Read percentages, nhs number and DOB from form
    [dp form_type form_idx reader_id debug_out] = read_chitra_form(form, debug_flag);
    
    %Save data into structures
    rows = 5*(form_idx-1) + (1:5);
    if form_type == 1
        density_percents_1(rows, :, reader_id) = dp;
        readers_1(reader_id) = 1;
        
        if ~isempty(excel_file);
            xls_range = ['B' num2str(rows(1)+1) ':B' num2str(rows(5)+1)];
            xlswrite(excel_file, {'Image order'}, reader_id, 'A1');
            xlswrite(excel_file, {'Density'}, reader_id, 'B1');
            
            xlswrite(excel_file, dp, reader_id, xls_range);
        end
            
    elseif form_type == 4
        form_type = form_type - 2;
        density_percents_2(rows, :, reader_id) = reshape(dp, 2, 5)';
        readers_2(reader_id) = 1;
        
        if ~isempty(excel_file);
            xls_range = ['E' num2str(rows(1)+1) ':F' num2str(rows(5)+1)];
            xlswrite(excel_file, {'Image order'}, reader_id, 'D1');
            xlswrite(excel_file, {'Density A'}, reader_id, 'E1');
            xlswrite(excel_file, {'Density B'}, reader_id, 'F1');
            
            xlswrite(excel_file, reshape(dp, 2, 5)', reader_id, xls_range);
        end
    else
        error('Bugger. Wrong form type.');
    end
    
    %Check for any warnings
    ww = [];
    
    %Were the page markers correctly found
    if debug_out.warning
        ww = [ww 'Form may not have been aligned correctly; ']; %#ok
    end
    
    %Were any density readings missing
    if any(debug_out.no_mark)
        ww = [ww 'No density reading detected for VAS line: ']; %#ok
        
        if form_type == 1
            num_lines = 5;
        else
            num_lines = 10;
        end
        for jj = 1:num_lines
            if debug_out.no_mark(1) 
                ww = [ww num2str(jj) '; ']; %#ok
            end
        end
    end
    if ~isempty(ww)
        ww = ['Warning in form ' num2str(form_idx) ' model ' num2str(form_type) ' for reader ' num2str(reader_id) ww]; %#ok
        display(ww);
    end

end

%If we've been given an excel file, write data to it
% if ~isempty(excel_file);
%     
%     %Check what we've read in data for
%     for ii = 1:15
%         if readers_1(ii)
%             xlswrite(excel_file, {'Image order'}, ii, 'A1');
%             xlswrite(excel_file, {'Density'}, ii, 'B1');
%             
%             xlswrite(excel_file, density_percents_1(:,:,ii), ii, 'B2:B101');
%         end
%         if readers_2(ii)
%             xlswrite(excel_file, {'Image order'}, ii, 'D1');
%             xlswrite(excel_file, {'Density A'}, ii, 'E1');
%             xlswrite(excel_file, {'Density B'}, ii, 'F1');
%             
%             xlswrite(excel_file, density_percents_2(:,:,ii), ii, 'E2:F51');
%         end
%     end
% end

%Now ask the user if they want to delete the form .tiff files
% answer = questdlg(...
%             ['Do you want to delete the .tiff files associated with' ...
%             ' the forms (the PDF file will still be saved)'],...
%             'Density reading complete','Yes', 'No', 'Yes');
% if strcmpi(answer, 'yes')
for ii = 1:num_forms
    delete([form_dir, form_list(ii).name]);
end
% end

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
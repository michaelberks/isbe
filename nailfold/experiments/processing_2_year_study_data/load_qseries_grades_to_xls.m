function [] = load_qseries_grades_to_xls(grade_dir, xls_filename, subject_ids, top_older)
%LOAD_QSERIES_GRADES_TO_XLS *Insert a one line summary here*
%   [] = load_qseries_grades_to_xls(grade_dir, xls_filename)
%
% Inputs:
%      grade_dir - *Insert description of input variable here*
%
%      xls_filename - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 09-Apr-2014
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
num_reasons = 13;

if exist(xls_filename, 'file')
    num = xlsread(xls_filename);
        
    if ~isempty(num)        
        answer = questdlg(...
            'Add grade readings to existing data or overwrite?',...
            'Excel file contains data','Add', 'Overwrite', 'Add');
        if strcmpi(answer, 'add')
            start_row = size(num,1)+2;
            do_headings = false;
        else
            answer = questdlg(...
            'Are you sure you want to overwrite the existing data?',...
            'Excel file contains data','Yes', 'No', 'No');
            if strcmpi(answer, 'yes')
                start_row = '2';
                do_headings = true;
            else
                do_headings = false;
                start_row = size(num,1)+2;
            end
        end
    else
        start_row = 2;
        do_headings = true;
    end
else
    start_row = 2;
    do_headings = true;
end
    

if do_headings
    headings1 = {
        'Subject ID' ...
        'Top image' ...
        'Images show disease' ...
        'Progression'};
    
    headingsr = cell(1,num_reasons);
    for i_r = 1:num_reasons
        headingsr{i_r} = ['Reason ' num2str(i_r)];
    end
    
    headings2 = {
        'Observer' ...
        'Timestamp' ...
        'Time taken'...
        'Display order'};
    
    headings = [headings1 headingsr headings2];
        
    xlswrite(xls_filename, headings, 1, 'A1'); 
end

grades_list = dir([grade_dir '\*grades.txt']);
num_grades = length(grades_list);

xls_data = cell(num_grades, length(headings));
curr_row = 1;
for i_g = 1:num_grades
    
    grade_name = grades_list(i_g).name;
    grade = load_qseries_grade_file([grade_dir '\' grade_name]);
    
    V_idx = find(grade_name == 'V', 1);
    if isempty(V_idx);
        display(['Error, unexpected filename, skipping: ' grade_name]);
        continue;
    end
        
    curr_col = 1;
    
    %1: Do subject ID
    subject_id = grade_name(1:V_idx-1);   
    xls_data{curr_row, curr_col} = subject_id;
    curr_col = curr_col+1;
    %-------------------
    
    %2: Mark which was the top image
    s_idx = find(strcmpi(subject_id, subject_ids)); 
    if top_older(s_idx)
        top_image = 'V1';
    else
        top_image = 'V6';
    end   
    xls_data{curr_row, curr_col} = top_image;
    curr_col = curr_col+1;
    %-------------------
    
    %3: Were the images normal
    xls_data{curr_row, curr_col} = double(grade.images_abnormal);
    curr_col = curr_col+1;
    %-------------------
    
    %4: If the images were abnormal, write progression
    if grade.images_abnormal
        %If top image was V1, make progression negative
        if top_older(s_idx)
            grade.progression = -1*grade.progression;
        end
        xls_data{curr_row, curr_col} = grade.progression;
    end
    curr_col = curr_col+1;
    %-------------------
    
    %5: If progression, mark the reasons
    for i_r = 1:num_reasons
        if grade.images_abnormal && grade.progression            
            xls_data{curr_row, curr_col} = double(grade.reasons(i_r));
        end
        curr_col = curr_col+1;
    end
    %-------------------
            
    %6: Observer
    xls_data{curr_row, curr_col} = grade.observer;
    curr_col = char(curr_col+1);
    %-------------------
    
    %7: Timestamp
    xls_data{curr_row, curr_col} = grade.timestamp;
    curr_col = curr_col+1;
    %-------------------
    
    %8: Time taken
    xls_data{curr_row, curr_col} = grade.time_taken;
    curr_col = curr_col+1;
    %-------------------
    
    %9: Display order
    xls_data{curr_row, curr_col} = s_idx;
    %-------------------
    
    display(['Successfully loaded ' grade_name ' to xls data']);
    curr_row = curr_row + 1;
end

xlswrite(xls_filename, xls_data, 1, ['A' num2str(start_row)]);

%     curr_col = 'A';
%     
%     %1: Do subject ID
%     subject_id = {grade_name(1:V_idx-1)};   
%     xlswrite(xls_filename, subject_id, 1, [curr_col num2str(curr_row)]);
%     curr_col = char(curr_col+1);
%     %-------------------
%     
%     %2: Mark which was the top image
%     s_idx = find(strcmpi(subject_id, subject_ids)); 
%     if top_older(s_idx)
%         top_image = {'V1'};
%     else
%         top_image = {'V6'};
%     end   
%     xlswrite(xls_filename, top_image, 1, [curr_col num2str(curr_row)]);
%     curr_col = char(curr_col+1);
%     %-------------------
%     
%     %3: Were the images normal
%     xlswrite(xls_filename, double(grade.images_abnormal), 1, [curr_col num2str(curr_row)]);
%     curr_col = char(curr_col+1);
%     %-------------------
%     
%     %4: If the images were abnormal, write progression
%     if grade.images_abnormal
%         %If top image was V1, make progression negative
%         if top_older(s_idx)
%             grade.progression = -1*grade.progression;
%         end
%         xlswrite(xls_filename, grade.progression, 1, [curr_col num2str(curr_row)]);
%     end
%     curr_col = char(curr_col+1);
%     %-------------------
%     
%     %5: If progression, mark the reasons
%     if grade.images_abnormal && grade.progression
%         xlswrite(xls_filename, double(grade.reasons'), 1, [curr_col num2str(curr_row)]);
%     end
%     curr_col = char(curr_col+num_reasons);
%     %-------------------
%             
%     %6: Observer
%     xlswrite(xls_filename, {grade.observer}, 1, [curr_col num2str(curr_row)]);
%     curr_col = char(curr_col+1);
%     %-------------------
%     
%     %7: Timestamp
%     try
%         xlswrite(xls_filename, cellstr(grade.timestamp), 1, [curr_col num2str(curr_row)]);
%     catch %#ok
%         display('Coudln''t write time stamp');
%     end
%     curr_col = char(curr_col+1);
%     %-------------------
%     
%     %8: Time taken
%     xlswrite(xls_filename, grade.time_taken, 1, [curr_col num2str(curr_row)]);
%     curr_col = char(curr_col+1);
%     %-------------------
%     
%     %9: Display order
%     xlswrite(xls_filename, s_idx, 1, [curr_col num2str(curr_row)]);
%     %-------------------
%     
%     display(['Successfully written ' grade_name ' to ' xls_filename]);
%     curr_row = curr_row + 1;

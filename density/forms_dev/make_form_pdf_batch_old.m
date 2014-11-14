function make_form_pdf_batch_old(xls_file, form_folder, sheet)

if nargin < 3
    sheet = 1;
end

%Read in the data from the excel spreadsheet
[dummy, patient_data] = xlsread(xls_file, sheet);

%Need to extract name, nhs_id and dob from each cell - we need to be careful
%because we can't assume the exact form the data will be in, other than 
% i) it will be in the order name, nhs_id, dob
% ii) there will be spaces between name-nhs_id and nhs_id-dob
for pp = 1:length(patient_data)
    
    %Find the spaces in the data
    spaces = find(patient_data{pp} == ' ');
    
    %Take the first name as the string before the first space
    first_name = patient_data{pp}(1:spaces(1)-1);
    
    %Similarly the second string
    second_name = patient_data{pp}(spaces(1)+1:spaces(2)-1);
    
    %Take the dob as the string following the final space and check this is
    %a valid dob
    dob = patient_data{pp}(spaces(end)+1:end);
    if valid_date(dob)

    else
        display('DOB is not valid. Try something else');
    end
    
    %Now for all the stuff in between, we first need to make sure there are
    %no more names, we then assume what remains is the patient nhs_id, which
    %may or may not have spaces in it
    nhs_id = [];
    for jj = 2:length(spaces)-1
        
        %Check whether this portion is a number
        if isempty(str2num(patient_data{pp}(spaces(jj)+1:spaces(jj+1)-1))) %#ok
            %Portion is not a number so append to second name (include the
            % initial space)
            second_name = [second_name, patient_data{pp}(spaces(jj):spaces(jj+1)-1)];%#ok
        else
            %Is a number so append string to number (exclude the space)
            nhs_id = [nhs_id patient_data{pp}(spaces(jj)+1:spaces(jj+1)-1)];%#ok
        end
    end
    %Convert the nhs_id string into a number
    nhs_id = str2num(nhs_id); %#ok
    
    if isempty(nhs_id)
        %We have a problem with the nhs_id
        display('ID is not valid');
    else
        display(['Name: ', second_name, ', ', first_name, ' ID: ', num2str(nhs_id), ' DOB: ', dob]);
        make_form_pdf(second_name, first_name, nhs_id, dob, form_folder)
    end
end


function date_ok = valid_date(date_string)

    date_ok =...
        (length(date_string) == 8 && date_string(3) == '/' && date_string(6) == '/') ||...
        (length(date_string) == 10 && date_string(3) == '/' && date_string(6) == '/');
        

    
    
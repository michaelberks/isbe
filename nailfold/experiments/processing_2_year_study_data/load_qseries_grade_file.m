function [grade] = load_qseries_grade_file(filename)
%LOAD_QSERIES_GRADE_FILE *Insert a one line summary here*
%   [grade] = load_qseries_grade_file(filename)
%
% Inputs:
%      filename - *Insert description of input variable here*
%
%
% Outputs:
%      grade - *Insert description of input variable here*
%
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


fid = fopen(filename);
while 1
    tline = fgets(fid);
    if ~ischar(tline)
        break;    
    else
        start_idx = find(tline ~= ' ', 1);
        
        if isempty(start_idx)
            display('Empty line found')
            continue;
        end
        tline = tline(start_idx:end - 2);   
        
        if tline(1) == '}'
            break;
        end
        
        colon_idx = find(tline == ':', 1);
        
        if isempty(colon_idx) %This should never happen but just in case
            display(['Unexpected line found: ' tline]);
            continue;
        end
        
        tag = tline(1:colon_idx -1);
        value = tline(colon_idx + 2:end);
        %display(value);
        
        switch tag
            
            case 'ncm_qseries_grade'
                %Do nothing, its just the title
                
            case 'version'
                grade.version = str2double(value);
                
            case 'observer'
                grade.observer = value;
                
            case 'timestamp'
                grade.timestamp = value;
                
            case 'time_taken'
                grade.time_taken = str2double(value);
                
            case 'ims_gradeable'
                grade.images_gradeable = value(1) == '1';
                
            case 'ims_normal'
                grade.images_abnormal = value(1) == '0';
                
                if ~grade.images_abnormal
                    %We can break as no further info needs collecting
                    break;
                end
                
            case 'progression_level'
                grade.progression = str2double(value);
                
                if ~grade.progression
                    %If no progression we can break as no reasons to give
                    break
                end
                %Othwerwise pre-allocate reasons
                grade.reasons = false(13,1);
                
            otherwise
                %Should be a reason now, but check just in case
                if tag(1) ~= 'r'
                    display(['Unexepcted tag found: ' tline]);
                    
                else
                    %Get the reason number
                    
                    %Try a digit number
                    r_idx = str2double(tag(2:3));
                    if isnan(r_idx)
                        r_idx = str2double(tag(2));
                    end
                    
                    grade.reasons(r_idx+1) = value(1) == '1';
                end
        end                                          
    end
end
fclose(fid);

%Some old grade file may not have the gradeable tag, so if the field hasn't
%been set, we just have to assume the images were gradeable
if ~isfield(grade, 'images_gradeable')
    grade.images_gradeable = 1;
end
        


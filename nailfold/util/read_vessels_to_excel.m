function [] = read_vessels_to_excel(vessel_data_file, excel_file, sheet_name)

% f_debug = (nargin == 0 && nargout == 0);
% if (f_debug)
% 	vessel_data_file = 'U:\projects\nailfold\images\tmp.txt';
% end
% 
% fid = fopen(vessel_data_file);

% Get header information
% %find observer
% loop = true;
% while loop
%     observer = fgetl(fid);
%     loop = isempty(strfind(observer, 'observer'));
% end
    
% %loop through rest of file
% vessel_num = 0;
% total_apices = 0;
% vessels = [];
% while ~feof(fid)
%     
%     %Find next vessel
%     while isempty(strfind(fgetl(fid), 'ncm_vessel:')) && ~feof(fid)
%     end
%     if feof(fid); break; end
%     
%     %We have a vessel, so increment vessel count
%     vessel_num = vessel_num + 1;
%     
%     %Get vessel properties
%     vessels(vessel_num).is_distal = 0;
%     vessels(vessel_num).size = '';
%     vessels(vessel_num).shape = '';
%     
%     %Jump to anchor point
%     while isempty(strfind(fgetl(fid), 'anchor:'))
%     end
%     anchor = textscan(fid, '%f%f');
%     vessels(vessel_num).anchor = [anchor{1} anchor{2}];
%     
%     %Jump to vessel outline points
%     while isempty(strfind(fgetl(fid), 'points:'))
%     end
%     vessel_str = textscan(fid,'%[^}]', 'bufsize', 1e5, 'HeaderLines', 1);
%     if size(vessel_str{1}, 1)
%         vessels(vessel_num).vessel_pts = str2num(vessel_str{1}{1});
%     else
%         vessels(vessel_num).vessel_pts = [];
%     end
%     
%     %Jump to apices
%     apex_num = 0;
%     apex_loop = true;
%     while apex_loop
%         next_line = fgetl(fid);
%         if ~isempty(strfind(next_line, 'inner_point:'))
%             apex_num = apex_num + 1;
%             total_apices = total_apices + 1;
%             apex = textscan(fid, '%f%f');
%             vessels(vessel_num).apex(apex_num).inner = [apex{1} apex{2}];
%             
%         elseif ~isempty(strfind(next_line, 'outer_point:'))
%             apex = textscan(fid, '%f%f');
%             vessels(vessel_num).apex(apex_num).outer = [apex{1} apex{2}];
%             
%         elseif ~isempty(strfind(next_line, '// apices'))
%             apex_loop = false; %break loop
%         end
%     end
% end
% fclose(fid);

[markup] = read_markup_from(vessel_data_file);
vessels = markup.vessels;

num_vessels = length(vessels);

is_distal = false(num_vessels, 1);
apex_data = zeros(num_vessels, 7);
vessel_shapes = cell(num_vessels, 1);
vessel_sizes = cell(num_vessels, 1);

for i_v = 1:num_vessels
    
    is_distal(i_v) = vessels(i_v).ncm_vessel_properties.is_distal;
    vessel_shapes{i_v} = vessels(i_v).ncm_vessel_properties.shape;
    vessel_sizes{i_v} = vessels(i_v).ncm_vessel_properties.size;
    
    if is_distal(i_v) && ~isempty(vessels(i_v).apices) && ...
            isfield(vessels(i_v).apices(1), 'outer_point') && ~isempty(vessels(i_v).apices(1).outer_point);
        
        apex_data(i_v,1:2) = vessels(i_v).apices(1).outer_point;
        apex_data(i_v,3:4) = vessels(i_v).apices(1).inner_point;
    end
end



apex_data(:,5) = (apex_data(:,1) + apex_data(:,3))/2;
apex_data(:,6) = (apex_data(:,2) + apex_data(:,4))/2;
apex_data(:,7) = sqrt(...
    (apex_data(:,1) - apex_data(:,3)).^2 + ...
    (apex_data(:,2) - apex_data(:,4)).^2 );

if nargin > 1
    warning('off', 'MATLAB:xlswrite:AddSheet');
    
    xlswrite(excel_file, {'Distal vessel'}, sheet_name, 'A1');
    xlswrite(excel_file, {'Apex inner'}, sheet_name, 'B1');
    xlswrite(excel_file, {'x'}, sheet_name, 'B2');
    xlswrite(excel_file, {'y'}, sheet_name, 'C2');
    xlswrite(excel_file, {'Apex outer'}, sheet_name, 'D1');
    xlswrite(excel_file, {'x'}, sheet_name, 'D2');
    xlswrite(excel_file, {'y'}, sheet_name, 'E2');
    xlswrite(excel_file, {'Apex centre'}, sheet_name, 'F1');
    xlswrite(excel_file, {'x'}, sheet_name, 'F2');
    xlswrite(excel_file, {'y'}, sheet_name, 'G2');
    xlswrite(excel_file, {'Apex width'}, sheet_name, 'H1');
    xlswrite(excel_file, {'Vessel shape'}, sheet_name, 'I1');
    xlswrite(excel_file, {'Vessel size'}, sheet_name, 'J1');
    
    if num_vessels
        xlswrite(excel_file, is_distal, sheet_name, 'A3');
        xlswrite(excel_file, apex_data, sheet_name, 'B3');
        xlswrite(excel_file, vessel_shapes, sheet_name, 'I3');
        xlswrite(excel_file, vessel_sizes, sheet_name, 'J3');
    end
end
    
    
    
    
    
    

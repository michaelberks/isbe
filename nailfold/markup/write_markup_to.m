function [] = write_markup_to(markup, filename, resize_factor)
% Parse the *_markup.txt file output by the ncm_qmarkup.exe application
% into a structure whose fields have corresponding names.   

if ~exist('resize_factor', 'var')
    resize_factor = 1;
end

% Open file ID
fid = fopen(filename, 'wt');
 
fprintf(fid, '%s \n', 'ncm_annotation: {');

%Write out general details
fprintf(fid, '  %s %d \n', 'version:', markup.version);
fprintf(fid, '  %s %s \n', 'observer:', markup.observer);
fprintf(fid, '  %s %s \n', 'timestamp:', markup.timestamp);
fprintf(fid, '  %s %d \n', 'time_taken:', markup.time_taken);
fprintf(fid, '  %s %s \n', 'image_grade:', markup.image_grade);

%Loop through vessels
num_vessels = length(markup.vessels);

fprintf(fid, '  %s \n', 'vessels: {');
for i_v = 1:num_vessels
    write_vessel(markup.vessels(i_v), fid, resize_factor);
end
fprintf(fid, '  %s \n', '} // vessels');

%Loop through haemorrhages
num_haems = length(markup.haemorrhages);

fprintf(fid, '  %s \n', 'haemorrhages: {');
for i_h = 1:num_haems
    write_haemorrhage(markup.haemorrhages(i_h), fid, resize_factor);
end
fprintf(fid, '  %s \n', '} // haemorrhages');

%Close up opening ncm_annotation brace and close file ID
fprintf(fid, '%s \n', '} // ncm_annotation');
fclose(fid);
%--------------------------------------------------------------------------
%END OF MAIN FUNCTION

%--------------------------------------------------------------------------
function [] = write_vessel(vessel, fid, resize_factor)

fprintf(fid, '    %s \n', 'ncm_vessel: {');
fprintf(fid, '      %s %d \n', 'version: ', vessel.version);

%Write vessel properties
fprintf(fid, '      %s \n', 'ncm_vessel_properties: {');
fprintf(fid, '        %s %d \n', 'version: ', vessel.ncm_vessel_properties.version);
if vessel.ncm_vessel_properties.is_distal
    fprintf(fid, '        %s \n', 'is_distal: yes');
else
    fprintf(fid, '        %s \n', 'is_distal: no');
end
fprintf(fid, '        %s %s \n', 'size: ', vessel.ncm_vessel_properties.size);
fprintf(fid, '        %s %s \n', 'shape: ', vessel.ncm_vessel_properties.shape);
fprintf(fid, '      %s \n', '} // ncm_vessel_properties');

%Write anchor coordinates
fprintf(fid, '      %s \n', 'anchor: {'); 
fprintf(fid, '        %5.3f %5.3f \n', vessel.anchor(1)*resize_factor, vessel.anchor(2)*resize_factor);
fprintf(fid, '      %s \n', '} // anchor');

%write vessel points coordinates
fprintf(fid, '      %s \n', 'points: {');
fprintf(fid, '      %s \n', '// (x, y, width)');
num_points = size(vessel.points, 1);
for i_pt = 1:num_points
    fprintf(fid, '        %5.3f %5.3f %5.3f \n',...
        vessel.points(i_pt,1)*resize_factor, vessel.points(i_pt,2)*resize_factor, vessel.points(i_pt,3));
end
fprintf(fid, '      %s \n', '} // points');

%Write out apices
fprintf(fid, '      %s \n', 'apices: {');
num_apices = length(vessel.apices);
for i_ap = 1:num_apices
    write_apex(vessel.apices(i_ap), fid, resize_factor);
end
fprintf(fid, '      %s \n', '} // apices');

%Close up the vessel brace
fprintf(fid, '    %s \n', '} // ncm_vessel');     

%--------------------------------------------------------------------------
function [] = write_apex(apex, fid, resize_factor)
if isempty(apex.version); return; end

fprintf(fid, '       %s \n', ' ncm_apex: {');
fprintf(fid, '          %s %d \n', 'version: ', apex.version);
fprintf(fid, '          %s \n', 'inner_point: {');
fprintf(fid, '          %5.3f %5.3f \n', apex.inner_point(1)*resize_factor, apex.inner_point(2)*resize_factor);
fprintf(fid, '          %s \n', '} // inner_point');
fprintf(fid, '          %s \n', 'outer_point: {');
fprintf(fid, '          %5.3f %5.3f \n', apex.outer_point(1)*resize_factor, apex.outer_point(2)*resize_factor); 
fprintf(fid, '          %s \n', '} // outer_point');
fprintf(fid, '       %s \n', '} // ncm_apex');

%--------------------------------------------------------------------------
function [] = write_haemorrhage(haemorrhage, fid, resize_factor)

fprintf(fid, '    %s \n', 'ncm_haemorrhage: {');
fprintf(fid, '      %s %d \n', 'version: ', haemorrhage.version);
fprintf(fid, '      %s \n', 'anchor: {'); 
fprintf(fid, '        %5.3f %5.3f \n', haemorrhage.anchor(1)*resize_factor, haemorrhage.anchor(2)*resize_factor);
fprintf(fid, '      %s \n', '} // anchor');
fprintf(fid, '    %s \n', '} // ncm_haemorrhage');



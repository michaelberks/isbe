function [] = write_automeasures_to(markup, filename, resize_factor)
% Parse the *_markup.txt file output by the ncm_qmarkup.exe application
% into a structure whose fields have corresponding names.   

if ~exist('resize_factor', 'var')
    resize_factor = 2;
end

% Open file ID
fid = fopen(filename, 'wt');
 
fprintf(fid, '%s \n', 'ncm_annotation: {');

%Write out general details
if isfield(markup, 'version')
    fprintf(fid, '  %s %d \n', 'version:', markup.version); %No version yet
else
    fprintf(fid, '  %s %d \n', 'version:', 97);
end
fprintf(fid, '  %s %s \n', 'observer:', 'matlab_software');
if isfield(markup, 'timestamp')
    fprintf(fid, '  %s %s \n', 'timestamp:', markup.timestamp);
else
    fprintf(fid, '  %s %s \n', 'timestamp:', datestr(now));
end
fprintf(fid, '  %s %d \n', 'time_taken:', 0);
fprintf(fid, '  %s %s \n', 'image_grade:', 'auto');

%Loop through vessels

fprintf(fid, '  %s \n', 'vessels: {');

% 1) distal vessels 
num_vessels = size(markup.distal.apex_xy,1);
for i_v = 1:num_vessels
    write_vessel(markup.distal, fid, i_v, true, resize_factor);
end

% 1) non-distal vessels 
num_vessels = size(markup.nondistal.apex_xy,1);
for i_v = 1:num_vessels
    write_vessel(markup.nondistal, fid, i_v, false, resize_factor);
end

fprintf(fid, '  %s \n', '} // vessels');

%Loop through haemorrhages
fprintf(fid, '  %s \n', 'haemorrhages: {');
if isfield(markup, 'haemorrhages')
    num_haems = length(markup.haemorrhages);
    for i_h = 1:num_haems
        write_haemorrhage(markup.haemorrhages(i_h), fid);
    end
end
fprintf(fid, '  %s \n', '} // haemorrhages');

%Close up opening ncm_annotation brace and close file ID
fprintf(fid, '%s \n', '} // ncm_annotation');
fclose(fid);
%--------------------------------------------------------------------------
%END OF MAIN FUNCTION

%--------------------------------------------------------------------------
function [] = write_vessel(vessels, fid, i_v, is_distal, resize_factor)

fprintf(fid, '    %s \n', 'ncm_vessel: {');
fprintf(fid, '      %s %d \n', 'version: ', 97);

%Write vessel properties
fprintf(fid, '      %s \n', 'ncm_vessel_properties: {');
fprintf(fid, '        %s %d \n', 'version: ', 97);
if is_distal
    fprintf(fid, '        %s \n', 'is_distal: yes');
else
    fprintf(fid, '        %s \n', 'is_distal: no');
end
   
fprintf(fid, '        %s %s \n', 'size: ', 'auto');
fprintf(fid, '        %s %s \n', 'shape: ', 'auto');
fprintf(fid, '      %s \n', '} // ncm_vessel_properties');

%Write anchor coordinates
fprintf(fid, '      %s \n', 'anchor: {'); 
fprintf(fid, '        %5.3f %5.3f \n', vessels.apex_xy(i_v,1)*resize_factor, vessels.apex_xy(i_v,2)*resize_factor);
fprintf(fid, '      %s \n', '} // anchor');

%write vessel points coordinates
fprintf(fid, '      %s \n', 'points: {');
fprintf(fid, '      %s \n', '// (x, y, width)');

if isfield(vessels, 'apex_aam_fit')
    num_points = size(vessels.apex_aam_fit(:,:,i_v), 1);
    for i_pt = 1:num_points
        fprintf(fid, '        %5.3f %5.3f %5.3f \n',...
            vessels.apex_aam_fit(i_pt,1,i_v)*resize_factor, vessels.apex_aam_fit(i_pt,2,i_v)*resize_factor, 0);
    end
    fprintf(fid, '      %s \n', '} // points');
end

a_theta = angle(vessels.base_orientation(i_v)) / 2;
nx = sin(a_theta);
ny = cos(a_theta);

apex_xy = [...
    vessels.apex_xy(i_v,1) + vessels.width_at_apex(i_v)*[-1 1]*nx/2;
    vessels.apex_xy(i_v,2) + vessels.width_at_apex(i_v)*[-1 1]*ny/2 ]';
apex_xy = apex_xy * resize_factor;

%Write out apices
fprintf(fid, '      %s \n', 'apices: {');
write_apex(apex_xy, fid);
fprintf(fid, '      %s \n', '} // apices');

feature_names = {...
    'width_at_apex', 'apex_width';
    'mean_weighted_width', 'mean_width';
    'min_width', 'min_width';
    'max_width', 'max_width';
    'std_width', 'std_width';
    'base_orientation', 'apex_orientation';
    'apex_aam_score', 'apex_aam_score';
    'candidate_class_probs', 'apex_class_prob';
    'candidate_displacements', 'apex_displacement';};

for i_f = 1:length(feature_names)
    fprintf(fid, '        %s%s %5.3f \n', feature_names{i_f,2}, ': ', vessels.(feature_names{i_f,1})(i_v));
end
ori_entropy = mb_entropy(vessels.orientation_hist(i_v,:));
fprintf(fid, '        %s %5.3f \n', 'orientation_entropy: ', ori_entropy);

%Close up the vessel brace
fprintf(fid, '    %s \n', '} // ncm_vessel');     

%--------------------------------------------------------------------------
function [] = write_apex(apex_xy, fid)

fprintf(fid, '       %s \n', ' ncm_apex: {');
fprintf(fid, '          %s %d \n', 'version: ', 97);
fprintf(fid, '          %s \n', 'inner_point: {');
fprintf(fid, '          %5.3f %5.3f \n', apex_xy(1,1), apex_xy(1,2));
fprintf(fid, '          %s \n', '} // inner_point');
fprintf(fid, '          %s \n', 'outer_point: {');
fprintf(fid, '          %5.3f %5.3f \n', apex_xy(2,1), apex_xy(2,2)); 
fprintf(fid, '          %s \n', '} // outer_point');
fprintf(fid, '       %s \n', '} // ncm_apex');

%--------------------------------------------------------------------------
function [] = write_haemorrhage(haemorrhage, fid)

fprintf(fid, '    %s \n', 'ncm_haemorrhage: {');
fprintf(fid, '      %s %d \n', 'version: ', haemorrhage.version);
fprintf(fid, '      %s \n', 'anchor: {'); 
fprintf(fid, '        %5.3f %5.3f \n', haemorrhage.anchor(1), haemorrhage.anchor(2));
fprintf(fid, '      %s \n', '} // anchor');
fprintf(fid, '    %s \n', '} // ncm_haemorrhage');



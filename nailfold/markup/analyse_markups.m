function analyse_markups()

rootdir = 'U:\projects\nailfold\remote_folder\markup';

d = dir(rootdir);
d = d([d(:).isdir]);

% Use 'all' for all images, 'common' for only those images marked up by
% every rater.
[gradernames_com, imagenames_com] = get_names(rootdir, d, 'common');
[gradernames, imagenames] = get_names(rootdir, d, 'all');

% Use only the nonpilot data
imagenames = setdiff(imagenames, imagenames_com);

n_graders = length(gradernames);
n_images = length(imagenames);
% n_graders = 2; gradernames = gradernames(1:n_graders);

times_taken = cell(n_images, n_graders);
grades_given = cell(n_images, n_graders);
size_count = num2cell(zeros(n_images, n_graders, 4));
shape_count = num2cell(zeros(n_images, n_graders, 4));
num_vessels = cell(n_images, n_graders);
num_paths = cell(n_images, n_graders);
num_images = cell(n_graders, 1);
vessels = {};

size_stats = nan(5, 3, n_graders);

% Cooccurrence of size/shape labels
size_shape_cooccurrence = zeros(4,4);

giant_fid = fopen(fullfile(rootdir, 'giants.log'), 'w');
for grader_index = 1:n_graders
    fprintf(giant_fid, '%s:\r\n', gradernames{grader_index});
    g_count = 0;
    
    filedir = fullfile(rootdir, gradernames{grader_index});
    files = dir(fullfile(filedir,'*c_markup.txt'));
    
    num_images{grader_index} = length(files);
    
    % Go in reverse order to get most recently marked up images first.
    vessel_mat = {};
    vessel_sizes = cell(1,3);
    for ifile = length(files):-1:1
        [s,t] = strtok(files(ifile).name,'#');
        imagename = t(2:end);
        image_index = find(strcmp(imagenames, imagename));
        
        % Ignore files that aren't on the list.
        if isempty(image_index)
            continue;
        end

        % Continue if we've already read a markup for this file.
        if ~isempty(times_taken{image_index, grader_index})
            continue;
        end

        markup = read_markup_from(fullfile(filedir, files(ifile).name));
        
        times_taken{image_index, grader_index} = double(markup.time_taken);
        grades_given{image_index, grader_index} = markup.image_grade;

        % Continue if there are no vessels
        if isempty(markup.vessels)
            num_vessels{image_index, grader_index} = 0;
            continue;
        end
        props = [markup.vessels(:).ncm_vessel_properties];
        
        num_vessels{image_index, grader_index} = numel(markup.vessels);
        
        % Discard any nondistal vessels
        is_distal = find([props(:).is_distal] == 1);
        props = props(is_distal);
                
        giant_count = [0 0 0];

        n_vessels = length(props); 
        v_mat = cell(n_vessels, 6);
        for iv = 1:n_vessels
            vessel = markup.vessels(is_distal(iv));
            
            % Store apex position and estimated vessel size.
            v_mat{iv, 1} = imagename;
            v_mat{iv, 2} = vessel.anchor(1);
            v_mat{iv, 3} = vessel.anchor(2);

            % Update vessel count for each size
            switch (props(iv).size)
                case {'Normal'},    size_index = 1;
                case {'Enlarged'},  size_index = 2;
                case {'Giant'},     size_index = 3;
                case {'Undefined'}, size_index = 4;
            end
            size_count{image_index, grader_index, size_index} = ...
                size_count{image_index, grader_index, size_index} + 1;
            
            % Update vessel count for each shape
            switch (props(iv).shape)
                case {'Normal'},        shape_index = 1;
                case {'Angiogenic', ...
                      'Angiogenetic'},  shape_index = 2;
                case {'Nonspecific', ...
                      'Non-specific', ...
                      'Meandering'},    shape_index = 3;
                case {'Undefined'},     shape_index = 4;
                otherwise, disp(['Unknown shape: ', props(iv).shape]);
            end
            shape_count{image_index, grader_index, shape_index} = ...
                shape_count{image_index, grader_index, shape_index} + 1;

            size_shape_cooccurrence(size_index, shape_index) = ...
                size_shape_cooccurrence(size_index, shape_index) + 1;
                
            if (size_index == 3)
                giant_count(1) = giant_count(1) + 1;
            end
            if (size_index == 2)
                if (shape_index == 2)
                    giant_count(2) = giant_count(2) + 1;
                elseif (shape_index == 3)
                    giant_count(3) = giant_count(3) + 1;
                end
            end
            
            % Store width of first labelled apex (if it exists).
            if ~isempty(vessel.apices)
                apex = vessel.apices(1);
                apex_vec = apex.inner_point - apex.outer_point;
                if ~isempty(apex_vec)
                    width = norm(apex_vec);
                    v_mat{iv, 4} = width;
                    
                    switch (props(iv).size)
                        case {'Normal', 'Enlarged', 'Giant'},
                            vessel_sizes{size_index}(end+1) = width;
                    end
                end
            end            
            
            if ~isempty(vessel.points)
                num_paths{image_index, grader_index} = ...
                    num_paths{image_index, grader_index} + 1;
            end
            
            v_mat{iv, 5} = props(iv).size;
        end % vessel
        
        vessel_mat{image_index} = v_mat;
        
        if (sum(giant_count) > 0)
            g_count = g_count + 1;
            fprintf(giant_fid, '    %10s %2i %2i %2i\r\n', ...
                    imagename, giant_count(1), giant_count(2), giant_count(3));
        end
    end % image
    
    fprintf(giant_fid, '#images = %i\r\n\r\n', g_count);
    
    % Compute mean and standard deviation of apex width for each labelled
    % class.
    for i = 1:3
        vsi = vessel_sizes{i};
        if ~isempty(vsi)
            size_stats(:, i, grader_index) = [length(vsi);
                                              mean(vsi); std(vsi);
                                              min(vsi); max(vsi)];
        else
            size_stats(:, i, grader_index) = [0; NaN; NaN; NaN; NaN];
        end            
    end
    
    vessel_mat = cat(1,vessel_mat{:});
    [n_vessels, n_columns] = size(vessel_mat);
    vessels(1:n_vessels, 1:n_columns, grader_index) = vessel_mat;
    
%     figure(); clf; hold on;
%         title(gradernames{grader_index});
%         x = linspace(0, 200, 101)';
%         px = zeros(length(x), 3);
%         for i = 1:3
%             mn = size_stats(1, i, grader_index);
%             sd = size_stats(2, i, grader_index);
%             xn = x - mn;
%             px(:,i) = (2*pi*sd*sd)^-0.5 * exp(-0.5 * (xn.*xn) / (sd*sd));
%         end
%         plot(x, px);
end % grader
fclose(giant_fid);

save(fullfile(rootdir,'markup_data.mat'),...
     'gradernames', 'imagenames', ...
     'num_images', 'times_taken', 'grades_given', ...
     'num_vessels', 'vessels', ...
     'size_count', 'size_stats', 'shape_count');

size_names = {'Normal';'Enlarged';'Giant';'Undefined'};
shape_names = {'Normal';'Angiogenic';'Nonspecific';'Undefined'};
 
outname = fullfile(rootdir,'markups.xlsx');
n_images = [gradernames(:) num_images(:)];
xlswrite(outname, n_images, 'nImages');

times_taken = [{''} gradernames; 
               imagenames times_taken];
xlswrite(outname, times_taken, 'TimeTaken');

grades_given = [{''} gradernames;
                imagenames grades_given];
xlswrite(outname, grades_given, 'GradeGiven');

num_vessels = [{''} gradernames; 
               imagenames num_vessels];
xlswrite(outname, num_vessels, 'VesselCount');

path_count = [{''} gradernames; 
              imagenames num_paths];
xlswrite(outname, path_count, 'PathCount');

size_count_mat = sum(cell2mat(size_count), 1);
size_count_mat = num2cell(reshape(size_count_mat, [n_graders,4])');
size_count_mat = [{''} gradernames; 
                  size_names size_count_mat];
xlswrite(outname, size_count_mat, 'SizeCount');

shape_count_mat = sum(cell2mat(shape_count), 1);
shape_count_mat = num2cell(reshape(shape_count_mat, [n_graders,4])');
shape_count_mat = [{''} gradernames; 
                   shape_names shape_count_mat];
xlswrite(outname, shape_count_mat, 'ShapeCount');

size_shape_mat = num2cell(size_shape_cooccurrence);
size_shape_mat = [{''} shape_names'; 
                  size_names size_shape_mat];
xlswrite(outname, size_shape_mat, 'SizeShapeCooccurrence');

name_mat = cell(1, n_columns, n_graders);
name_mat(1,1,:) = gradernames;
label_mat = repmat({'Image','Vessel X','Vessel Y','Width','Class',''}, ...
                   [1,1,n_graders]);
vessels = [name_mat; label_mat; vessels];
vessels = vessels(:,:);
xlswrite(outname, vessels, 'Vessels');

header = {'','Normal','','Enlarged','','Giant','';
          '','Mean','StdDev','Mean','StdDev','Mean','StdDev'};
name_mat = cell(1, 2, n_graders);
name_mat(1,1,:) = gradernames;
cols = 2:3;
n_cols = length(cols)*3;
size_stats = reshape(size_stats(cols,:,:), [n_cols, n_graders])';
size_stats = mat2cell(size_stats, ones(1,n_graders), ones(1,n_cols));
size_stats = [header; 
              gradernames(:) size_stats];
xlswrite(outname, size_stats, 'SizeStats');



%% Get the grader names and image names
function [gradernames, imagenames] = get_names(rootdir, d, mode)

if ~exist('mode','var'), mode = 'all'; end

gradernames = {};
imagenames = {};

for id = 1:length(d)
    filedir = fullfile(rootdir, d(id).name);
    files = dir(fullfile(filedir,'*c_markup.txt'));
    if (length(files) < 25) || ...
       (any(strcmp(d(id).name, {'mberks', 'ptresadern', 'gdinsdale'})))
        continue;
    end

    gradernames{end+1} = d(id).name;

    % Get filenames.
    names = cell(length(files),1);
    for ifile = 1:length(files)
        [s,t] = strtok(files(ifile).name,'#');
        names{ifile} = t(2:end);
    end
    
    if isempty(imagenames)
        imagenames = names;
    end
    
    if strcmp(mode, 'all')
        imagenames = union(imagenames, names);
    else
        imagenames = intersect(imagenames, names);
    end
end

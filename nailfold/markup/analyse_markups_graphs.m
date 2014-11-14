clc; 
clear;
close all;

rootdir = 'U:\projects\nailfold\remote_folder\markup';
datafile = fullfile(rootdir, 'markup_data.mat');

outdir = 'S:\projects\mammography\matlab\papers\2013embc (annotation)\fig';

f_export = true;

if ~exist(datafile, 'file')
    return
end

load(datafile);

n_graders = length(gradernames);
n_images = length(imagenames);

size_list = {'Normal','Enlarged','Giant','Undefined'};
n_sizes = length(size_list);

shape_list = {'Normal','Angiogenic','Nonspecific','Undefined'};
n_shapes = length(shape_list);

grade_list = {'Normal','Early','Active','Late','Ungradeable_Extreme',...
              'Ungradeable_Quality','Undefined'};
grade_list_short = {'N','E','A','L','UE','UQ','UD'};
n_grades = length(grade_list);


%% Spearman's rank correlation for the number of vessels in each image

sc = sum(cell2mat(size_count), 3);
rho = [];
for g1 = 1:n_graders
    nv1 = sc(:,g1);
    for g2 = g1+1:n_graders
        nv2 = sc(:,g2);
        
        rho(end+1) = spearmans_corr(nv1, nv2);
    end
end

outstr = sprintf('$\\rho$ = $%04.2f$ (min = $%04.2f$, max = $%04.2f$, n = $%i$)', ...
                 mean(rho), min(rho), max(rho), length(rho));
  
if f_export
    fid = fopen(fullfile(outdir, '../spearmans.tex'), 'w');
    if (fid ~= 0)
        fprintf(fid, '%s', outstr);
        fclose(fid);
    end
end


%% Confusion (or cooccurrence?) matrix for image grades

confusion_mat = zeros(n_grades, n_grades);
for i = 1:n_images
    for g1 = 1:n_graders
        index1 = find(strcmp(grades_given{i, g1}, grade_list));
        
        for g2 = g1+1:n_graders
            index2 = find(strcmp(grades_given{i, g2}, grade_list));
            inds = sort([index1, index2]);
            confusion_mat(inds(1), inds(2)) = ...
                confusion_mat(inds(1), inds(2)) + 1;
        end
    end
end

if f_export
    fid = fopen(fullfile(outdir, '../confusion.tex'), 'w');
    if (fid ~= 0)
        % First row
        for i2 = 1:n_grades
            fprintf(fid, '%8s & ', grade_list_short{i2});
        end
        fprintf(fid, '%8s \\\\ \n', '');
        fprintf(fid, '\\hline\n');
        
        % Subsequent rows
        for i1 = 1:n_grades
            for i2 = 1:n_grades
                if (i1 > i2)
                    fprintf(fid, '%8s & ', '');
                else
                    fprintf(fid, '%8i & ', confusion_mat(i1,i2));
                end
            end
            fprintf(fid, '%8s \\\\ \n', grade_list_short{i1});
        end
        fclose(fid);
    end
end


%% 'Correlation' between abnormal images and abnormal vessels
grade_size = zeros(n_sizes, n_grades);
% grade_shape = zeros(n_shapes, n_grades);
for g = 1:size(vessels,3)
    for v = 1:size(vessels,1)
        if isempty(vessels{v, 1, g})
            break;
        end
        
        image_index = find(strcmp(imagenames, vessels{v, 1, g}));
        image_grade = grades_given{image_index, g};
        grade_index = find(strcmp(grade_list, image_grade));
        
        size_index = find(strcmp(size_list, vessels{v, 5, g}));
        grade_size(size_index, grade_index) = ...
            grade_size(size_index, grade_index) + 1;
        
%         shape_index = find(strcmp(size_list, vessels{v, ???, g}));
%         grade_shape(shape_index, grade_index) = ...
%             grade_shape(shape_index, grade_index) + 1;
    end
end

if f_export
    fid = fopen(fullfile(outdir, '../grade_vs_size.tex'), 'w');
    if (fid ~= 0)
        % First row
        fprintf(fid, '  %12s ', '');
        for i2 = 1:n_grades
            fprintf(fid, '& %8s ', grade_list_short{i2});
        end
        fprintf(fid, '\\\\ \n');
        fprintf(fid, '\\hline\n');
        
        % Subsequent rows
        for i1 = 1:n_sizes
            fprintf(fid, '  %12s ', size_list{i1});
            for i2 = 1:n_grades
                fprintf(fid, '& %8i ', grade_size(i1,i2));
            end
            fprintf(fid, '\\\\ \n');
        end
        fclose(fid);
    end
end


%% Bar chart of image grades
grade_count = zeros(n_graders, n_grades);
valid_images = zeros(n_images, n_graders);
for g = 1:n_graders
    for i = 1:n_images
        grade_index = find(strcmp(grades_given{i, g}, grade_list));
        grade_count(g, grade_index) = grade_count(g, grade_index) + 1;
        
        % If normal, early, active or late then count as valid
        if (1 <= grade_index) && (grade_index <= 4)
            valid_images(i, g) = 1;
        end
    end
end

if (f_export)
    figure(1);
        bar(grade_count, 'stacked');
        xlabel('Rater'); ylabel('Frequency');
        title('Image grade');
        legend(latex2str(grade_list), 'location','bestoutside');
        graph(gcf); set(gcf,'paperposition',[0,0,18,9]);
        exportfig(fullfile(outdir,'graph_image_grade'));
end


%% Bar charts of vessel shape categories

% Bar charts of vessel shapes.
shape_count_mat = cell2mat(shape_count);
shape_count_mat = reshape(sum(shape_count_mat, 1), [n_graders, n_shapes])';
if (f_export)
    figure(1);
        bar(shape_count_mat', 'stacked');
        legend(shape_list, 'location','southeast');
        xlabel('Rater'); ylabel('Frequency');
        title('Vessel shape labelling');
        graph(gcf); exportfig(fullfile(outdir,'graph_shape_frequency'));
end
row_sum = sum(shape_count_mat, 1);
shape_count_mat = shape_count_mat ./ row_sum(ones(1,n_sizes),:);
if (f_export)
    figure(1);
        bar(shape_count_mat', 'stacked');
        legend(shape_list, 'location','southeast');
        xlabel('Rater'); ylabel('Proportion');
        title('Vessel shape labelling');
        ylim([0, 1.05]);
        graph(gcf); exportfig(fullfile(outdir,'graph_shape_proportion'));
end
clear('row_sum');


%% Bar charts of vessel size categories

% Bar charts of vessel sizes.
size_count_mat = cell2mat(size_count);
size_count_mat = reshape(sum(size_count_mat, 1), [n_graders, n_sizes])';
if (f_export)
    figure(1);
        bar(size_count_mat', 'stacked');
        legend(size_list, 'location','southeast');
        xlabel('Rater'); ylabel('Frequency');
        title('Vessel size labelling');
        graph(gcf); exportfig(fullfile(outdir,'graph_size_frequency'));
end
row_sum = sum(size_count_mat, 1);
size_count_mat = size_count_mat ./ row_sum(ones(1,n_sizes),:);
if (f_export)
    figure(1);
        bar(size_count_mat', 'stacked');
        legend(size_list, 'location','southeast');
        xlabel('Rater'); ylabel('Proportion');
        title('Vessel size labelling');
        ylim([0, 1.05]);
        graph(gcf); exportfig(fullfile(outdir,'graph_size_proportion'));
end


%% Parzen windows estimates of width distribution for the three vessel sizes

sizes = size_list(1:3);
all_widths = cell(1,3);
width_range = 0:0.1:150;
prob_vecs = zeros(3, length(width_range));

for g = 1:n_graders
    v = vessels(1:row_sum(g), :, g);
    size_labels = v(:,5);
    width_cells = cell(1,3);
    
    for i = 1:length(sizes)
        inds = find(strcmp(size_labels, sizes{i}));
        widths = cell2mat(v(inds, 4));
        widths(widths==0) = [];
        
        width_cells{1,i} = widths;
        all_widths{1,i} = [all_widths{1,i}; widths];
        
        prob_vecs(i, :) = ncm_parzen(widths, 3, width_range);
    end
end

for i = 1:length(sizes)
    prob_vecs(i, :) = ncm_parzen(all_widths{1,i}, 2, width_range);
end

comparison = sign(prob_vecs(2,:) - prob_vecs(1,:));
crossing_ind = find(diff(comparison) > 0);
crossing_width1 = mean(width_range(crossing_ind(1)+[1,2]));

comparison = sign(prob_vecs(3,:) - prob_vecs(2,:));
crossing_ind = find(diff(comparison) > 0);
crossing_width2 = mean(width_range(crossing_ind(1)+[1,2]));

disp([crossing_width1, crossing_width2]);

if (f_export)
    figure(1); clf; hold on;
        plot(width_range, prob_vecs');
        xlim(width_range([1,end]));
        xlabel('Width (pixels)'); ylabel('Probability');
        title('All raters');
        legend(sizes, 'location','northeast');
        ax = axis;
        plot(crossing_width1([1,1]), 1e6*[-1,1], '--', 'color', 0.75*[1,1,1]);
        h = plot(crossing_width2([1,1]), 1e6*[-1,1], '--', 'color', 0.75*[1,1,1]);
        axis(ax);
    graph(gcf); exportfig(fullfile(outdir,'graph_widths'));
end


%% Parzen windows estimates of width distribution for enlarged vessels,
%  broken down by image grade

all_widths = cell(1,3);
width_range = 0:0.1:150;
prob_vecs = zeros(n_grades, length(width_range));

v_enlarged = [];
for g = 1:n_graders
    v = vessels(1:row_sum(g), :, g);
    inds = strcmp(v(:,5),'Enlarged');

    size_labels = v(:,5);
    width_cells = cell(1,3);

    image_index = find(strcmp(imagenames, vessels{v, 1, g}));
    image_grade = grades_given{image_index, g};
    grade_index = find(strcmp(grade_list, image_grade));

    size_index = find(strcmp(size_list, vessels{v, 5, g}));
    grade_size(size_index, grade_index) = ...
        grade_size(size_index, grade_index) + 1;
    
    
    for i = 1:length(sizes)
        inds = find(strcmp(size_labels, sizes{i}));
        widths = cell2mat(v(inds, 4));
        widths(widths==0) = [];
        
        width_cells{1,i} = widths;
        all_widths{1,i} = [all_widths{1,i}; widths];
        
        prob_vecs(i, :) = ncm_parzen(widths, 3, width_range);
    end
end

for i = 1:length(sizes)
    prob_vecs(i, :) = ncm_parzen(all_widths{1,i}, 2, width_range);
end

comparison = sign(prob_vecs(2,:) - prob_vecs(1,:));
crossing_ind = find(diff(comparison) > 0);
crossing_width1 = mean(width_range(crossing_ind(1)+[1,2]));

comparison = sign(prob_vecs(3,:) - prob_vecs(2,:));
crossing_ind = find(diff(comparison) > 0);
crossing_width2 = mean(width_range(crossing_ind(1)+[1,2]));

disp([crossing_width1, crossing_width2]);

if (f_export)
    figure(1); clf; hold on;
        plot(width_range, prob_vecs');
        xlim(width_range([1,end]));
        xlabel('Width (pixels)'); ylabel('Probability');
        title('All raters');
        legend(sizes, 'location','northeast');
        ax = axis;
        plot(crossing_width1([1,1]), 1e6*[-1,1], '--', 'color', 0.75*[1,1,1]);
        h = plot(crossing_width2([1,1]), 1e6*[-1,1], '--', 'color', 0.75*[1,1,1]);
        axis(ax);
    graph(gcf); exportfig(fullfile(outdir,'graph_widths'));
end


%% Box plots of vessel width for the three defined sizes
sizes = size_list(1:3);
for i = 1:length(sizes)
    width_mat = zeros(0,2);
    
    for g = 1:n_graders
        v = vessels(1:row_sum(g), :, g);
        size_labels = v(:,5);
        inds = find(strcmp(size_labels, sizes{i}));
        widths = cell2mat(v(inds, 4));
        width_mat = [width_mat; 
                     widths g*ones(size(widths))];
    end
    
    boxplot(width_mat(:,1), width_mat(:,2));
    xlabel('Rater'); ylabel('Vessel width (pixels)'); title(sizes{i});
    filename = sprintf('graph_boxplot_%s', sizes{i});
    ylim([-5, 155]);
    graph(gcf); exportfig(fullfile(outdir,filename));
end


function parse_flow_stats

clc;

statsroot = 'U:\tmp\ncm\synthesis\halfsize\flow';

legend_loc = 'southeast';
% ylims = [0, 1.6];
ylims = [];

% errornames = {'dv_mean_abs', 'dv_rms', 'MAE(v)', 'mean(std(v))'};
errornames = {'MAE(v)'};
allkeynames = {'brightness_sigma', '\sigma_{brightness}';
               'contrast_sigma', '\sigma_{contrast}';
               'jitter_sigma', '\sigma_{jitter}';
               'sigma_n', '\sigma_{noise}';
               'flowrate', 'Relative flow rate';
               'n_cells', 'Number of cells';
               'contrast0', 'Baseline contrast'};
% allkeynames = {'flowrate', 'Relative flow rate'};

% allkeynames = {'contrast0'};

for i = 1:size(allkeynames,1)
    keynames = allkeynames(i,:);
        
    statsdir = fullfile(statsroot, keynames{1});
    [keyvalues, errors] = parse_files(statsdir, keynames(:,1), errornames);
    table = sortrows([keyvalues errors], 1);
    % [mae_u, mae_v] = make_table_1d(keyvalues, errors);
    % table = [mae_u nan(size(mae_u,1),1) mae_v];
%     xlswrite(fullfile(statsdir, 'stats.xls'), table);
    figure(1); clf; hold on;
        plot(table(:, 1), table(:, 2:end), '.-', ...
             'linewidth', 1.5, ...
             'markersize', 8);
        legend(latex2str(errornames), 'location', legend_loc);
        xlabel(keynames{2}); 
        ylabel('Error (pixels/frame)');
        if ~isempty(ylims), ylim(ylims); end
        ax = axis; ax([1]) = 0; axis(ax);

    fig_width = 8;
    graph(gcf, fig_width); 
    exportfig(fullfile(statsdir, keynames{1}));
end

return

d = dir(statsroot);
d = d([d(:).isdir]);
for i = 3:length(d)
    filename = fullfile(d(i).name, '*.png')
    copyfile(fullfile(statsroot, filename), statsroot);
end

return

statsdir = fullfile(statsroot, 'c_sigma');
[keyvalues, errors] = parse_files(statsdir, {'contrast_sigma'});
[mae_u, mae_v] = make_table_1d(keyvalues, errors);
table = [mae_u nan(size(mae_u,1),1) mae_v];
xlswrite(fullfile(statsdir, 'stats.xls'), table);
figure(1); clf; hold on;
    plot(mae_u(2:end,1), mae_u(2:end,2), 'r--');
    plot(mae_v(2:end,1), mae_v(2:end,2), 'b-');
    legend({'MAE(u)', 'MAE(v)'},'location', legend_loc);
    xlabel('Contrast standard deviation'); ylabel('Error (pixels/frame)');
    ylim(ylims);
graph(gcf); exportfig(fullfile(statsdir,'contrast_sigma'));

statsdir = fullfile(statsroot, 'j_sigma');
[keyvalues, errors] = parse_files(statsdir, {'jitter_sigma'});
[mae_u, mae_v] = make_table_1d(keyvalues, errors);
table = [mae_u nan(size(mae_u,1),1) mae_v];
xlswrite(fullfile(statsdir, 'stats.xls'), table);
figure(1); clf; hold on;
    plot(mae_u(2:end,1), mae_u(2:end,2), 'r--');
    plot(mae_v(2:end,1), mae_v(2:end,2), 'b-');
    legend({'MAE(u)', 'MAE(v)'},'location', legend_loc);
    xlabel('Translation standard deviation (pixels)'); ylabel('Error (pixels/frame)');
    ylim(ylims);
graph(gcf); exportfig(fullfile(statsdir,'translation_sigma'));

statsdir = fullfile(statsroot, 'sigma_n');
[keyvalues, errors] = parse_files(statsdir, {'sigma_n'});
[mae_u, mae_v] = make_table_1d(keyvalues, errors);
table = [mae_u nan(size(mae_u,1),1) mae_v];
xlswrite(fullfile(statsdir, 'stats.xls'), table);
figure(1); clf; hold on;
    plot(mae_u(2:end,1), mae_u(2:end,2), 'r--');
    plot(mae_v(2:end,1), mae_v(2:end,2), 'b-');
    legend({'MAE(u)', 'MAE(v)'},'location', legend_loc);
    xlabel('Added noise standard deviation (pixels)'); ylabel('Error (pixels/frame)');
    ylim(ylims);
graph(gcf); exportfig(fullfile(statsdir,'sigma_n'));

statsdir = fullfile(statsroot, 'flowrate+n_cells');
[keyvalues, errors] = parse_files(statsdir, {'flowrate', 'n_cells'});
[mae_u, mae_v] = make_table_2d(keyvalues, errors);
table = [mae_u nan(size(mae_u,1),1) mae_v];
xlswrite(fullfile(statsdir, 'stats.xls'), table);
row = 3; col = 6;
figure(1); clf; hold on;
    plot(mae_u(2:end,1), mae_u(2:end,col), 'r--');
    plot(mae_v(2:end,1), mae_v(2:end,col), 'b-');
    legend({'MAE(u)', 'MAE(v)'},'location', legend_loc);
    xlabel('Relative flow rate'); ylabel('Error (pixels/frame)');
    ylim(ylims);
graph(gcf); exportfig(fullfile(statsdir,'flowrate_640cells'));
figure(1); clf; hold on;
    plot(mae_u(1,2:end), mae_u(row,2:end), 'r--');
    plot(mae_v(1,2:end), mae_v(row,2:end), 'b-');
    legend({'MAE(u)', 'MAE(v)'},'location', legend_loc);
    xlabel('#cells'); ylabel('Error (pixels)');
    ylim(ylims);
graph(gcf); exportfig(fullfile(statsdir,'1flow_ncells'));


function [keyvalues, errorvalues] = parse_files(statsdir, keynames, errornames)

d = dir(fullfile(statsdir, '*_stats.txt'));

keyvalues = zeros(length(d), length(keynames));
errorvalues = zeros(length(d), length(errornames));

for i = 1:length(d)
    fid = fopen(fullfile(statsdir, d(i).name));
    
    if (fid ~= 0)
        while ~feof(fid)
            s = fgetl(fid);
            
            % Find the key
            [str,tok] = strtok(s, ':');
            
            key = strtrim(str);
            value = str2double(tok(2:end));
            
            match = find(strcmp(key, keynames));
            if ~isempty(match)
                keyvalues(i, match) = value;
                continue;
            end
            
            % Get the corresponding errors
            [str,tok] = strtok(s, '=');
            
            key = strtrim(str);
            value = str2double(tok(2:end));
            
            match = find(strcmp(key, errornames));
            if ~isempty(match)
                errorvalues(i, match) = value;
                continue;
            end
        end
        
        fclose(fid);
    end
end


function [mae_u, mae_v] = make_table_1d(value_list, errors)

param1list = sort(unique(value_list(:,1)));
n_params = length(param1list);
n_errors = size(errors, 2);

blank = [ NaN,              Na; 
          param1list(:),    zeros(n_params, n_errors) ];
mae_u = blank;
mae_v = blank;

for i = 1:size(value_list, 1)
    param1 = value_list(i, 1);
    i1 = find(param1list == param1) + 1;
    
    mae_u(i1, 2) = errors(i, 1);
    mae_v(i1, 2) = errors(i, 2);
end


function [mae_u, mae_v] = make_table_2d(value_list, errors)

param1list = sort(unique(value_list(:,1)));
param2list = sort(unique(value_list(:,2)));

blank = [ NaN,              param2list(:)'; 
          param1list(:),    zeros(length(param1list), length(param2list)) ];
mae_u = blank;
mae_v = blank;

for i = 1:size(value_list, 1)
    param1 = value_list(i, 1);
    i1 = find(param1list == param1) + 1;
    param2 = value_list(i, 2);
    i2 = find(param2list == param2) + 1;
    
    mae_u(i1, i2) = errors(i, 1);
    mae_v(i1, i2) = errors(i, 2);
end


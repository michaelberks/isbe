data = xlsread('K:\isbe\misc\renal Audit 2010-updated.xlsx', 1, 'i3:j59');
%
means = naNmean(data);
%%
f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [324 269 557 676],...
    'PaperPositionMode','auto');

h = boxplot(data,...
    'labels', {'', ''},...
    'notch', 'on',...
    'positions', [1 2],...
    'widths', 0.8,...
    'whisker', 2); hold on;
for ii = 1:length(h)
    set(h(ii,1), 'linewidth', 2);
    set(h(ii,2), 'linewidth', 2);
end
plot([0.79 1.22], [means(1) means(1)], 'g', 'linewidth', 2);
plot([1.62 2.38], [means(2) means(2)], 'g', 'linewidth', 2);
%title('eGFR data for 2009 and 2010');
set(gca, 'fontsize', 14, 'xtick', [1 2], 'xticklabel', {'2009', '2010'});
ylabel('eGFR (ml/min/1.73m^2)');

print('-dtiff', '-noui', '-painters', f1, '-r300', 'C:\isbe\misc\renal_audit_egr.tif');

%%
display(iqr(data));
display(prctile(data, [5 25 50 75 95]));
display(naNmedian(data));
%%
idx = isnan(data(:,1));
data_pairs = data;
data_pairs(idx,:) = [];
data_pairs = sortrows(data_pairs);

f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [1 29 1280 929],...
    'PaperPositionMode','auto');
hold on;
for ii = 1:size(data_pairs,1);
    plot(ii, data_pairs(ii,1), 'bo', 'markersize', 12);
    plot(ii, data_pairs(ii,2), 'rx', 'markersize', 12);
    if ii == 1
        legend({'2009', '2010'}, 'location', 'southeast');
    end
    plot([ii ii], data_pairs(ii,:), 'k', 'linewidth', 2);
    plot(ii, data_pairs(ii,1), 'bo', 'markersize', 12);
    plot(ii, data_pairs(ii,2), 'rx', 'markersize', 12);
    
end
set(gca, 'fontsize', 18, 'xticklabel', []);
ylabel('eGFR (ml/min/1.73m^2)');
print('-dtiff', '-noui', '-painters', f1, '-r300', 'C:\isbe\misc\renal_audit_egr2.tif');
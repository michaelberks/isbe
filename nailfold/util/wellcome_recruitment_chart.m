study_dirs = [
    dir('N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\0*')
    dir('N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\1*')];

num_subs = length(study_dirs);
recruitment_dates = zeros(num_subs,1);

for i_sub = 1:num_subs; 
    recruitment_dates(i_sub) = study_dirs(i_sub).datenum; 
end;
recruitment_dates(43) = recruitment_dates(44);

end_date = datenum('30-Jul-2015 00:00:00');
start_date = recruitment_dates(1);
curr_date = recruitment_dates(end);
projected_date = start_date + 120*(curr_date-start_date)/num_subs;

if (projected_date < end_date)
    lcolor = 'g';
else
    lcolor = 'r';
end
    
figure('WindowStyle', 'normal'); hold all;
plot(recruitment_dates, 1:num_subs, 'bo', 'markersize', 8, 'markerfacecolor', 'c');
plot(projected_date, 120, [lcolor 'p'], 'markersize', 10, 'markerfacecolor', lcolor);
plot([start_date projected_date], [1 120], [lcolor '--'], 'linewidth', 2);
plot([start_date end_date], [1 120], 'k--', 'linewidth', 2);
legend({...
    'Date subjects imaged',...
    'Projected end date',...
    'Current recruitment rate',...
    'Target rate (30^{th} Jul 2015)'}, ...
    'location', 'southeast');

plot(recruitment_dates, 1:num_subs, 'b-');
plot([projected_date projected_date], [1 120], [lcolor '-.']);


xticks = [...
    datenum('01-Feb-2015')
    datenum('01-Mar-2015')
    datenum('01-Apr-2015')
    datenum('01-May-2015')
    datenum('01-Jun-2015')
    datenum('01-Jul-2015')
    datenum('01-Aug-2015')];
%
xticklabels = {...
    '1st Feb';
    '1st Mar';
    '1st Apr';
    '1st May';
    '1st Jun';
    '1st Jul';
    '1st Aug';};

set(gca, 'xtick', xticks, 'xticklabel', xticklabels, 'FontSize', 12);
title({'Projected vs target recruitment for Wellcome nailfold study';
    ['Recruited: ' num2str(num_subs) ', target: 120']});
ylabel('Number of subjects recruited');
exportfig('N:\Nailfold Capillaroscopy\Wellcome\recruitment_chart.png');
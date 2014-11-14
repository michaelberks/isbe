%Load in the data from XLS spreadsheet
xls_name = 'K:\isbe\students_and_teaching\MPhys\stepwedge_markers_good one.xls';

[~,small_r_names] = xlsread(xls_name, 1, 'A2:A43');
[~,large_r_names] = xlsread(xls_name, 1, 'A45:A61');

[~,small_l_names] = xlsread(xls_name, 1, 'A78:A121');
[~,large_l_names] = xlsread(xls_name, 1, 'A63:A76');

[small_r_data] = xlsread(xls_name, 1, 'B2:Q43');
[large_r_data] = xlsread(xls_name, 1, 'B45:Q61');

[small_l_data] = xlsread(xls_name, 1, 'B78:Q121');
[large_l_data] = xlsread(xls_name, 1, 'B63:Q76');

%%
%Plot each type (small, large, left, right)
f1 = figure; axis equal ij; hold on;
plot(small_r_data(:,1:6), small_r_data(:,9:14), 'rx');

f2 = figure; axis equal ij; hold on;
plot(small_l_data(:,1:6), small_l_data(:,9:14), 'rx');

f3 = figure; axis equal ij; hold on;
plot(large_r_data(:,1:8), large_r_data(:,9:16), 'rx');

f4 = figure; axis equal ij; hold on;
plot(large_l_data(:,1:8), large_l_data(:,9:16), 'rx');
%%
%Separate out the 2 clusters in the small mammos. For rights, it looks like
%thresholding the 3rd marker at 3000 will do it, for lefts, thresholding
%the first marker at 1000px
set1 = small_r_data(:,3) > 3000;
figure(f1); plot(small_r_data(set1,1:6), small_r_data(set1,9:14), 'bx');

set2 = small_l_data(:,1) < 1000;
figure(f2); plot(small_l_data(set2,1:6), small_l_data(set2,9:14), 'bx');

%Set which mammograms these are
small_r_names{set1} %#ok
small_l_names{set2} %#ok

%In each case, it's the 4 views for one case that are in the odd one out
%cluster - suggests the sheet was probably was upside down, or maybe
%back-to-front. You should be able to work out which...






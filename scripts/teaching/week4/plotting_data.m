%--------------------------------------------------------------------------
% ------------------------- Plotting data -----------------------------
%--------------------------------------------------------------------------
pred_dir = 'P:\MatlabTutorial\data\retinograms_roc\retinograms_roc\';

%Load in some data
load([pred_dir 'set1\roc\roc_stats.mat']);

%% We have 2D array roc_pts that gives the coordinates to plot an ROC curve

%1) Create a new figure window
figure; %What happens if we don't do this?

%2) Plot the data using the plot function
plot(roc_pts(:,1), roc_pts(:,2));

%3) Add markers to show each operating point in the ROC curve
plot(roc_pts(:,1), roc_pts(:,2), 'rx');

%% Oops, we've wiped out are line... try again using hold...
figure; 
plot(roc_pts(:,1), roc_pts(:,2));
hold on; %Means future plots are added to the current axes
plot(roc_pts(:,1), roc_pts(:,2), 'rx');

%% 4) Now lets label the axes, and add a title
title('ROC curves for vessel segmentation, DRIVE database');
xlabel('1 - specificity');
ylabel('1 - sensitivity');

%% 5) Force the axes to have the same scale, and limit at 0 and 1
axis equal; 
axis([0 1 0 1]); 

%% ------------------------------------------------------------------------
%Now let add multiple ROC curves
data_names = {'set1', 'set2', 'set3', 'set4'};

figure;
hold all; %What effect does using 'all' instead of 'on' have?
for i_data = 1:4
    load([pred_dir data_names{i_data} '\roc\roc_stats.mat'])
    plot(roc_pts(:,1), roc_pts(:,2), '-', 'linewidth', 2.0);    
end

title('ROC curves for vessel segmentation, DRIVE database');
xlabel('1 - specificity');
ylabel('1 - sensitivity');

axis equal; 
axis([0 1 0 1]);
%% How do we know which line is which? Add a legend...
legend({'Set 1', 'Set 2', 'Set 3', 'Set 4'}, 'location', 'southeast');

%% What if we want to zoom in to just the top left portion?
axis([0 0.5 0.5 1]);




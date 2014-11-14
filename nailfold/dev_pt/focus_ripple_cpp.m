clc; 

filename = 'u:/projects/nailfold/tmp/focus_log.txt';
dat = dlmread(filename);

titles = {'z-pos','sharpness','cyclical', 'weighted-dz', 'speed', 'total'};
ncols = length(titles);

for i = 1:ncols
    subplot(ncols,1,i); 
    plot(dat(:,i),'b-');
    ylabel(titles{i});
end


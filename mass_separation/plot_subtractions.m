function plot_subtractions(mass_files, start, fin, mass_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% author:   Michael Berks
% date:     19/11/2007  15:00
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    mass_path = 'C:\isbe\dev\annotations\';
end
if mass_path(end) ~= filesep
    mass_path(end+1) = filesep;
end

for ii = start:fin
    temp = load([mass_path, mass_files(ii).name]);
    mass = temp.mass; clear temp;
    
    cols = size(mass.mass_ROI, 2);
    y_level = round(mean(mass.mass_outline(:,2)));
    
    prof1 = improfile(mass.mass_ROI,...
        [51 cols-50], [y_level y_level], cols-100);
    prof2 = improfile(mass.mass_spline,...
        [51 cols-50], [y_level y_level], cols-100);
    prof3 = improfile(mass.mass_spline_it,...
        [51 cols-50], [y_level y_level], cols-100);
    prof4 = improfile(double(mass.mass_ROI) - mass.mass_sub_it,...
        [51 cols-50], [y_level y_level], cols-100);
    
    a1 = [51, cols-50, min(prof1), max(prof1)];
    
    figure('Name', ['Mass ', num2str(ii)]);
    subplot(2, 4, 1); 
    imagesc(mass.mass_ROI);
    hold on; axis image; colormap(gray(256));
    plot([1 cols], [y_level y_level], 'y:', 'LineWidth', 1.5);
    
    subplot(2, 4, 2);
    imagesc(mass.mass_spline);
    hold on; axis image; colormap(gray(256));
    
    subplot(2, 4, 3);
    imagesc(mass.mass_spline_it);
    hold on; axis image; colormap(gray(256));
    
    subplot(2, 4, 4);
    imagesc(double(mass.mass_ROI) - mass.mass_sub_it);
    hold on; axis image; colormap(gray(256));
     
    subplot(2, 4, 5); 
    plot(51:cols-50, prof1);
    axis(a1);
    
    subplot(2, 4, 6); 
    plot(51:cols-50, prof2);
    axis(a1);
    
    subplot(2, 4, 7); hold on;
    plot(51:cols-50, prof2, 'r:');
    plot(51:cols-50, prof3);
    axis(a1);
    
    subplot(2, 4, 8); 
    plot(51:cols-50, prof4);
    axis(a1);
end
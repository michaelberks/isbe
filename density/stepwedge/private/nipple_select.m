function [nipple_position] = nipple_select(IMAGE)

% user select the position of the nipple
%disp('click to select the location of the nipple, then press enter');
    
% use impixel to select point
nipple_fig = figure('Name','Select nipple position');
imagesc(max(IMAGE(:))-IMAGE); axis image; colormap(gray(256));
title('Select nipple position');
set(gca, 'xticklabel', [], 'yticklabel', []);

[xi,yi,P] = impixel; %#ok
close(nipple_fig);

if length(xi) ~= 1
    warndlg('Please select 1, and only 1 point', 'Wrong number of points selected', 'modal');
    [nipple_position] = nipple_select(IMAGE);
else
    nipple_position(1) = xi(1);
    nipple_position(2) = yi(1);
end

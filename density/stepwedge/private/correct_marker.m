function [xc, yc] = correct_marker(mammo_patch)
%FIND_MARKER given an initial approximation to a breastt hickness marker
%centre, find the precise centre automatically and estimate the radius of
%the marker
%
% Arguments:
%
% 'mammo'
%   - Mammogram containing breast thickness markers
%
% 'xi,yi'
%   - Initial estimate of marker centre
%
% 'reject'
%   - if 0 (i.e. for first attempt) marker will be found fully
%   automatically. If 1 (i.e. user has marker first attempt as reject),
%   user is asked to select points on edge to detect marker
%
% 'marker_fig'
%   - handle to the figure of marker images (e.g. generated in calling function
%   marker detect) so we can draw ouput in this figure
%
%
% Outputs:
%
% 'xc, yc'
%   - x,y coordinates of marker centre
%
%
% See also: STEPWEDGE, MARKERDETECT
%
% Created: 09-May-2006
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5125 
% Copyright: (C) University of Manchester
    
reject = 1;        
    
%if we have to manually select the points do that now
while reject

    % display the contrast enhanced image
    manual_fig = figure(...
        'Name', 'Manually mark points',...
        'WindowStyle', 'normal',...
        'MenuBar', 'none',...
        'Units', 'normalized',...
        'OuterPosition', [0.25 0.25 0.5 0.5]);

    imagesc(mammo_patch); axis image; colormap(gray(256)); hold all;
    title('Double-click on the marker centre or if no marker is present, press enter');
    set(gca, 'xticklabel', [], 'yticklabel', []);

    %get the user to select points on the edge
    [x_select, y_select, P] = impixel; %#ok
    
    if isempty(x_select)
        answer = questdlg(...
            'No marker is present in this region present','No marker present','Continue', 'Select marker', 'Continue');
        if strcmpi(answer, 'Continue')
            reject = 0;
        end
        xc = [];
        yc = [];
    else
        xc = x_select(end);
        yc = y_select(end);
        plot(xc, yc, 'go', 'MarkerSize', 20);
        plot(xc, yc, 'gx', 'MarkerSize', 20);
        
        answer = questdlg(...
            'Is the selected marker in the correct position?','Marker Point','Yes', 'No', 'Yes');
        if strcmpi(answer, 'yes')
            reject = 0;
        end
    end
end
close(manual_fig);

    
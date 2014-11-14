function [breast_border, breast_air, errorcheck] = accept_breast_border(mammo_fig, mammo, breast_border, breast_air, mammo_name, debug_mode)
%ACCEPT_BREAST_BORDER Interactive user check that breast segmentation is
% acceptable
%   [breast_border,breast_air,errorcheck] = accept_breast_border(mammo,segmentation)
%
% Inputs:
%      mammo - Mammogram
%
%      segmentation - Structure containing outline of breast
%
%
% Outputs:
%      breast_border - Coordinates of breast border
%
%      breast_air - Indices to coordinates lying on the breast air boundary
%      (as opposed to the chest wall/pectoral side)
%
%      errorcheck - Flag if border not acceptable
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Sep-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if ~exist('debug_mode', 'var')
    debug_mode = 0;
end
if ~exist('mammo_name', 'var')
    mammo_name = [];
end

if ~exist('breast_border', 'var') || isempty(breast_border)
   breast_border = [0 0];
end
if ~exist('breast_air', 'var') || isempty(breast_air)
   breast_air = 1;
end

if ~exist('mammo_name', 'var')
    mammo_name = [];
end

if exist('mammo_fig', 'var') && ~isempty(mammo_fig) && ishandle(mammo_fig)
    figure(mammo_fig);

elseif exist('mammo', 'var') && ~isempty(mammo);
    mammo_fig = figure(...
        'Name', ['Breast segmentation: ' mammo_name],...
        'Units', 'normalized',...
        'OuterPosition', [0 0 1 1]);
    imagesc(mammo); axis image; colormap(gray(256)); hold on;
end

plot(breast_border(:,1), breast_border(:,2), 'y', 'linewidth', 2,...
    'XDataSource', 'breast_border(:,1)',...
    'YDataSource', 'breast_border(:,2)');
plot(breast_border(breast_air,1), breast_border(breast_air,2), 'g', 'linewidth', 2,...
   'XDataSource', 'breast_border(breast_air,1)',...
   'YDataSource', 'breast_border(breast_air,2)');

errorcheck = true;
while errorcheck
    answer = questdlg(...
        'Is the breast border acceptable?','Breast border','Yes', 'No', 'Yes');
    if strcmpi(answer, 'yes')
        errorcheck = 0;
    else
        figure(mammo_fig);
        mam_ax = gca;
        title('Starting at the beginning of the breast air boundary, manually segment the breast');
        [~, xi, yi] = roipoly;
        xi(end) = [];
        yi(end) = [];
        d = [0; cumsum(sqrt(diff(xi).^2 + diff(yi).^2))];
        breast_border = interp1(d, [xi yi], linspace(0, d(end), round(d(end)/10)), 'linear');
        breast_border(breast_border<1) = 1;
        breast_border(breast_border(:,1)>size(mammo,2),1) = size(mammo,2);
        breast_border(breast_border(:,2)>size(mammo,1),2) = size(mammo,1);
        
        breast_air = 1; %#ok
        refreshdata(mam_ax, 'caller');
        
        title('Select the start and end of the breast/air boundary');
        points_ok = false;
        while ~points_ok
            figure(mammo_fig);
            [xi,yi,~] = impixel;
            
            if length(xi) ~= 2
                wd = warndlg({['You selected ' num2str(length(xi)) ' point(s) instead of 2 point'];...
                    'Please select again'}, 'Incorrect selection');
                uiwait(wd);
            else
                points_ok = true;
            end
        end
        
        d = (breast_border(:,1)-xi(1)).^2+ (breast_border(:,2)-yi(1)).^2;
        [~, start_idx] = min(d);

        d = (breast_border(:,1)-xi(2)).^2+ (breast_border(:,2)-yi(2)).^2;
        [~, end_idx] = min(d);
        
        if start_idx > end_idx
            breast_air = (start_idx:end_idx)';
        else
            breast_air = (end_idx:-1:start_idx)';
        end
        refreshdata(mam_ax, 'caller');
    end
end

if ~debug_mode
    close(mammo_fig);
end
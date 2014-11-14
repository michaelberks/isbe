function [] = display_radial_hist(radial_bands, g_width)
%DISPLAY_RADIAL_HIST *Insert a one line summary here*
%   [] = display_radial_hist(radial_bands)
%
% Inputs:
%      radial_bands - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 28-Oct-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if nargin < 2
    g_width = 0;
end
[rows cols num_angles] = size(radial_bands);

if g_width
    g_filt = fspecial('gaussian', 6*g_width+1, g_width);
else
    g_filt = 1;
end

figure;
a1 = axes;
i1 = imagesc(imfilter(sum(radial_bands,3), g_filt)); axis image; colormap(jet(256));

set(i1, 'ButtonDownFcn', @display_hist);
figs = inf(4,1);
fig_idx = 1;
max_val = 0.5*max(radial_bands(:));
phi = linspace(0, pi, 100);

    function display_hist(hObject, eventdata) %#ok

        XYZ = get(a1, 'CurrentPoint');
        x = round(XYZ(1,1));
        y = round(XYZ(1,2));

        if (x > 1) && (x <= cols) && (y > 1) && (y <= rows)
            if ishandle(figs(fig_idx))
                figure(figs(fig_idx));
                hold off;
            else
                figs(fig_idx) = figure('windowstyle', 'normal');
            end
            fig_idx = mod(fig_idx,4)+1;
            
            roi = reshape(radial_bands(y-3*g_width:y+3*g_width, x-3*g_width:x+3*g_width,:), [], num_angles);
            rad_vals = g_filt(:)' * roi;
            [theta rho] = polar_bar([rad_vals 0 0], pi*(-1:num_angles)/num_angles);
            
            %Plot radial axes marker
            plot([-max_val max_val], [0 0], 'k--'); hold on;
            plot([0 0], [0 max_val], 'k--');
            for ii = 2:2:10
                plot(ii*max_val*cos(phi)/10, ii*max_val*sin(phi)/10, 'k--');
            end
            
            plot(rho.*cos(theta), rho.*sin(theta)); axis equal;
            axis([-max_val max_val min(rho.*sin(theta)) max_val]);
            set(gca, 'xticklabel', []);
            title(['Orientation histogram for x = ' num2str(x), ' y = ' num2str(y)]);
            xlabel(['Sum = ' num2str(sum(rad_vals),3), ', S.D. = ' num2str(std(rad_vals),3) ', Entropy = ' num2str(mb_entropy(rad_vals),3)]);  
        end

    end
end



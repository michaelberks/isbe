function [error] = subtract_mass_test(test_list,test_path, n1, n2, spacing, sigma, f_method, green)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% author:   Michael Berks
% date:     09/06/2006  11:50
%
% function: Use TPS to interpolate grey-levels from landmark pts in 
% ROI surrounding mass and subtract mass background. Ouput structure of ROIs
% of subtracted mass gray levels
%
%inputs:
%   data_in: structure containing border vector and ROI for each mass
%   n1: no. of pixels from shape border to start of landmark pts
%   n2: no. of pixels from shape border to end of landmark pts
%   spacing: no. of pixels between landmark pts in x and y directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 8
    green = 'biharmTPS'; %set default interpolant as thin-plate
end

% if plot_flag
%     figure
% end

gauss_filt = fspecial('gaussian', 5*sigma, sigma);
N = length(test_list);
error = zeros(N, 2);

for ii = 1:N
    test_data = u_load([test_path, test_list(ii).name]);
    
    display(['Processing mass ', num2str(ii)]);
    
    %extract shape border and ROI from input structure
    shape_vec = test_data.mass_outline;
    shape_ROI = test_data.mass_ROI;
    real_bg = test_data.real_bg;
    
    shape_bw = roipoly(shape_ROI, shape_vec(:,1),...
        shape_vec(:,2));
    
    %display('start iterating');
    under_ssd = 3; iter = 1;
    while under_ssd > 1
        switch f_method
            case 1
                smooth_ROI = imfilter(shape_ROI, gauss_filt, 'symmetric');
                
            case 2
                smooth_ROI = medfilt2(shape_ROI, [sigma sigma], 'symmetric');
                
            case 3
                smooth_ROI = wiener2(shape_ROI, [sigma sigma], 'symmetric');
                
        end

        [bg_estimates p_list] = spline_estimation(smooth_ROI, shape_bw,...
            n1, n2, spacing, green);
        
        %Compute errors
        if iter == 1
            error(ii,1) = sqrt(sum((double(real_bg(:)) - bg_estimates(:)).^2) / length(bg_estimates));
            error(ii,2) = error(ii,1);
            iter = 2;
        end
        %Subtract interpolated values from background and create new ROI of
        %subtracted gray-levels (zero outside of shape)
        break;
        bg_diff = double(shape_ROI(p_list)) - bg_estimates';
        total_ssd = sum(bg_diff.^2) / length(p_list);
        under_ssd = sum(bg_diff(bg_diff > 0).^2) / length(p_list);
        display(['total = ', num2str(total_ssd), ' under = ', num2str(under_ssd)]);
        
        if under_ssd > 0.5*total_ssd   
            shape_ROI(p_list) = uint8(bg_estimates);
            error(ii,2) = sqrt(sum((double(real_bg(:)) - bg_estimates(:)).^2) / length(bg_estimates));
        else
            break;
        end
    end
    %display('end iterating');
    %display('');
    
end

display(['Finished subtracting background:', ...
                ' n1 = ', num2str(n1), ...
                ' n2 = ', num2str(n2), ...
                ' spacing = ', num2str(spacing), ...
                ' sigma = ', num2str(sigma)]);
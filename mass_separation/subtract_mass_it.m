function [data_out] = subtract_mass_it(data_in, path,...
    n1, n2, spacing, sigma, plot_flag, green, f_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% author:   Michael Berks
% date:     11/04/2006  11:50
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

if nargin < 7
    plot_flag = 0; %set default of no plotting
end
if nargin < 8
    green = 'biharmTPS'; %set default interpolant as thin-plate
end
if nargin < 9
    f_method = 3;
end

if f_method == 1;
    gauss_filt = fspecial('gaussian', 5*sigma, sigma);
end

N = length(data_in);

for ii = 1:N
    temp = load([path, data_in(ii).name]);
    mass = temp.mass; clear temp;
    
    %extract shape border and ROI from input structure
    shape_vec = mass.mass_outline;
    shape_ROI = mass.mass_spline;
    
    shape_bw = roipoly(shape_ROI, shape_vec(:,1),...
        shape_vec(:,2));
    
    display('start iterating');
    under_ssd = 3;
    while under_ssd > 1
        switch f_method
            case 1
                smooth_ROI = imfilter(shape_ROI, gauss_filt, 'symmetric');
                f_string = 'Gaussian';
            case 2
                smooth_ROI = medfilt2(shape_ROI, [sigma sigma], 'symmetric');
                f_string = 'Median';
            case 3
                smooth_ROI = wiener2(shape_ROI, [sigma sigma], 'symmetric');
                f_string = 'Wiener';
        end

        [bg_estimates p_list] = spline_estimation(smooth_ROI, shape_bw,...
            n1, n2, spacing, green);
    
        %Subtract interpolated values from background and create new ROI of
        %subtracted gray-levels (zero outside of shape)
        
        bg_diff = double(shape_ROI(p_list)) - bg_estimates';
        total_ssd = sum(bg_diff.^2) / length(p_list);
        under_ssd = sum(bg_diff(bg_diff > 0).^2) / length(p_list);
        display(['total = ', num2str(total_ssd), ' under = ', num2str(under_ssd)]);
        
        if under_ssd > 0.5*total_ssd   
            shape_ROI(p_list) = uint8(bg_estimates);
        else
            break;
        end
    end
    display('end iterating');
    display('');
    
    if plot_flag
        figure('Name', ['Mass', zerostr(ii,3)]);
        imagesc(shape_ROI); axis image; colormap(gray(256));
        hold on;
        %plot(l_c, l_r, 'rx');        
    end
    mass.mass_list = p_list;
    mass.mass_spline_it = shape_ROI;
    clear shape_ROI bg_estimates;
    
    save([path, data_in(ii).name], 'mass');
        
end

data_out = mass;

display(['Finished subtracting background:', ...
                ' n1 = ', num2str(n1), ...
                ' n2 = ', num2str(n2), ...
                ' spacing = ', num2str(spacing), ...
                ' sigma = ', num2str(sigma),...
                ' filter = ', f_string]);
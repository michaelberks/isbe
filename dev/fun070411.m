function [data_out, error] = fun070411(data_in, path,...
    n1, n2, spacing, sigma, plot_flag, green)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make cut and paste mass images
%
if nargin < 3
    n1 = 20; %set defaults of n1, n2, spacing as used by Caulkin
    n2 = 20;
    spacing = 5;
    sigma = 5;
end
if nargin < 7
    plot_flag = 0; %set default of no plotting
end
if nargin < 8
    green = 'biharmTPS'; %set default interpolant as thin-plate
end


gauss_filt = fspecial('gaussian', 5*sigma, sigma);
N = length(data_in);
error = zeros(N, 1);
data_in(N).target_spline = [];
save(path, 'data_in');

for ii = 1:N
    
    %extract shape border and ROI from input structure
    shape_vec = data_in(ii).mass_outline;
    shape_ROI = data_in(ii).target_ROI;    
    clear data_in;
    
    shape_bw = fliplr(roipoly(shape_ROI, shape_vec(:,1),...
        shape_vec(:,2)));
    
    [bg_estimates p_list] = spline_estimation(shape_ROI, shape_bw,...
        n1, n2, spacing, gauss_filt, green);
    
    %Subtract interpolated values from background and create new ROI of
    %subtracted gray-levels (zero outside of shape)
    %subtract_tex = shape_ROI(p_list) - uint8(f_z);
    
    target_spline = shape_ROI;
    target_spline(p_list) = uint8(bg_estimates);
    clear shape_ROI bg_estimates;
    
    if plot_flag
        figure('WindowStyle', 'docked', 'Name', ['Mass', zerostr(ii,3)]);
        imagesc(target_spline); axis image; %colormap gray;
        hold on;      
    end

    load(path);
    data_in(ii).target_spline = target_spline; clear target_spline;
    save(path, 'data_in');    
end
data_out = data_in;

display(['Finished subtracting background:', ...
                ' n1 = ', num2str(n1), ...
                ' n2 = ', num2str(n2), ...
                ' spacing = ', num2str(spacing), ...
                ' sigma = ', num2str(sigma)]);
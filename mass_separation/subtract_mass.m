function [data_out] = subtract_mass(data_in, path,...
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
if nargin < 9
    f_method = 3;
end

if f_method == 1;
    gauss_filt = fspecial('gaussian', 5*sigma, sigma);
end

N = length(data_in);
%error = zeros(N, 1);
%data_in(N).mass_sub = [];
%data_in(N).mass_spline = [];
%save(path, 'data_in');

for ii = 1:N
    temp = load([path, data_in(ii).name]);
    mass = temp.mass; clear temp;
    
    %extract shape border and ROI from input structure
    %shape_vec = data_in(ii).mass_outline;
    %shape_ROI = data_in(ii).mass_ROI;
    shape_vec = mass.mass_outline;
    shape_ROI = mass.mass_ROI;
    
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
    
    shape_bw = roipoly(shape_ROI, shape_vec(:,1),...
        shape_vec(:,2));
    
    [bg_estimates p_list] = spline_estimation(smooth_ROI, shape_bw,...
        n1, n2, spacing, green);
    
    %Subtract interpolated values from background and create new ROI of
    %subtracted gray-levels (zero outside of shape)
    
    mass_spline = shape_ROI;
    mass_spline(p_list) = uint8(bg_estimates);
    
    subtract_tex = double(shape_ROI(p_list)) - bg_estimates';
    subtract_tex2 = double(smooth_ROI(p_list)) - bg_estimates';
    
    mass_sub = zeros(size(shape_ROI));
    mass_sub(p_list) = subtract_tex;
    mass_sub_smooth = zeros(size(shape_ROI));
    mass_sub_smooth(p_list) = subtract_tex2;
    mass.mass_list = p_list;
    mass.mass_sub = mass_sub;
    mass.mass_sub_smooth = mass_sub_smooth;
    mass.mass_spline = mass_spline;
    
    
    save([path, data_in(ii).name], 'mass');
    
    if plot_flag
        %figure('WindowStyle', 'docked', 'Name', ['Mass', zerostr(ii,3)]);
        %imagesc(mass_sub); axis image; %colormap gray;
        %hold on;
        %figure('WindowStyle', 'docked', 'Name', ['Mass', zerostr(ii,3)]);
        %imagesc(mass_spline); axis image; %colormap gray;
        %hold on;
        %plot(l_c, l_r, 'rx');
        temp = smooth_ROI;
        temp(p_list) = bg_estimates;
        figure('WindowStyle', 'docked', 'Name', ['Mass', zerostr(ii,3)]);
        surface(double(temp(1:5:end, 1:5:end))); clear temp;
    end
    clear shape_ROI subtract_tex* bg_estimates mass_sub* mass_spline;
end
%data_out = data_in;
data_out = 0;

display(['Finished subtracting background:', ...
                ' n1 = ', num2str(n1), ...
                ' n2 = ', num2str(n2), ...
                ' spacing = ', num2str(spacing), ...
                ' sigma = ', num2str(sigma),...
                ' filter = ', f_string]);
function [fatty_er streaky_er parenchyma_er dense_er] = ...
    subtract_mass_experiment_mean(data_in, path, ...
    bgs, n1, n2, spacing, sigma, plot_flag)
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

if nargin < 4
    n1 = 20; %set defaults of n1, n2, spacing as used by Caulkin
    n2 = 20;
    spacing = 5;
    sigma = 5;
end
if nargin < 8
    plot_flag = 0; %set default of no plotting
end

gauss_filt = fspecial('gaussian', 5*sigma, sigma);
N = length(data_in);

fatty_er.ssd = zeros(N, 1);
fatty_er.over_pts = zeros(N, 1);
fatty_er.under_pts = zeros(N, 1);
fatty_er.equal_pts = zeros(N, 1);
fatty_er.over_error = zeros(N, 1);
fatty_er.under_error = zeros(N, 1);

streaky_er.ssd = zeros(N, 1);
streaky_er.over_pts = zeros(N, 1);
streaky_er.under_pts = zeros(N, 1);
streaky_er.equal_pts = zeros(N, 1);
streaky_er.over_error = zeros(N, 1);
streaky_er.under_error = zeros(N, 1);

parenchyma_er.ssd = zeros(N, 1);
parenchyma_er.over_pts = zeros(N, 1);
parenchyma_er.under_pts = zeros(N, 1);
parenchyma_er.equal_pts = zeros(N, 1);
parenchyma_er.over_error = zeros(N, 1);
parenchyma_er.under_error = zeros(N, 1);

dense_er.ssd = zeros(N, 1);
dense_er.over_pts = zeros(N, 1);
dense_er.under_pts = zeros(N, 1);
dense_er.equal_pts = zeros(N, 1);
dense_er.over_error = zeros(N, 1);
dense_er.under_error = zeros(N, 1);

for ii = 1:N
    temp = load([path, data_in(ii).name]);
    mass = temp.mass; clear temp;
    
    fatty = bgs.fatty;
    fatty(mass.p_list1500) = fatty(mass.p_list1500) + uint8(mass.mass_sub);
    streaky = bgs.streaky;
    streaky(mass.p_list1500) = streaky(mass.p_list1500) + uint8(mass.mass_sub);
    parenchyma = bgs.parenchyma;
    parenchyma(mass.p_list1500) = parenchyma(mass.p_list1500) + uint8(mass.mass_sub);
    dense = bgs.dense;
    dense(mass.p_list1500) = dense(mass.p_list1500) + uint8(mass.mass_sub);
    
    mass_bw = poly2mask(mass.mass_outline(:,1)-mass.offset1500(1)+750,...
                        mass.mass_outline(:,2)-mass.offset1500(2)+750,...
                        1500, 1500);
    
    %
    % Compute fatty errors
    %%%%%%%%%%%%
    [bg_estimates p_list] = average_estimation(fatty, mass_bw,...
        n1, n2, spacing, gauss_filt);    
    if plot_flag
        fatty_spline = fatty;
        fatty_spline(p_list) = uint8(bg_estimates);
        figure('WindowStyle', 'docked', 'Name', ['Mass', zerostr(ii,3)]);
        image(fatty_spline); axis image; colormap(gray(256));
        hold on;
    end
    fatty_er.ssd(ii) = sqrt(sum((double(bgs.fatty(p_list))...
        - bg_estimates').^2) / length(bg_estimates));
    
    subtract_error = double(bgs.fatty(p_list)) - bg_estimates';
    
    fatty_er.over_pts(ii) = sum(subtract_error < 0);
    fatty_er.under_pts(ii) = sum(subtract_error > 0);
    fatty_er.equal_pts(ii) = sum(~subtract_error);
    fatty_er.over_error(ii) = -sum(subtract_error(subtract_error < 0))/...
        fatty_er.over_pts(ii);
    fatty_er.under_error(ii) = sum(subtract_error(subtract_error > 0))/...
        fatty_er.under_pts(ii);    
    clear bg_estimates subtract_error;
    
    %
    % Compute streaky errors
    %%%%%%%%%%%%
    
    [bg_estimates p_list] = average_estimation(fatty, mass_bw,...
        n1, n2, spacing, gauss_filt);
    if plot_flag
        streaky_spline = streaky;
        streaky_spline(p_list) = uint8(bg_estimates);
        figure('WindowStyle', 'docked', 'Name', ['Mass', zerostr(ii,3)]);
        image(streaky_spline); axis image; colormap(gray(256));
        hold on;
    end
    streaky_er.ssd(ii) = sqrt(sum((double(bgs.streaky(p_list))...
        - bg_estimates').^2) / length(bg_estimates));
    
    subtract_error = double(bgs.streaky(p_list)) - bg_estimates';
    
    streaky_er.over_pts(ii) = sum(subtract_error < 0);
    streaky_er.under_pts(ii) = sum(subtract_error > 0);
    streaky_er.equal_pts(ii) = sum(~subtract_error);
    streaky_er.over_error(ii) = -sum(subtract_error(subtract_error < 0))/...
        streaky_er.over_pts(ii);
    streaky_er.under_error(ii) = sum(subtract_error(subtract_error > 0))/...
        streaky_er.under_pts(ii);    
    clear bg_estimates subtract_error;
    
    %
    % Compute parenchyma errors
    %%%%%%%%%%%%
    
    [bg_estimates p_list] = average_estimation(fatty, mass_bw,...
        n1, n2, spacing, gauss_filt);
    if plot_flag
        parenchyma_spline = parenchyma;
        parenchyma_spline(p_list) = uint8(bg_estimates);
        figure('WindowStyle', 'docked', 'Name', ['Mass', zerostr(ii,3)]);
        image(parenchyma_spline); axis image; colormap(gray(256));
        hold on;
    end
    parenchyma_er.ssd(ii) = sqrt(sum((double(bgs.parenchyma(p_list))...
        - bg_estimates').^2) / length(bg_estimates));
    
    subtract_error = double(bgs.parenchyma(p_list)) - bg_estimates';
    
    parenchyma_er.over_pts(ii) = sum(subtract_error < 0);
    parenchyma_er.under_pts(ii) = sum(subtract_error > 0);
    parenchyma_er.equal_pts(ii) = sum(~subtract_error);
    parenchyma_er.over_error(ii) = -sum(subtract_error(subtract_error < 0))/...
        parenchyma_er.over_pts(ii);
    parenchyma_er.under_error(ii) = sum(subtract_error(subtract_error > 0))/...
        parenchyma_er.under_pts(ii);    
    clear bg_estimates subtract_error;
    
    %
    % Compute dense errors
    %%%%%%%%%%%%
    
    [bg_estimates p_list] = average_estimation(dense, mass_bw,...
        n1, n2, spacing, gauss_filt);
    if plot_flag
        dense_spline = streaky;
        dense_spline(p_list) = uint8(bg_estimates);
        figure('WindowStyle', 'docked', 'Name', ['Mass', zerostr(ii,3)]);
        image(dense_spline); axis image; colormap(gray(256));
        hold on; 
    end
    dense_er.ssd(ii) = sqrt(sum((double(bgs.dense(p_list))...
        - bg_estimates').^2) / length(bg_estimates));
    
    subtract_error = double(bgs.dense(p_list)) - bg_estimates';
    
    dense_er.over_pts(ii) = sum(subtract_error < 0);
    dense_er.under_pts(ii) = sum(subtract_error > 0);
    dense_er.equal_pts(ii) = sum(~subtract_error);
    dense_er.over_error(ii) = -sum(subtract_error(subtract_error < 0))/...
        dense_er.over_pts(ii);
    dense_er.under_error(ii) = sum(subtract_error(subtract_error > 0))/...
        dense_er.under_pts(ii);    
    clear bg_estimates subtract_error;
    
    save([path, data_in(ii).name], 'mass'); clear mass;   
end

display(['Finished subtracting background:', ...
                ' n1 = ', num2str(n1), ...
                ' n2 = ', num2str(n2), ...
                ' spacing = ', num2str(spacing), ...
                ' sigma = ', num2str(sigma)]);
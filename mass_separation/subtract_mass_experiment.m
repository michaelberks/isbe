function [fatty_er streaky_er parenchyma_er dense_er] = subtract_mass_experiment(data_in, path, ...
    bgs, n1, n2, spacing, sigma, plot_flag, green, f_method)
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
if nargin < 9
    green = 'biharmTPS'; %set default interpolant as thin-plate
end
if nargin < 10
    f_method = 1;
end

if f_method == 1;
    gauss_filt = fspecial('gaussian', 5*sigma, sigma);
end
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
    
    if plot_flag
        dense_spline = dense;
        fatty_spline = fatty;
        streaky_spline = streaky;
        parenchyma_spline = parenchyma;
    end
    
    switch f_method
        case 1
            dense = imfilter(dense, gauss_filt);
            fatty = imfilter(fatty, gauss_filt);
            parenchyma = imfilter(parenchyma, gauss_filt);
            streaky = imfilter(streaky, gauss_filt);
        case 2
            dense = medfilt2(dense, [sigma sigma]);
            fatty = medfilt2(fatty, [sigma sigma]);
            parenchyma = medfilt2(parenchyma, [sigma sigma]);
            streaky = medfilt2(streaky, [sigma sigma]);
        case 3
            dense = wiener2(dense, [sigma sigma]);
            fatty = wiener2(fatty, [sigma sigma]);
            parenchyma = wiener2(parenchyma, [sigma sigma]);
            streaky = wiener2(streaky, [sigma sigma]);
        case 4
            dense = alternative_smooth(dense, mass.p_list1500, sigma);
            fatty = alternative_smooth(fatty, mass.p_list1500, sigma);
            parenchyma = alternative_smooth(parenchyma, mass.p_list1500, sigma);
            streaky = alternative_smooth(streaky, mass.p_list1500, sigma);
    end
    %
    % Compute fatty errors
    %%%%%%%%%%%%
    
    [bg_estimates p_list] = spline_estimation(fatty, mass_bw,...
        n1, n2, spacing, green);    
    if plot_flag
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
    
    [bg_estimates p_list] = spline_estimation(streaky, mass_bw,...
        n1, n2, spacing, green);
    if plot_flag   
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
    
    [bg_estimates p_list] = spline_estimation(parenchyma, mass_bw,...
        n1, n2, spacing, green);
    if plot_flag    
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
    
    [bg_estimates p_list] = spline_estimation(dense, mass_bw,...
        n1, n2, spacing, green);
    if plot_flag
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
    
    %save([path, data_in(ii).name], 'mass'); clear mass;   
end

display(['Finished subtracting background:', ...
                ' n1 = ', num2str(n1), ...
                ' n2 = ', num2str(n2), ...
                ' spacing = ', num2str(spacing), ...
                ' sigma = ', num2str(sigma)]);
            
end

function smooth_ROI = alternative_smooth(bg_ROI, px_id, sigma)
    mask = ones(size(bg_ROI));
    mask(px_id) = 0;
    filt = ones(sigma*5)/(25*sigma^2);
    mask = imfilter(mask, filt, 'replicate');
    mask(px_id) = 1;
    gf = fspecial('gaussian', 5*sigma, sigma);
    smooth_ROI = bg_ROI;
    smooth_ROI(px_id) = 0;
    smooth_ROI = double(imfilter(smooth_ROI, gf, 'symmetric')) ./ mask;
    smooth_ROI = uint8(smooth_ROI);
    figure; imagesc(smooth_ROI); axis image; colormap(gray(256));
end
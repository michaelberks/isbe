function [bg_estimates, p_list] = median_estimation...
    (bg_ROI, mass_mask, n1, n2, spacing, filt_par, green)
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

    if nargin < 3
        n1 = 20; %set defaults of n1, n2, spacing as used by Caulkin
        n2 = 20;
        spacing = 10;
        filt_par = 5;
    end
    if nargin < 7
        green = 'biharmTPS'; %set default interpolant as thin-plate
    end
    
    [r c] = size(bg_ROI);
    
    %dilate the shape mask by n1
    mask1 = mass_mask;
    for jj = 1:n1
        mask1 = imdilate(mask1, strel('disk', 1));
    end
    
    %dilate the shape mask by n2
    mask2 = mass_mask;
    for jj = 1:(n1 + n2)
        mask2 = imdilate(mask2, strel('disk', 1));
    end
    
    %subtract mask1 from mask2, result is ring in shape of mass, n1 pixels
    % from the shape border, n2-n1 pixels wide
    mask3 = mask2 - mask1; clear mask2 mask1;
    
    %create a mask of equally spaced dots (as defined by input: spacing) 
    mask4 = zeros([r c]);
    mask4(1:spacing:end, 1:spacing:end) = 1;
    
    %mask of landmark pts is intersection of spaced dots and ring of shape
    mask5 = mask3 & mask4; clear mask3 mask4
    landmark_pts = find(mask5);
    
    fade_mask = mass_mask;
    for jj = 1:n1
        fade_mask = imdilate(fade_mask, strel('disk', 1));
    end
    
    %Define points to be interpolated by TPS - as row vectors
    mass_pl = regionprops(bwlabel(fade_mask, 4), 'PixelIdxList');
    clear mass_mask fade_mask;
    p_list = sort(mass_pl.PixelIdxList);
    [i_y i_x] = ind2sub([r c], p_list);
    i_x = double(i_x)';
    i_y = double(i_y)';
    
    %Define source points for TPS - as row vectors
    [l_r l_c] = ind2sub([r c], landmark_pts);
    s_x = double(l_c)';
    s_y = double(l_r)';
    
    %Define source values to be interpolated by TPS
    smooth_ROI = medfilt2(bg_ROI, [filt_par filt_par]);
    s_z = double(smooth_ROI(landmark_pts))';
    clear smooth_ROI;
    
    try
        %Compute interpolation
        L_inv = tps_weights(s_x, s_y, green);
        f_z = tps_warp(L_inv, s_x, s_y, s_z, i_x, i_y, green);
    catch
    
        f_z(1:length(i_x)) = mean(s_z);
        display('Problem estimating background:');
        err = lasterror;
        display(err.message);
    end
    
    bg_estimates = f_z;            
end
function [data_tex, error] = subtract_mass_bg(data_in, temp_in, temp_out,...
    n1, n2, spacing, sigma, plot_flag, test_flag, green)
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
    test_flag = 0; %set default of no testing
end
if nargin < 10
    green = 'biharmTPS'; %set default interpolant as thin-plate
end

if plot_flag
    figure
end

error = [];
gauss_filt = fspecial('gaussian', 5*sigma, sigma);
N = length(data_in);
save(temp_in, 'data_in'); clear data_in;

for ii = 1:N
    
    %extract shape border and ROI from input structure
    load(temp_in);
    shape_vec = data_in(ii).mass_outline...
        + repmat(data_in(ii).mass_centroid,...
        size(data_in(ii).mass_outline,1), 1);
    
    %size_shape_vec = length(shape_vec);
    
    mass_ROI = data_in(ii).mass_ROI;
    
    [r c] = size(mass_ROI);
    clear data_in;
    
    %compute BW mask and list of belonging to the shape
    shape_bw = roipoly(mass_ROI, shape_vec(:,1), shape_vec(:,2));
    shape_pl = regionprops(bwlabel(shape_bw, 4), 'PixelList');
    tex_vec = shape_pl.PixelList;
    
    %dilate the shape mask by n1
    mask1 = shape_bw;
    for jj = 1:n1
        mask1 = imdilate(mask1, strel('disk', 1));
    end
    
    %dilate the shape mask by n2
    mask2 = shape_bw;
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
    
    %Define source points for TPS - as row vectors
    [l_r l_c] = ind2sub([r c], landmark_pts);
    s_x = double(l_c)';
    s_y = double(l_r)';
    
    %Define points to be interpolated by TPS - as row vectors
    i_x = double(tex_vec(:,1))';
    i_y = double(tex_vec(:,2))';
    
    %Define source values to be interpolated by TPS
    smooth_ROI = imfilter(mass_ROI, gauss_filt);
    s_z = double(smooth_ROI(landmark_pts))';
    clear smooth_ROI;
    
    try
        
        %Compute interpolation
        L_inv = tps_weights(s_x, s_y, green);
        f_z = tps_warp(L_inv, s_x, s_y, s_z, i_x, i_y, green);
    catch
        
        f_z(1:length(i_x)) = mean(s_z);
        err = lasterror;
        display(err.message);
        display(['Problem subtracting background:', ...
                ' shape number ', num2str(ii), ...
                ' n1 = ', num2str(n1), ...
                ' n2 = ', num2str(n2), ...
                ' spacing = ', num2str(spacing), ...
                ' sigma = ', num2str(sigma)]);
    end
    
    if test_flag
        error(ii) = sum((double(mass_ROI(sub2ind([r c], i_y, i_x)))...
            - f_z).^2) / length(f_z);
    end
    
    %Subtract interpolated values from background and create new ROI of
    %subtracted gray-levels (zero outside of shape)
    subtract_tex = ...
        max(mass_ROI(sub2ind([r c], i_y, i_x)) - uint8(f_z), 0);
    clear mass_ROI f_z;
    subtract_ROI = zeros([r c]);
    subtract_ROI(sub2ind([r c], i_y, i_x)) = subtract_tex; 
    clear subtract_tex;
    
    if plot_flag
        subplot(4, 5, mod(ii - 1, 20)+1);
        imagesc(subtract_ROI); colormap gray; axis image;
        hold on;
        plot(l_c, l_r, 'rx');        
    end

    if ii > 1;
        load(temp_out);
    end
    data_tex(ii).subtract_ROI = subtract_ROI; clear subtract_ROI;
    save(temp_out, 'data_tex');
    if ii < N
        clear data_tex;
        pack
    end
    
end

display(['Finished subtracting background:', ...
                ' n1 = ', num2str(n1), ...
                ' n2 = ', num2str(n2), ...
                ' spacing = ', num2str(spacing), ...
                ' sigma = ', num2str(sigma)]);
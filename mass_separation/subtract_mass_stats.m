function [fatty_er streaky_er parenchyma_er dense_er] = subtract_mass_stats(data_in, path, ...
    bgs, n1, n2, spacing, sigma, f_method)
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
%   n2: no. of pixels from start to end of landmark pts
%   spacing: no. of pixels between landmark pts in x and y directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if f_method == 1;
    gauss_filt = fspecial('gaussian', 5*sigma, sigma);
end
N = length(data_in);

fatty_er.diff = zeros(N, 1);
streaky_er.diff = zeros(N, 1);
parenchyma_er.diff = zeros(N, 1);
dense_er.diff = zeros(N, 1);

fatty_er.increase = zeros(N, 1);
streaky_er.increase = zeros(N, 1);
parenchyma_er.increase = zeros(N, 1);
dense_er.increase = zeros(N, 1);

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
    
    switch f_method
        case 1
            dense_smooth = imfilter(dense, gauss_filt);
            fatty_smooth = imfilter(fatty, gauss_filt);
            parenchyma_smooth = imfilter(parenchyma, gauss_filt);
            streaky_smooth = imfilter(streaky, gauss_filt);
        case 2
            dense_smooth = medfilt2(dense, [sigma sigma]);
            fatty_smooth = medfilt2(fatty, [sigma sigma]);
            parenchyma_smooth = medfilt2(parenchyma, [sigma sigma]);
            streaky_smooth = medfilt2(streaky, [sigma sigma]);
        case 3
            dense_smooth = wiener2(dense, [sigma sigma]);
            fatty_smooth = wiener2(fatty, [sigma sigma]);
            parenchyma_smooth = wiener2(parenchyma, [sigma sigma]);
            streaky_smooth = wiener2(streaky, [sigma sigma]);
    end
    
    mass_bw = poly2mask(mass.mass_outline(:,1)-mass.offset1500(1)+750,...
                        mass.mass_outline(:,2)-mass.offset1500(2)+750,...
                        1500, 1500);
                    
    r = 1500; c = 1500;
    
    %dilate the shape mask by n1
    mask1 = mass_bw;
    for jj = 1:n1
        mask1 = imdilate(mask1, strel('disk', 1));
    end
    
    %dilate the shape mask by n2
    mask2 = mask1;
    for jj = 1:n2
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
    landmark_pts = find(mask5); clear mask5;
    
    fatty_er.increase(ii) = (mean(fatty_smooth(landmark_pts)) - ...
        mean(fatty(landmark_pts)));
    streaky_er.increase(ii) = (mean(streaky_smooth(landmark_pts)) - ...
        mean(streaky(landmark_pts)));
    parenchyma_er.increase(ii) = (mean(parenchyma_smooth(landmark_pts)) - ...
        mean(parenchyma(landmark_pts)));
    dense_er.increase(ii) = (mean(dense_smooth(landmark_pts)) - ...
        mean(dense(landmark_pts)));
    
    fatty_er.diff(ii) = sum((fatty_smooth(landmark_pts) - ...
        fatty(landmark_pts)).^2) / length(landmark_pts);
    streaky_er.diff(ii) = sum((streaky_smooth(landmark_pts) - ...
        streaky(landmark_pts)).^2) / length(landmark_pts);
    parenchyma_er.diff(ii) = sum((parenchyma_smooth(landmark_pts) - ...
        parenchyma(landmark_pts)).^2) / length(landmark_pts);
    dense_er.diff(ii) = sum((dense_smooth(landmark_pts) - ...
        dense(landmark_pts)).^2) / length(landmark_pts);
end

display(['Finished subtracting background:', ...
                ' n1 = ', num2str(n1), ...
                ' n2 = ', num2str(n2), ...
                ' spacing = ', num2str(spacing), ...
                ' sigma = ', num2str(sigma)]);
            
end
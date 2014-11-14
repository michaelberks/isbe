function [enhanced_profile] = ...
    enhance_spicules(s_profile, range, width, method)

x = [-range:range]';

if nargin < 4
    method = 'localContrast';
end

switch method
    case 'localContrast'
        %local contrast filtering
        
        localcon_x = zeros(size(s_profile));
        for a = 1:width;

            f_x = exp(-(x.^2) / (2*a^2)) / sqrt(2*pi*a);
            f_x = f_x / sum(f_x);

            if_x = s_profile - imfilter(s_profile, f_x, 'replicate');
            localcon_x = localcon_x + if_x;    
        end
        localcon_x(localcon_x < 0) = 0;
        %figure, imagesc(localcon_x); colormap gray; axis image;
        enhanced_profile = localcon_x;
        
    case 'topHat'
        %top hat filtering - not sure its working!
        
        tophat_x = zeros(size(s_profile));
        figure;
        for a = 1:width/2;
            f_x(abs(x) > a, :) = -1;
            f_x(abs(x) <= a, :) = 1;
            subplot(4, 5, a);
            plot(f_x);

            if_x = imfilter(s_profile, f_x, 'replicate');
            if_x(if_x < 0) = 0;
            tophat_x = tophat_x + if_x;    
        end
        %tophat_x(tophat_x < 0) = 0;
        figure, imagesc(tophat_x); colormap gray; axis image;
        enhanced_profiles.tophat_x = tophat_x;
        
    case 'gaussian2nd'
        %2nd order Gaussian filtering - not sure its working!
        
        gauss2nd_x = zeros(size(s_profile));
        for a = 1:width;

            f_x = ((x.^2 - a^2) / (a^4)) .* exp(-(x.^2) / (2*a^2)) / sqrt(2*pi*a);
            f_x = f_x / sum(f_x);
            subplot(4, 5, a);
            plot(f_x);

            if_x = imfilter(s_profile, f_x, 'replicate');
            %if_x(if_x < 0) = 0;

            gauss2nd_x = gauss2nd_x + if_x;    
        end   
        gauss2nd_x(gauss2nd_x < 0) = 0;
        figure, imagesc(gauss2nd_x); colormap gray; axis image;
        enhanced_profiles.gauss2nd_x = gauss2nd_x;
        
    otherwise
        display([method, ' method not recognised'])
        enhanced_profile = [];
end
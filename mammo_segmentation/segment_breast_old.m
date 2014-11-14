function [ image_out ] = segment_breast_old(image_in)
% segment_breast
% created: 09/05/2006 16:34
% author: Michael Berks
% function: segment the breast border using Ferrari et al (2004)
%           method
i1 = image_in;

[im_r, im_c] = size(i1);
%i2 = uint8(floor(mat2gray(log(1 + double(i1)))*256));
%show contrast enhanced image
%figure, imagesc(i2); axis image; colormap(gray(256));


i3 = lloyd_max(i1); %clear i2;
%show binary image after lloyd-max binarisation
figure;



% disk_size = 10;

if 0
    % morpological opening until only breast region remains
    do = 1; 
    while do
        disk_size = disk_size + 5;
        i4 = imopen(i3, strel('disk', disk_size, 0));
        [ll, no_objects] = bwlabel(i4, 4);
        do = no_objects - 1;
    end
end

% i4 = imopen(i3, strel('disk', disk_size, 0)); clear i3;

%Throw away all but the largest region
[ll, no_objects] = bwlabel(i3, 4);
[nn] = hist(ll(:), no_objects + 1); clear no_objects;
[mm ind] = max(nn(2:end)); clear max;
i4 = not(ll - ind); clear ind ll;
subplot(1,2,1); imagesc(i4); axis image; colormap(gray(256));

%Fill in any holes in this main region
i4 = bwmorph(i4, 'fill');

%We have the largest object, now find the two points (1 in the top half,
% 1 in the bottom) of minimum width
widths = inf(im_r, 1);
for ii = 1:im_r
    w = find(~i4(ii,:), 1);
    if ~isempty(w)
        widths(ii) = w;
    end
end
[x_top, y_top] = min(widths(1:floor(im_r/2)));
[x_bot, y_bot] = min(widths(ceil(im_r/2):end));
y_bot = y_bot + ceil(im_r/2) - 1;

%Use these points as cutoffs for the main breast region
i4([y_top y_bot], :) = 0;

%Now trace the boundary of the main breast, using a start point based on
%the top cutoff
breast_border = bwtraceboundary(i4, [y_top + 1, 1], 'E');

hold on;
plot(x_top, y_top, 'rx', 'MarkerSize', 10)
plot(x_bot, y_bot, 'gx', 'MarkerSize', 10)

%skin_air is in (x,y) co-ordinates not (r, c)
skin_air(:, 1) = breast_border(:, 2);
skin_air(:, 2) = breast_border(:, 1);

%remove any boundary points with repeated y-coordinates
%skin_air = sortrows(skin_air, 2);
%skin_air = skin_air(logical([1; diff(skin_air(:,2)) > 0]),:);

if 0
    %attempt to fit spline to border
    skin_air_y = skin_air(:,2);
    skin_air_x = spline(skin_air(1:20:end, 2), skin_air(1:20:end, 1), skin_air_y);
    clear skin_air;
    skin_air = [skin_air_x, skin_air_y];
    plot(skin_air(:,1), skin_air(:,2),'b');
end

if 0 %improve border using normal profiles
    
    [fx, fy] = gradient(skin_air);

    %normalise fy
    fy = fy ./ [sqrt(sum(fy.^2, 2)), sqrt(sum(fy.^2, 2))];
    
    for ii = 1:length(fy(:,1)) %= number of rows in skin_air
        n1_x = skin_air(ii, 1) - 10*fy(ii, 2);
        n1_y = skin_air(ii, 2) + 10*fy(ii, 1);
        n2_x = skin_air(ii, 1) + 40*fy(ii, 2);
        n2_y = skin_air(ii, 2) - 40*fy(ii, 1);

        if (n1_x >= 1 && n1_x < im_c && n2_x >= 1 && n2_x < im_c && n1_y >= 1 && n1_y < im_r && n2_y >= 1 && n2_y < im_r)
            [cx, cy, cp] = improfile(i1, [n1_x, n2_x], [n1_y, n2_y], 50);
            normal_p(ii, :) = cp'; %#ok
            normal_x(ii, :) = cx'; %#ok
            normal_y(ii, :) = cy'; %#ok

            [hist_n, hist_x] = hist(cp, 16);
            hist_mode = hist_x(find(hist_n == max(hist_n(1:2)), 1) + 1);
            border_idx = find(cp <= hist_mode, 1);
            skin_air(ii,1) = cx(border_idx);
            skin_air(ii,2) = cy(border_idx);
        end

        if 0%not(rem(ii, 10))
            %plot(cx, cy, 'b');
            plot(n1_x, n1_y, 'gx');
            plot(n2_x, n2_y, 'rx');
        end
    end

    %remove any boundary points with repeated y-coordinates
    skin_air = sortrows(skin_air, 2);
    skin_air = skin_air(logical([1; diff(skin_air(:,2)) > 0]),:);
    plot(skin_air(:,1), skin_air(:,2), 'g-');
end

if 0
    %attempt to remove outliers
    
    dd = diff(skin_air);
    ss = std(dd);
    mm = mean(dd);
    test1 = (dd(:,1) < mm(1) + ss(1)) & (dd(:,1) > mm(1) - ss(1));
    test2 = (dd(:,2) < mm(2) + ss(2)) & (dd(:,2) > mm(2) - ss(2));
    test3 = test1 & test2;
    skin_air = skin_air(test3, :);
    
end

%[r, c] = size(i4);
%i5 = roipoly(i4, [1; skin_air(:,1); 1], [1; skin_air(:,2); r]); 
%imshow(i5);

if 0
    %attempt to fit spline to border
    skin_air_y = skin_air(:, 2);
    skin_air_x = spline(skin_air(1:5:end, 2), skin_air(1:5:end, 1), skin_air_y);
    clear skin_air;
    skin_air = [skin_air_x, skin_air_y];
end

if 0
    %attempt to fit polynomial border
    p_co = polyfit(skin_air(:,2), skin_air(:,1), 5);
    skin_air(:,1) = polyval(p_co, skin_air(:,2));
end

if 0
    % attempt fourier smooth
    freq_x = fft(skin_air(:,1), 1024); 
    idx = round(0.1*length(freq_x));
    freq_x(idx:end) = 0;
    ifreq_x = ifft(freq_x, 1024);
    skin_air(:,1) = ifreq_x(1:length(skin_air(:,1)));
end

if 0
    %attempt snake
    plot(skin_air(:,1), skin_air(:,2), 'b-');
    skin_air = active_contour(skin_air(1:10:end,:), i5, 10);
end
%plot(skin_air2(:,1), skin_air2(:,2), 'r-');

%show original image
subplot(1,2,2); imagesc(i1); axis image; colormap(gray(256));
hold on;
plot(skin_air(:,1), skin_air(:,2), 'b-');
image_out = skin_air;

%
%
%
%
%any other rubbish
%
%create list of objects and get object properties
%[ll, no_objects] = bwlabel(i4, 4);
%s  = regionprops(ll, 'PixelList');
%start_pix = s(1).PixelList(1, 1:2);




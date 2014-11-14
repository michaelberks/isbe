function mask_c = breast_edge(IMAGE,sw_lookup,leftright)

%VERSION USED FOR 2002/2003 ANALYSIS BY HENNA!!!

% this program is to try and identify breast periphery areas

% assume sample image is already read in as IMAGE, and has been reduced to 
% managable size.

% find original IMAGE size - so can return final mask at same size
dim_orig = size(IMAGE);


% reduce image a bit further anyway
IMAGE_smaller = imresize(IMAGE,0.5);

% mediam filter the image:
% median filter
IMAGE_smaller=medfilt2(IMAGE_smaller);

figure(8);
imshow(IMAGE_smaller);

% select a region containing only the breast - avoid markers, wedge, etc...
disp('Select ROI containing breast region, but no confusing markers, etc..');
disp('On chest-wall side define breast region tightly');
disp('On nipple side define region loosely - program will autofind here.');
[mask,xi,yi] = roipoly(IMAGE_smaller);
figure;imshow(mask);
text(5,25,'mask chosen');

ymin = min(yi);
ymax = max(yi);

IMAGE_masked = double(mask) .* double(IMAGE_smaller);
% transform mask to be 255 in background and zero in foreground
mask = -255 * (double(mask)-1);

IMAGE_masked = uint8(IMAGE_masked + mask);

figure;imshow(IMAGE_masked);
text(5,25,'image masked');

%leftright = input('is this left or right? (input "l" or "r")','s');
if(leftright == 0)
    % flip right images
    mask = flipdim(mask,2);
    IMAGE_smaller = flipdim(IMAGE_smaller,2);
    IMAGE_masked = flipdim(IMAGE_masked,2);
    %    imshow(IMAGE_smaller);    
end

% first apply filter:
filter = fspecial('gaussian',15,5);
%IMAGE_filt = imfilter(IMAGE_masked,filter);
IMAGE_filt = uint8(round(filter2(filter,IMAGE_masked)));  % imfilter missing from matlab version at christie ???

%figure(2);
%imshow(IMAGE_filt);

dim = size(IMAGE_masked);

mask_b = double(mask>9999);

% fall off will contain the profiles of the fall off region for each line,
% scaled to length 100 pixels
falloff = zeros(dim(1),100);

edge_points=zeros(1,dim(1));
periph_points=zeros(1,dim(1));

periph_val=0;
n_periph=0;
edge_val=0;
n_edge=0;

% first try to define breast area, work line by line:

for j=1:dim(1)

    if(mod(j,100)==0)
        disp(j);
    end
    profile = double(IMAGE_masked(j,:));

    %figure(3);
    %plot(profile);

    edge_point=1;
    start_point=1;
    
    % find where to start looking
    for i=dim(2):-1:1
        if(mask(j,i) < 200)
            start_point = i;
            break;
        end
    end
    
    if(start_point > 7 )
    
      % first look for breast edge
        for i=start_point-6:-1:7
            local_av = mean(profile(i-5:i+5));
            local_std = std(profile(i-5:i+5));
    
            if(profile(i-6) < local_av - 2*local_std - 2)
                edge_point = i-6;
                edge_val = edge_val + profile(i-6);
                n_edge = n_edge + 1;
                break;
            end
        end

        profile = double(IMAGE_filt(j,:));
    
        % now look for change in gradient, indicating periphery
  
        periph_point = 0;
    
        for ii = i-8:-1:6
            % find change in signal over 10 pixels around this point
            %delta = profile(ii-5)-profile(ii+5);
            delta = sum(profile(ii-5:ii))-sum(profile(ii+1:ii+6));
            %disp([ii profile(ii-10) profile(ii+10) delta]);
            if(delta > -1)
                periph_point = ii;
                periph_val = periph_val + profile(ii); % also keep a log of the average pixel val here
                n_periph=n_periph+1;
                break;
            end
        end
        
        mask_b(j,1:max(periph_point,1))=1;
%        mask_b(j,max([periph_point 1]):max([edge_point periph_point 1]))=0.5;
        mask_b(j,max([periph_point edge_point 1]):dim(2))=0;
           
        scalevals = double(IMAGE_filt(j,max([periph_point 1]):max([edge_point periph_point 1])));
        maxi = max(scalevals);
        mini = min(scalevals);
        
        dimscale = size(scalevals);
        if(maxi-mini < 0.5)
            % make sure we don't introduce any NaN's
            maxi=255;
            mini=0;
        end
        scalevals = (scalevals-mini)/(maxi-mini);
            
        if(dimscale(2) > 1.5 & periph_point > 1.5)
            scale_resized = imresize(scalevals,100/dimscale(2));
            % check that we got the right number of bins after the resize
            dimtemp = size(scale_resized);
            if(dimtemp(2)==100)
                %do nothing
            else
                temp = ones(1,100);
                temp(1,1:min(100,dimtemp(2)))=scale_resized(1,:);
                scale_resized = temp;
            end
            falloff(j,:)=1-scale_resized(1,:);
        end
        mask_b(j,max([periph_point 1]):max([edge_point periph_point 1]))=1-scalevals;
        
        figure(8);
        if(leftright==1)
            drawcross(edge_point,j,2,'b');
            if(periph_point > 0.5)
                drawcross(periph_point,j,2,'r');
            end
        else
            drawcross(dim(2)-edge_point,j,2,'b');
            if(periph_point > 0.5)
                drawcross(dim(2)-periph_point,j,2,'r');
            end
        end
        periph_points(j)=periph_point;
        edge_points(j)=edge_point;
    end
end

%figure(4);
%imshow(mask_b);

% try finding the average width of the periphery

diff_points=edge_points-periph_points;
av_diff = round(mean(diff_points(round(ymin)+20:round(ymax)-20)));
mask_c = mask_b > 9999;

% need to include some smoothing of edge point positions
% median smooth to remove spikes
edge_points = medfilt1(edge_points,10);
% mean smooth to remove larger wobbles
edge_points = round(smooth(edge_points,25));

% sum the fall off profiles to find the average shape of the falloff:
av_fall = sum(falloff);
% scale the average fall off to lie between 0 and 1:
maxi = max(av_fall);
mini = min(av_fall);
av_fall = (av_fall-mini)/(maxi-mini);

%figure(11);
%plot(av_fall);
%av_fall_old = av_fall;

periph_val = periph_val/n_periph; % this is now the average grey-level at the inner edge of the margin
edge_val = edge_val/n_edge; % this is now the average grey level at the outer edge
% scale the av_fall to lie between periph_val and edge_val
av_fall = round( (1-av_fall)*(edge_val-periph_val) + periph_val );

av_fall = sw_lookup(av_fall); % av_fall is now in terms of stepwedge thickness
av_fall = av_fall/sw_lookup(round(periph_val)); % scale av_fall so it lies between 1 and 0

av_fall = av_fall';
%av_fall_old

%figure(10);
%plot(av_fall);

for j=1:dim(1)
 %   if(mod(j,100)==0)
 %       disp(j);
 %   end
    numpoints = max([1 edge_points(j)]) - max([1 edge_points(j)-av_diff]) + 1;
    
    % vals is just linear fall off to zero in periphery region
    if(numpoints == 1)
        vals = 1;
    else
        % resize the average fall to be the length of the average falloff distance:
        vals = imresize(av_fall,numpoints/100);
        vals = vals(1,:);
        % imresize is sometimes a bit rubbish - so have to check here that we 
        % have the right number of points
        dimvals = size(vals);
        if(dimvals(2)==numpoints)
            % do nothing
        else
            temp = zeros(1,numpoints);
            temp(1:dimvals(2))=vals;
            vals=temp;
        end
    end
    
 %   if(mod(j,100)==0)
 %       vals
 %       numpoints
 %   end
    
    mask_c(j,1:max([1 edge_points(j)-av_diff]))=1;
    mask_c(j,max([1 edge_points(j)-av_diff]):max([1 edge_points(j)])) = vals;
    mask_c(j,max([1 edge_points(j)]):dim(2))=0;
end

% combine all this fancy stuff with the original mask:
mask = (255-mask)/255;
mask_c = mask_c .* mask;

if(leftright==0)
    % flip images back
    mask_c = flipdim(mask_c,2);
    IMAGE_smaller = flipdim(IMAGE_smaller,2);
    %figure(2);
    %imshow(IMAGE_smaller);
end



figure(5);
imshow(mask_c);

mask_c = imresize(mask_c,dim_orig);


function [mask_b,edgex,edgey] = breast_edge(IMAGE, coarse_edgex, coarse_edgey, leftright, debug_mode)

% this function is to try and identify breast periphery areas

if nargin < 5
    debug_mode = 0;
end

% find original IMAGE size - so can return final mask at same size
dim_orig = size(IMAGE);


% reduce image a bit further anyway
further_resize = 0.5;
IMAGE_smaller = imresize(IMAGE,further_resize);
dim_small = size(IMAGE_smaller);

% median filter
IMAGE_smaller=medfilt2(IMAGE_smaller);

%Build the mask from the user defined coarse edges
mask = poly2mask(coarse_edgex, coarse_edgey, dim_small(1), dim_small(2));

IMAGE_masked = IMAGE_smaller;
IMAGE_masked(~mask) = 4095;

%leftright = input('is this left or right? (input "l" or "r")','s');
if(leftright == 0)
    % flip right images
    mask = flipdim(mask,2);
    IMAGE_masked = flipdim(IMAGE_masked,2);
end

%MB we never seem to use the filtered image so I've commented this out
% first apply filter:
%filter = fspecial('gaussian',15,5);
%IMAGE_filt = uint8(round(filter2(filter,IMAGE_masked)));  % imfilter missing from matlab version at christie ???
%IMAGE_filt = imfilter(IMAGE_masked, filter);

dim = size(IMAGE_masked);
edgex = zeros(1,dim(1));

% first try to define breast area, work line by line in y (horizontally):
for j = 1:dim(1)

    if debug_mode && (mod(j,100)==0)
        disp(j);
    end
    profile = double(IMAGE_masked(j,:));

    edge_point = 1;    
    start_point = find(mask(j,:), 1, 'last');
    
    if(start_point > 7 )
    
      % first look for breast edge
        for i=start_point-6:-1:7
            local_av = mean(profile(i-5:i+5));
            local_std = std(profile(i-5:i+5));
    
            if(profile(i-6) < local_av - 2*local_std - 2)
                edge_point = i-6;
                break;
            end
        end
        edgex(j)=edge_point;
    end
end

if debug_mode
    figure('Name', 'Reduced size mammogram with unsmoothed edge points fitted');
    imagesc(-double(IMAGE_masked)); axis image; colormap(gray(256)); hold on;
    plot(edgex, 1:dim(1), 'r.');
end

% need to include some smoothing of edge point positions
% median smooth to remove spikes
edgex = medfilt1(edgex,10);
% mean smooth to remove larger wobbles
edgex = round(smooth(edgex,25));

%Scale points back to size of input image
edgex = edgex / further_resize;
edgey = (1:dim(1)) / further_resize;

% Need to take intersection of user generated mask and automatic mask -
% user mask must be elarged first
% automatic mask needs points [1,1] and [1, dim_orig(1)] at start/end
mask_b = poly2mask([1 edgex 1], [1 edgey dim_orig(1)], dim_orig(1), dim_orig(2));
mask = imresize(mask, dim_orig);
mask_b = mask_b & mask;

%Throw away any edge points on the chest wall at start/end of boundary
not_chestwall = find(edgex > 1);
edgex([1:not_chestwall(1), not_chestwall(end):end]) = [];
edgey([1:not_chestwall(1), not_chestwall(end):end]) = [];

% flip images and edge points back if processing right hand breast image
if(leftright==0)
    % flip images back
    mask_b = flipdim(mask_b,2);
    % flip edge points too
    edgex = edgex*-1. + dim_orig(2) + 1;
end

%Plot the fitted edge on the original image
if debug_mode
    figure('Name', 'Original mammogram with edge fitted');
    imagesc(-double(IMAGE)); axis image; colormap(gray(256)); hold on;
    plot(edgex, edgey, 'g', 'LineWidth', 1.5);
end

%--------------------------------------------------------------------------
% MB:
% I don't get any of the following steps. We've found what we think is the
% boundary pixel on each row of the mammogram, giving us a set of (x,y)
% coordinates as the boundary. All we need to do now is enlarge these
% co-ordinates by the rescale factor. As I have done above.
% Instead we do the following:
% 1) Create a binary mask from the boundary co-ordinates
% 2) Resize the binary mask
% 3) Locate the perimeter of this mask (which does some strange things!)
% 4) Convert this back to an edge
%
% As a result I've commented out the code below, and will delete upon
% approval later
%
%--------------------------------------------------------------------------
% % % mask_b will be new mask taking edge finding into account
% % mask_b = zeros(dim);
% % 
% % for j=1:dim(1)
% %     
% %     %MB change here to get rid of edge points points from chest wall
% %     %boundary
% %     %mask_b(j,1:max([1 edge_points(j)]))=1;
% %     mask_b(j,1:edge_points(j))=1;
% % end
% % 
% % %resize mask back to size of image passed in, without applying roi
% % %correction (to find edge)
% % mask_b_noroi = imresize(mask_b, dim_orig);
% % 
% % %work out edge points on resized mask by using bwperim and sparse
% % breastperim = bwperim(mask_b_noroi);
% % 
% % %zero points on edge of image (not part of breast edge)
% % breastperim(:,1)=0;
% % 
% % %MB don't bother with these images
% % %figure(517);
% % %imshow(breastperim,[]);
% % 
% % [edgey, edgex] = find(sparse(breastperim));
% % %numedgepts = size(edgey);
% % 
% % % need to sort edge points so they are in x order, so edge_profile can work
% % % along the edge
% % edgematrix = [edgey edgex];
% % edgematrix = sortrows(edgematrix,1);
% % edgey = edgematrix(:,1);
% % edgex = edgematrix(:,2);
% % 
% % %save these to data file in case we need them again
% % save(filename, 'xi', 'yi', 'edgey', 'edgex');
% % 
% % % combine all the new mask with the original mask:
% % %mask = (255-mask)/255;
%--------------------------------------------------------------------------





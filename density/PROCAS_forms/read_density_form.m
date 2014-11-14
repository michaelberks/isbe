function [density_percent nhs_id, dob, debug_out] = ...
    read_density_form(image_in, debug_flag)
%READ_DENSITY_FORM Given a scanned image of a density form, read in the 4
% density percentages from the analogue scale and the patient ID number
%   [density_percent debug_out] = read_density_form(image_in, debug_flag)
%
% Inputs:
%      image_in- scanned image of density form
%
%      debug_flag- flag {0,1} to turn debug_flag output on or off
%
%
% Outputs:
%      patient_id_str - Patient ID read from form and coverted to a string
%
%      density_percent - Percentages read from each of the 4 line-scales
%
%      debug_out- Structure of information that may be useful to debug the
%      method in case of warnings or errors
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Feb-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% Set the default for the debug parameter to 0
if nargin < 2
    debug_flag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the co-ordinates the image markers and the line-scales in the
% original true form

%image_markers
markers_orig = [25 25; 25 825; 575 825];

%start and end points for the 4 line-scales
line_scale_orig = [ 101 313; 500 313;...
                    101 438; 500 438;...
                    101 563; 500 563;...
                    101 688; 500 688 ];

%Start and end points for the barcode
barcode_orig = [340 0; 585 0; 340 19; 585 19];

[box boy] = meshgrid(340:585, 1:20);
%[lox(:,:, loy] = 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make sure image_in is a double scaled from [0 1] and of optimum size
image_in = double(image_in);
image_in = (image_in - min(image_in(:))) / (max(image_in(:)) - min(image_in(:)));

rescale = 840/size(image_in, 2);
%width should 800 so
image_small = imresize(image_in, rescale, 'bilinear');

%Take binary map of images at 0.7 threshold
image_small = image_small > 0.7;
image_in = image_in > 0.7;

%load in markers (marker1, marker2, marker3)
load components.mat marker*

%Apply normalised cross-correlation to find the markers in the image
C1 = normxcorr2(marker1, image_small);
C2 = normxcorr2(marker2, image_small);

[debug_out.C1_xcorr, ind1] = max(C1(:));
[debug_out.C2_xcorr, ind2] = max(C2(:));

%if xcorr values are below some threshold (should empirically test what a
%suitable value is) - warn the user we may not have correctly found the
%image markers
debug_out.warning = 0;
if debug_out.C1_xcorr < 0.6 || debug_out.C2_xcorr < 0.6
    display('Warning: low cross-correlation scores. We may not have correctly found the image markers');
    debug_out.warning = 1;
end

%Find the marker positions and correct for the size of the marker template
[marker_pos1(2) marker_pos1(1)] = ind2sub(size(C1), ind1); %work in x,y space
marker_pos1 = marker_pos1 - [25 25];

[marker_pos2(2) marker_pos2(1)] = ind2sub(size(C1), ind2); %work in x,y space
marker_pos2 = marker_pos2 - [25 25];


%check the image is the correct way up, (i.e. the y component of the second marker
% is greater than the y component of the first). If not, rotate the image through 180
if  marker_pos2(2) < marker_pos1(2)
    
    %rotate image
    image_small = rot90(image_small, 2);
    
    %correct marker_pos1 and marker_pos2 for rotated image;
    [r c] = size(image_small);

    marker_pos1(1) = c + 1 - marker_pos1(1);
    marker_pos1(2) = r + 1 - marker_pos1(2);
    marker_pos2(1) = c + 1 - marker_pos2(1);
    marker_pos2(2) = r + 1 - marker_pos2(2);

end

%Now the image is the correct orientation, search for the 3rd marker and
%4th marker, depending on which we find, we've either got an old or a new
%form (note we don't need to search through the whole image, nor establish
%the marker location)
C3 = normxcorr2(marker3, image_small);
[debug_out.C3_xcorr ind3] = max(C3(:));

C4 = normxcorr2(marker4, image_small);
[debug_out.C4_xcorr ind4] = max(C4(:));

%New form if we've found a closer match to the 4th marker than the 3rd
new_form = debug_out.C3_xcorr < debug_out.C4_xcorr;

if new_form
    [marker_pos3(2) marker_pos3(1)] = ind2sub(size(C4), ind4); %work in x,y space
else
    [marker_pos3(2) marker_pos3(1)] = ind2sub(size(C3), ind3); %work in x,y space
end
marker_pos3 = marker_pos3 - [25 25];

%Compute the rigid transform (scaling, translation and rotation) of the
%original marker coordinates onto the marker positions found in this image
markers_image = [marker_pos1; marker_pos2; marker_pos3];
T = align_points_affine(markers_image', markers_orig');

%[d,Z,transform] = mb_procrustes(markers_image, markers_orig);

%If we're debugging plot the markers on the form
if debug_flag
    im_fig = figure; imagesc(image_small); axis image; colormap(gray(256)); hold on;
    plot(marker_pos1(1), marker_pos1(2), 'rx', 'MarkerSize', 10);
    plot(marker_pos2(1), marker_pos2(2), 'bx', 'MarkerSize', 10);
    plot(marker_pos3(1), marker_pos3(2), 'gx', 'MarkerSize', 10);
    
    plot(markers_image(1,1), markers_image(1,2), 'ro', 'MarkerSize', 20);
    plot(markers_image(2,1), markers_image(2,2), 'bo', 'MarkerSize', 20);
    plot(markers_image(3,1), markers_image(3,2), 'go', 'MarkerSize', 20);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the barcode - take an image profile at 10 locations through the
% bar_code
% Use the transform to calculate start and end points for the barcode
barcode_affine = transform_points_affine(barcode_orig', T)';
bo_affine = transform_points_affine([box(:)'; boy(:)'], T);

box_affine = reshape(bo_affine(1,:), size(box));
boy_affine = reshape(bo_affine(2,:), size(boy));

%barcode_affine = geom_transformpoints(barcode_orig', T)';
if debug_flag
    plot(...
        [barcode_affine(1,1) barcode_affine(2,1)...
         barcode_affine(3,1) barcode_affine(4,1) barcode_affine(1,1)],...
        [barcode_affine(1,2) barcode_affine(2,2)...
         barcode_affine(3,2) barcode_affine(4,2) barcode_affine(1,2)],...
        'g--');
    plot(box_affine, boy_affine, 'b.');
end

%Compute the barcode points with respect to the large image
box_affine_large = box_affine / rescale;
boy_affine_large = boy_affine / rescale;

%Now take a small region from the image to speed up interpolation
sc = floor(min(box_affine_large(:))) - 1;
ec = ceil(max(box_affine_large(:))) + 1;
sr = floor(min(boy_affine_large(:))) - 1;
er = ceil(max(boy_affine_large(:))) + 1;

barcode_region = image_in(sr:er, sc:ec);
box_affine_large = box_affine_large - sc + 1;
boy_affine_large = boy_affine_large - sr + 1;

barcode_region2 = interp2(barcode_region, box_affine_large, boy_affine_large, 'nearest');

%Transform barcode region into bar code
barcode_sum = sum(barcode_region2);
barcode_sum = barcode_sum > (size(barcode_region2,1) / 2);
start_c = find(~barcode_sum, 1);
end_c = find(~barcode_sum, 1, 'last');
barcode = barcode_sum(start_c:end_c);
%barcode_bin = imresize(barcode, [1 56]);
[barcode_bin debug_out.barcode_warning] = convert_barcode(barcode, 56);

if debug_flag == 2
    
    figure; 
    subplot(5,1,1); imagesc(barcode_region); axis image; colormap(gray(256)); hold on;
    plot(box_affine_large, boy_affine_large, 'b.');
    subplot(5,1,2); imagesc(barcode_region2); axis image; colormap(gray(256));
    subplot(5,1,3); imagesc(barcode_sum); colormap(gray(256));
    subplot(5,1,4); imagesc(barcode); colormap(gray(256));
    subplot(5,1,5); imagesc(barcode_bin); colormap(gray(256));
end
clear image_in;

nhs_bin = barcode_bin(2:35);
dd_bin = barcode_bin(36:40); 
mm_bin = barcode_bin(41:44); 
yy_bin = barcode_bin(45:55);

%Convert the binary numbers back to decimal numbers
nhs_id = 0;
for ii = 33:-1:0
    nhs_id = nhs_id + nhs_bin(34-ii)*2^ii;
end

dd = 0;
for ii = 4:-1:0
    dd = dd + dd_bin(5-ii)*2^ii;
end

mm = 0;
for ii = 3:-1:0
    mm = mm + mm_bin(4-ii)*2^ii;
end

yy = 0;
for ii = 10:-1:0
    yy = yy + yy_bin(11-ii)*2^ii;
end

%Finally group together dd, mm, yy to build the DOB string
dob = [zerostr(dd,2), '/', zerostr(mm,2), '/', num2str(yy)];

if debug_flag
    display(['NHS id: ', num2str(nhs_id), ' DOB: ', dob]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the transform to calculate start and end points for line profiles in
% the image
line_scales_image = zeros(8, 2, 8);
for ii = 1:4
    offset = 2*(ii+1);
    
    %Compute the start and end points in the original co-ordinate frame
    ls = line_scale_orig + repmat([0 offset], 8, 1);
    
    %Transform to the frame in this image
%     line_scales_image(:,:,ii) = ...
%         transform.b * ls * transform.T + repmat(transform.c(1,:), 8, 1);
    
    line_scales_image(:,:,ii) = transform_points_affine(ls', T)';
    
    %Repeat for the negative offset
    ls = line_scale_orig - repmat([0 offset], 8, 1);
%     line_scales_image(:,:,ii+4) = ...
%        transform.b * ls * transform.T + repmat(transform.c(1,:), 8, 1);

    line_scales_image(:,:,ii+4) = transform_points_affine(ls', T)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For each line, take line profiles at offsets parallel to the line-scale,
% these will be 1 where this is no mark and zero otherwise. To overcome
% possible speckle noise, we take the sum of the line profiles. Each offset
% then effectively votes whether a mark is present. If a mark is found at
% all 8 offsets the profile sum will be zero at this point. In general we
% are looking for the minimum point along the sum. Where several points are
% joint minima (this is likely if the mark is thick), we find the middle
% position of the mark.

%Pre-allocate space for the sum of the line profiles
profile_sum = zeros(400,4);
profile_min = zeros(4,1);
density_percent = zeros(4,1);


if debug_flag == 2
    prof_fig = figure;
end

debug_out.error = 0;
debug_out.no_mark = zeros(4,1);
for ll = 1:4   
    
    idx = [2*ll-1 2*ll];
    for ii = 1:8     
        
        %Take the line profile
        l_profile = improfile(image_small, line_scales_image(idx,1,ii), line_scales_image(idx,2,ii), 400, 'bilinear');
        
        %Check no points on the profile are NaN - this would imply a start
        %or end point is outside the range of the image, and would mean our
        %transforms are incorrect. In this case set the debug_out.error to
        %1
        if any(isnan(l_profile))
            display('Error: transformed start/end point out of range');
            debug_out.error = 1;
        else
        %If no NaNs, add the profile to the profile sum
            profile_sum(:,ll) = profile_sum(:,ll) + l_profile;
        end
        
        if debug_flag == 2
            figure(im_fig);
            plot(line_scales_image(idx,1,ii), line_scales_image(idx,2,ii), 'r:');
            figure(prof_fig); subplot(4, 1, ll); hold all;
            plot(1:400, l_profile);
        end
        
    end
    %Find the minimum in the profile sum
    profile_min(ll) = min(profile_sum(:,ll));
    
    if profile_min(ll) >= 7
        %This suggests no marker was found on the line - save the
        %percentage as zero and flag a warning to the user
        debug_out.no_mark(ll) = 1;
        density_percent(ll) = 0;
        display(['Warning: no marker was found on line ' num2str(ll),...
            '. Density percentage has been set to zero'...
            '. Run with debug_flag set to 1 to check the form image']);
    else  
        %Find the middle of the points that correspond to the minimum (note by
        %taking the median we further protect ourselves against a 'stray'
        %mark on the line) - divide by 4 to obtain the percentage
        density_percent(ll) = round(median(find(profile_sum(:,ll) == profile_min(ll))) / 4);
    end
    
    if debug_flag
        figure(im_fig);
        x_pos = line_scales_image(idx(1),1,8)...
            + density_percent(ll)*diff(line_scales_image(idx,1,8))/100;
        y_pos = max(line_scales_image(idx,2,8));
        text(x_pos, y_pos, [num2str(density_percent(ll)), '%']);
    end
end

%Save the profile sums for debugging
debug_out.profile_sum = profile_sum;

%On new forms we have swapped the order of lines so that they are now RML,
%LML, RCC, LCC. This perm required to go from the old lines to the new is:
% {1,2,3,4}->{3 1 4 2}
if new_form
    debug_out.profile_sum = debug_out.profile_sum(:,[3 1 4 2]);
    debug_out.no_mark = debug_out.no_mark([3 1 4 2]);
    density_percent = density_percent([3 1 4 2]);
end
%
end

%--------------------------------------------------------------------------
function transform = align_points_affine(target, source)

    M = size(target,2);
    X = [target' ones(M,1)];

    % just solve for the first two columns of T
    U = source';

    % We know that X * T = U
    if rank(X)>=3
        [Q,R]= qr(X,0);
        Tinv = R\(Q'*U);
    else
        error('GEOM:RankError:Affine','Affine Transform: Source points matrix has rank < 3');
    end

    % add third column
    Tinv(:,3) = [0 0 1]';

    transform = inv(Tinv);
    transform = transform(:,1:2);
    transform = transform(:);
end

%--------------------------------------------------------------------------
function points = transform_points_affine(pointsin, T)
    
    X = [pointsin' ones(size(pointsin,2),1)];
    points= (X * reshape(T,[3 2]))';

end


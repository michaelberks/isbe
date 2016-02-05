function [accum, axis_rho, axis_theta, varargout] = ...
    Hough_Grd_ANEWBEG(img, varargin)
%Detect lines in a grayscale image
%
%  [accum, axis_rho, axis_theta, lineprm, lineseg, dbg_label] =
%      Hough_Grd(img, grdthres, detsens)
%  Hough transform for line detection based on image's gradient field.
%  NOTE:    Operates on grayscale images, NOT B/W bitmaps.
%           NO loops involved in the implementation of Hough transform,
%               which makes the operation fast.
%           Able to detect the two ends of line segments.
%
%  INPUT: (img, grdthres, detsens)
%  img:         A 2-D grayscale image (NO B/W bitmap)
%  grdthres:    (Optional, default is 8, must be non-negative)
%               The algorithm is based on the gradient field of the
%               input image. A thresholding on the gradient magnitude
%               is performed before the voting process of the Hough
%               transform to remove the 'uniform intensity' (sort-of)
%               image background from the voting process. In other words,
%               pixels with gradient magnitudes smaller than 'grdthres'
%               are considered not belong to any line.
%  detsens:     (Optional, default is 0.08)
%               A value that controls the sensitivity of line detection.
%               The range of the value is from 0 to 1, open ends.
%               Although the physical meaning of this parameter is not
%               obvious, it controls the algorithm in a fuzzy manner.
%               The SMALLER the value is, the MORE features in the image
%               will be considered as lines.
%
%  OUTPUT: [accum, axis_rho, axis_theta, lineprm, lineseg, dbg_label]
%  accum:       The result accumulation array from the Hough transform.
%               The horizontal dimension is 'theta' and the vertical
%               dimension is 'rho'.
%  axis_rho:    Vector that contains the 'rho' values for the rows in
%               the accumulation array.
%  axis_theta:  Vector that contains the 'theta' values for the columns
%               in the accumulation array.
%  lineprm:     (Optional)
%               Parameters (rho, theta) of the lines detected. Is a N-by-2
%               matrix with each row contains the parameters (rho, theta) 
%               for a line. The definitions of 'rho' and 'theta' are as
%               the following:
%               'rho' is the perpendicular distance from the line to the
%               origin of the image. 'rho' can be negative since 'theta'
%               is constrained to [0, pi]. The unit of 'rho' is in pixels.
%               'theta' is the sweep angle from axis X (i.e. horizontal
%               axis of the image) to the direction that is perpendicular
%               to the line. The range of 'theta' is [0, pi].
%  lineseg:     (Optional)
%               Parameters (x1, x2, y1, y2) of line segments detected.
%               Lines given by (rho, theta) are infinite lines. This
%               function tries to detect line segments in the raw image
%               and output them as 'lineseg', which is a Ns-by-4 matrix
%               with each row contains the parameters (x1, x2, y1, y2)
%               that defines the two ends of a line segment.
%  dbg_label:   (Optional, for debug purpose)
%               Labelled segmentation map from the search of peaks in
%               the accumulation array
%
%  EXAMPLE #1:
%  rawimg = imread('TestHT_Chkbd1.bmp');
%  fltr4img = [1 2 3 2 1; 2 3 4 3 2; 3 4 6 4 3; 2 3 4 3 2; 1 2 3 2 1];
%  fltr4img = fltr4img / sum(fltr4img(:));
%  imgfltrd = filter2( fltr4img , rawimg );
%  tic;
%  [accum, axis_rho, axis_theta, lineprm, lineseg] = ...
%      Hough_Grd(imgfltrd, 6, 0.02);
%  toc;
%  figure(1); imagesc(axis_theta*(180/pi), axis_rho, accum); axis xy;
%  xlabel('Theta (degree)'); ylabel('Pho (pixels)');
%  title('Accumulation Array from Hough Transform');
%  figure(2); imagesc(rawimg); colormap('gray'); axis image;
%  DrawLines_2Ends(lineseg);
%  title('Raw Image with Line Segments Detected');
%
%  EXAMPLE #2:
%  rawimg = imread('TestHT_Pentagon.png');
%  fltr4img = [1 1 1 1 1; 1 2 2 2 1; 1 2 4 2 1; 1 2 2 2 1; 1 1 1 1 1];
%  fltr4img = fltr4img / sum(fltr4img(:));
%  imgfltrd = filter2( fltr4img , rawimg );
%  tic;
%  [accum, axis_rho, axis_theta, lineprm, lineseg] = ...
%      Hough_Grd(imgfltrd, 8, 0.05);
%  toc;
%  figure(1); imagesc(axis_theta*(180/pi), axis_rho, accum); axis xy;
%  xlabel('Theta (degree)'); ylabel('Pho (pixels)');
%  title('Accumulation Array from Hough Transform');
%  figure(2); imagesc(imgfltrd); colormap('gray'); axis image;
%  DrawLines_Polar(size(imgfltrd), lineprm);
%  title('Raw Image (Blurred) with Lines Detected');
%  figure(3); imagesc(rawimg); colormap('gray'); axis image;
%  DrawLines_2Ends(lineseg);
%  title('Raw Image with Line Segments Detected');
%
%  BUG REPORT:
%  This is an alpha version. Please send your bug reports, comments and
%  suggestions to pengtao@glue.umd.edu . Thanks.

%  Author:  Tao Peng
%           Department of Mechanical Engineering
%           University of Maryland, College Park, Maryland 20742, USA
%           pengtao@glue.umd.edu
%  Version: alpha       Revision: Dec. 05, 2005


%  Reserved for extension
%  accumres:	(Optional, default is [4, 1], minimum is [2, 0.5])
%               The desired resolutions in 'rho' and 'theta'. Is a
%               2-element vector in following format:
%               [resolution in 'rho' (pixels),
%               resolution in 'theta' (degrees)]

% Modified by
% Andreas Eilschou, 2010

%%%%%%%% Arguments and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Validation of arguments
if ndims(img) ~= 2 || ~isnumeric(img),
    error('Hough_Grd: ''img'' has to be 2 dimensional');
end
if ~all(size(img) >= 16),
    error('Hough_Grd: ''img'' has to be larger than 16-by-16');
end

% Parameters (default values)
prm_grdthres = 8;
prm_accumres = [4, 1];
prm_detsens = 0.08;

prm_fltraccum = true;

% Validation of arguments
vap_grdthres = 1;
if nargin > vap_grdthres,
    if isnumeric(varargin{vap_grdthres}) && ...
            varargin{vap_grdthres}(1) >= 0,
        prm_grdthres = varargin{vap_grdthres}(1);
    else
        error(['Hough_Grd: ''grdthres'' has to be ', ...
            'a non-negative number']);
    end
end
%{
vap_accumres = 3;
if nargin > vap_accumres,
    if numel(varargin{vap_accumres}) == 2 && ...
            isnumeric(varargin{vap_accumres}) && ...
            ( varargin{vap_accumres}(1) >= 2 && ...
            varargin{vap_accumres}(2) >= 0.5 ),
        prm_accumres = varargin{vap_accumres};
    else
        error(['Hough_Grd: ''accumres'' has to be a two-element ', ...
            'vector and no smaller than [2, 0.5]']);
    end
end
%}
vap_detsens = 2;
if nargin > vap_detsens,
    if isnumeric(varargin{vap_detsens}) && ...
            varargin{vap_detsens}(1) > 0 && varargin{vap_detsens}(1) < 1,
        prm_detsens = varargin{vap_detsens};
    else
        error('Hough_Grd: ''detsens'' has to be in the range (0, 1)');
    end
end

%%% MODIFIED %%%
prm_line_min = 125;
%{
vap_line_min = 3;
if nargin > vap_line_min,
    if isnumeric(varargin{vap_line_min}) && ...
            varargin{vap_line_min}(1) >= 0,
        prm_line_min = varargin{vap_line_min};
    else
        error('Hough_Grd: ''line_min'' has to be in the range (0, 1)');
    end
end
%}
%%% MODIFIED END %%%

func_compu_lineprm = ( nargout > 3 );

% Size of the accumulation array
imgsize = size(img);
coef_rhorng = [ -imgsize(2), sqrt(sum(imgsize.^2)) ];
coef_thetarng = [-pi/18, pi+pi/18];

prm_accumsize = [ ...
    round( (coef_rhorng(2)-coef_rhorng(1)) * (2/prm_accumres(1)) ) , ...
    round( (coef_thetarng(2)-coef_thetarng(1)) * (180/pi) * ...
    (2/prm_accumres(2)) ) ];

% Default filter for the accumulation array
prm_acmfltr_R = 4;
prm_acmfltr_w = [1 2 4 8];

fltr4accum = ones(2 * prm_acmfltr_R - 1) * prm_acmfltr_w(1);
for k = 2 : prm_acmfltr_R,
    fltr4accum(k:(2*prm_acmfltr_R-k), k:(2*prm_acmfltr_R-k)) = ...
        prm_acmfltr_w(k);
end
fltr4accum = fltr4accum / sum(fltr4accum(:));

% Parameters for the algorithm using local maximum filter
prm_fltrLM_R = 4;
prm_fltrLM_s = 1.3;
prm_fltrLM_r = ceil( prm_fltrLM_R * 0.6 );
prm_fltrLM_npix = 6;

% Reserved parameters
dbg_on = false;      % debug information
dbg_bfigno = 5;
%if nargout > 5,  dbg_on = true;  end


%%%%%%%% Building accumulation array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the gradient and the magnitude of gradient
img = double(img);
[grdx, grdy] = gradient(img);
grdmag = sqrt(grdx.^2 + grdy.^2);

% Clear the margins of the gradient field
prm_grdfldmgn = 4;
grdmag([1:prm_grdfldmgn, (end-prm_grdfldmgn+1):end], :) = 0;
grdmag(:, [1:prm_grdfldmgn, (end-prm_grdfldmgn+1):end]) = 0;

% Get the linear indices, as well as the subscripts, of the pixels
% whose gradient magnitudes are larger than the given threshold
grdmasklin = find(grdmag > prm_grdthres);
[grdmask_IdxI, grdmask_IdxJ] = ind2sub(imgsize, grdmasklin);

% Compute (the line parameter) 'theta' for all voting pixels
grdphs_vot = atan2( grdy(grdmasklin), grdx(grdmasklin) );

% -- Regulate 'grdphs_vot' from [-pi, pi] to [0, pi]
grdphs_vot = grdphs_vot + pi * (grdphs_vot < 0);

% Compute the 'theta'-subscript (to the accumulation array)
% for all voting pixels
coef_subtheta = prm_accumsize(2) / ...
    (coef_thetarng(2) - coef_thetarng(1)) * (1 - 1e-6);
sub_theta = ceil( (grdphs_vot - coef_thetarng(1)) * coef_subtheta );

% 'theta' vector for the accumulation array
axis_theta = (coef_thetarng(1) + 0.5 / coef_subtheta) : ...
    (1 / coef_subtheta) : coef_thetarng(2);

% Compute the 'rho' values for all voting pixels
% rho = (J - 0.5) * cos(theta) + (I - 0.5) * sin(theta)
rho_vot = (grdmask_IdxJ - 0.5) .* cos(grdphs_vot) + ...
    (grdmask_IdxI - 0.5) .* sin(grdphs_vot);

% Compute the 'rho'-subscript (to the accumulation array)
% for all voting pixels
coef_subrho = prm_accumsize(1) / ...
    (coef_rhorng(2) - coef_rhorng(1)) * (1 - 1e-6);
sub_rho = ceil( (rho_vot - coef_rhorng(1)) * coef_subrho );

% 'rho' vector for the accumulation array
axis_rho = (coef_rhorng(1) + 0.5 / coef_subrho) : ...
    (1 / coef_subrho) : coef_rhorng(2);


%%% MODIFIED START %%%

% The weight function for the accumulation array:
% w(g) = 1/A Int_0^g h(g) dg
% with h(g) the histogram of gradient magnitude values in the ROI and the
% normalization factor 1/A the area under h(g)

% Number of bins in the histogram
bins = 100;
% Number of elements in each bin of the histogram
%figure; hist(grdmag(grdmasklin),20);
[n_elements,bin_centers] = hist(grdmag(grdmasklin),bins);
% Bin width
bin_width = bin_centers(2)-bin_centers(1);
% Cumulative sum of these elements
c_elements = cumsum(n_elements);
% Divide by total area
grd_weight = c_elements ./ sum(n_elements); %sum(grdmag(grdmasklin));
% Gradient optima
miniMum = min(grdmag(grdmasklin));
maxiMum = max(grdmag(grdmasklin));
% Substitute gradient magnitude weight with the new weight
for i = grdmasklin'
    % Take care of upper boundary
    if grdmag(i) == maxiMum
        weight_index = bins;
    else
        weight_index = (round((grdmag(i)-miniMum)/bin_width + 1/2));
    end
    grdmag(i) = grd_weight(weight_index);
end
%figure; hist(grdmag(grdmasklin),20);

%%% MODIFIED END %%%


% Build the accumulation array
accum = accumarray( sub2ind(prm_accumsize, sub_rho, sub_theta), ...
    grdmag(grdmasklin) );
accum = [ accum ; ...
    zeros(prm_accumsize(1) * prm_accumsize(2) - numel(accum), 1) ];
accum = reshape( accum, prm_accumsize );


%%%%%%%% Locating peaks in the accumulation array %%%%%%%%%%%%%%%%%%%

% Stop if no need to estimate the parameters of the lines
if ~func_compu_lineprm,
    return;
end


%%% MODIFIED START %%%

%  All values with rho < 0 or 0.7 Pi < theta < 0.98 Pi are set to zero
axis_rho_mask = axis_rho < 0;
axis_theta_mask = (axis_theta < 0.03*pi) | (axis_theta > 0.45*pi);%0.2*pi);
%axis_theta_mask = (0.7*pi < axis_theta) & (axis_theta < 0.98*pi);
accum(axis_rho_mask,:) = 0;
accum(:,axis_theta_mask) = 0;

% A correction on the values taking the length of the projected line
% l(rho,theta) into account:
% H'(rho,theta) = H(rho,theta)/sqrt(L)
% with L = l for pectoral boundaries greater than prm_line_min and 
% with L = prm_line_min otherwise

% Only calculate for nonzero accum values
accum_nonzero = find(accum);
for k = accum_nonzero'
    [i,j] = ind2sub(prm_accumsize,k);
    line_rho = axis_rho(i);
    line_theta = axis_theta(j);
    
    [I1 J1 I2 J2] = polar2borderpoints(imgsize,line_rho,line_theta);

    % The length of the line
    l = sqrt((I1-I2)^2 + (J1-J2)^2);
    
    % The value correction
    if l > prm_line_min
        accum(i,j) = accum(i,j)/sqrt(l);
    else
        accum(i,j) = accum(i,j)/sqrt(prm_line_min);
    end
end

%%% MODIFIED END %%%


% Smooth the accumulation array
if prm_fltraccum,
    accum = filter2( fltr4accum, accum );
end

% Find the maximum value in the accumulation array
accum_max = max(accum(:));


%------- Algorithm 1: Local maximum filter (begin) -------------
% Build the local maximum filter
fltr4LM = zeros(2 * prm_fltrLM_R + 1);

[mesh4fLM_x, mesh4fLM_y] = meshgrid(-prm_fltrLM_R : prm_fltrLM_R);
mesh4fLM_r = sqrt( mesh4fLM_x.^2 + mesh4fLM_y.^2 );
fltr4LM_mask = ...
	( mesh4fLM_r > prm_fltrLM_r & mesh4fLM_r <= prm_fltrLM_R );
fltr4LM = fltr4LM - ...
	fltr4LM_mask * (prm_fltrLM_s / sum(fltr4LM_mask(:)));

if prm_fltrLM_R >= 4,
	fltr4LM_mask = ( mesh4fLM_r < (prm_fltrLM_r - 1) );
else
	fltr4LM_mask = ( mesh4fLM_r < prm_fltrLM_r );
end
fltr4LM = fltr4LM + fltr4LM_mask / sum(fltr4LM_mask(:));


% Whole accumulation array as the one AOI
accumAOI = [1, size(accum,1), 1, size(accum,2)];


% **** Debug code (begin)
if dbg_on && prm_useaoi,
    dbg_accumLM = zeros(size(accum));
end
% **** Debug code (end)

% For each of the AOIs selected, locate the local maxima
lineprm = zeros(0,3);
%lineprm = zeros(0,2);
for k = 1 : size(accumAOI, 1),
    aoi = accumAOI(k,:);    % just for referencing convenience

    % Apply the local maxima filter
    candLM = conv2( accum(aoi(1):aoi(2), aoi(3):aoi(4)) , ...
        fltr4LM , 'same' );

	% Thresholding of 'candLM' & 'accum'
    candLM_mask = ( candLM > (max(candLM(:))*prm_detsens)); % & ...
    %   accum(aoi(1):aoi(2),aoi(3):aoi(4)) > (accum_max*prm_detsens) );

    % Clear the margins of 'candLM_mask'
    candLM_mask([1:prm_fltrLM_R, (end-prm_fltrLM_R+1):end], :) = 0;
    candLM_mask(:, [1:prm_fltrLM_R, (end-prm_fltrLM_R+1):end]) = 0;

    % **** Debug code (begin)
    if dbg_on && prm_useaoi,
        dbg_accumLM(aoi(1):aoi(2), aoi(3):aoi(4)) = ...
            dbg_accumLM(aoi(1):aoi(2), aoi(3):aoi(4)) + candLM;
    end
    % **** Debug code (end)

    % Group the local maxima candidates by adjacency, compute the
    % centroid position for each group and take that as the center
    % of one circle detected
    [candLM_label, candLM_nRgn] = bwlabel( candLM_mask, 8 );

    max_peak_value = 0;
    for ilabel = 1 : candLM_nRgn,
        % Indices (to current AOI) of the pixels in the group
        candgrp_masklin = find( candLM_label == ilabel );
        [candgrp_IdxI, candgrp_IdxJ] = ...
            ind2sub( size(candLM_label) , candgrp_masklin );

        % Indices (to 'accum') of the pixels in the group
        candgrp_IdxI = candgrp_IdxI + ( aoi(1) - 1 );
        candgrp_IdxJ = candgrp_IdxJ + ( aoi(3) - 1 );
        candgrp_idx2acm = ...
            sub2ind( size(accum) , candgrp_IdxI , candgrp_IdxJ );

        % Minimum number of qulified pixels in the group
        if sum(numel(candgrp_masklin)) < prm_fltrLM_npix,
            continue;
        end
        
        % Compute the centroid position
        candgrp_acmsum = sum( accum(candgrp_idx2acm) );
        cc_rho = sum( candgrp_IdxI .* accum(candgrp_idx2acm) ) / ...
            candgrp_acmsum;
        cc_theta = sum( candgrp_IdxJ .* accum(candgrp_idx2acm) ) / ...
            candgrp_acmsum;
%         if max_peak_value < max(candLM(candgrp_masklin))
%            max_peak_value = max(candLM(candgrp_masklin));
%            lineprm = [cc_rho,cc_theta];
%         end
        lineprm = [lineprm; cc_rho, cc_theta, accum(round(cc_rho),round(cc_theta))];
    end
end

%output maximum peak
if size(lineprm,1) > 0
    varargout{2} = max(lineprm(:,3));
else
    varargout{2} = 0;
end


% **** Debug code (begin)
if dbg_on,
    figure(dbg_bfigno);
    if prm_useaoi,
        imagesc(dbg_accumLM); axis xy;
    else
        imagesc(candLM); axis xy;
    end
end
% **** Debug code (end)
%------- Algorithm 1: Local maximum filter (end) ---------------


% Convert 'lineprm' from pixel coordinates to (rho, theta)
lineprm(:,1:2) = [ (lineprm(:,1) - 0.5)/coef_subrho + coef_rhorng(1), ...
	(lineprm(:,2) - 0.5)/coef_subtheta + coef_thetarng(1) ];
%lineprm = [ (lineprm(:,1) - 0.5)/coef_subrho + coef_rhorng(1), ...



%%% MODIFIED START %%%

% Select the pectoral line that maximizes the product of the hough-value
% and the area of the pectoral line
lineprm_candidates = size(lineprm,1);
if lineprm_candidates > 1

lineprm = horzcat(lineprm,zeros(lineprm_candidates,1));
img_mask = logical(img);
for k = 1:size(lineprm,1)
    [I1 J1 I2 J2] = polar2borderpoints(size(img),lineprm(k,1),lineprm(k,2));
    if img_mask(floor(I1),floor(J1)) == 0 || img_mask(floor(I2),floor(J2)) == 0
        lineprm(k,3) = 0;
    end
    % Area of Pectoral
    lineprm(k,4) = abs(I1 - I2) * abs(J1 - J2) / 2;
end


% Find two highest peaks, choose the one with biggest peakvalue*area
[dummy, line_sort_index] = sort(lineprm(:,3),'descend');
lsi = line_sort_index(1:2);
[dummy, sort_max_index] = max(lineprm(lsi,3).*lineprm(lsi,4));
lineprm = lineprm(line_sort_index(sort_max_index),1:2);
elseif lineprm_candidates == 1
    lineprm = lineprm(1:2);
else
    lineprm = [];
end



%%%[~, line_max_index] = max(lineprm(:,3).*lineprm(:,4));
%%%lineprm = lineprm(line_max_index,1:2);
%[~, line_sort_index] = sort(lineprm(:,3),'descend');
%[~, sort_max_index] = max(lineprm(line_sort_index(1:2),4));
%lineprm = lineprm(line_sort_index(sort_max_index),1:2);

%{
% Select the pectoral line that maximizes the pectoral area.
% Restrict pectoral lines to not go across the background
figure;
lineprm = horzcat(lineprm,zeros(size(lineprm,1),1));
max_pect_area = [0,0];
for k = 1:size(lineprm,1)
    [area,pectoral_mask,outside_ROI] = pectoral_area(img,lineprm(k,1),lineprm(k,2));
    if ~outside_ROI && area > max_pect_area(1)
        max_pect_area = [area,k];
    end
    imshow(pectoral_mask,[])
    title(['No: ' num2str(k) ', (rho,theta) = (' num2str(lineprm(k,1)) ', ' num2str(lineprm(k,2)) ')'])
    pause
    lineprm(k,end) = area;
% subplot(1,size(lineprm,1)+1,k+1); imagesc(pect_area_mask); title(num2str(area));
end


if max_pect_area > 0 %(size(lineprm,1) > 0)
    lineprm = lineprm(max_pect_area(2),:);
else
    lineprm = [];
end
%}

%%% MODIFIED END %%%


% Output 'lineprm'
varargout{1} = lineprm;
if nargout <= 4,
    return;
end

%{
%%%%%%%% Locating the line segments in the raw image %%%%%%%%%%%%%%%%

% Parameters for locating the line segments
prm_ls_tTol = [-pi/7.5, pi/7.5];
prm_ls_pTol = [-3.5, 3.5];
prm_ls_minlen = 10;

% Locate the two ends for all lines detected
rgnmask = logical(zeros(imgsize));
lineseg = zeros(0, 4);

for k = 1 : size(lineprm, 1),
    % Compute the 'rho' values for all pixels voted,
    % assuming they contribute to the current line
    rho2_vot = (grdmask_IdxJ - 0.5) * cos( lineprm(k,2) ) + ...
        (grdmask_IdxI - 0.5) .* sin( lineprm(k,2) );

    % Find the pixels that belong to the line
    PixOnLn_votmask = ( grdphs_vot > (lineprm(k,2) + prm_ls_tTol(1)) & ...
        grdphs_vot < (lineprm(k,2) + prm_ls_tTol(2)) & ...
        rho2_vot > (lineprm(k,1) + prm_ls_pTol(1)) & ...
        rho2_vot < (lineprm(k,1) + prm_ls_pTol(2)) );

	PixOnLn_lin = grdmasklin( find(PixOnLn_votmask) );
    if isempty(PixOnLn_lin),
        continue;
    end
    [PixOnLn_IdxI, PixOnLn_IdxJ] = ind2sub(imgsize, PixOnLn_lin);

    % Find the axis-aligned bounding box for these pixels
    bndbox = [ min(PixOnLn_IdxI), max(PixOnLn_IdxI), ...
        min(PixOnLn_IdxJ), max(PixOnLn_IdxJ) ];

    % Seperate the line segments by adjacencies
    rgnmask( bndbox(1):bndbox(2), bndbox(3):bndbox(4) ) = 0;
    rgnmask( PixOnLn_lin ) = 1;
    [rgnlabel, nrgn] = bwlabel( ...
        rgnmask(bndbox(1):bndbox(2), bndbox(3):bndbox(4)), 8 );

    for k_rgn = 1 : nrgn,
        % Get the linear indices of pixels that belong to one segment
        seg_lin = find( rgnlabel == k_rgn );
        [seg_IdxI, seg_IdxJ] = ind2sub( size(rgnlabel), seg_lin );

        % Ignore if the line segment is too short
        segspan = [ min(seg_IdxI) + bndbox(1) - 1, ...
            max(seg_IdxI) + bndbox(1) - 1, ...
            min(seg_IdxJ) + bndbox(3) - 1, ...
            max(seg_IdxJ) + bndbox(3) - 1 ];
        if (segspan(2) - segspan(1) + 1) < prm_ls_minlen && ...
                (segspan(4) - segspan(3) + 1) < prm_ls_minlen,
            continue;
        end

        % Get the [x1 x2 y1 y2] structure for the line SEGMENT
        if lineprm(k,2) > pi/4 && lineprm(k,2) < 3*pi/4,
            % ls_base_x = 0;
            ls_base_y = lineprm(k,1) / sin(lineprm(k,2));
            lineseg = [ lineseg; segspan(3) - 0.5, segspan(4) - 0.5, ...
                ls_base_y - (segspan(3) - 0.5) / tan(lineprm(k,2)), ...
                ls_base_y - (segspan(4) - 0.5) / tan(lineprm(k,2)) ];
        else
            % ls_base_y = 0;
            ls_base_x = lineprm(k,1) / cos(lineprm(k,2));
            lineseg = [ lineseg; ...
                ls_base_x - (segspan(1) - 0.5) * tan(lineprm(k,2)), ...
                ls_base_x - (segspan(2) - 0.5) * tan(lineprm(k,2)), ...
                segspan(1) - 0.5, segspan(2) - 0.5 ];
        end
    end
end

% Output 'lineseg'
varargout{2} = lineseg;
if nargout > 5,
    varargout{3} = dbg_label;
end
%}
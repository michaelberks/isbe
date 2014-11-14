function [brief_codes, x_pts, y_pts] = compute_brief(image_in, f_rows, f_cols, sigma, x_pts, y_pts, S, nd, method, debug_flag)
%COMPUTE_BRIEF *Insert a one line summary here*
%   [brief_codes, x_pts, y_pts] = compute_brief(image_in, f_rows, f_cols, sigma, x_pts, y_pts, S, nd, method)
%
% Inputs:
%      image_in - Image to extract BRIEF codes from
%
%      f_rows - Nx1 vector of row subscripts for feature points
%
%      f_cols - Nx1 vector of column subscripts for feature points
%
%      sigma - SD of Gaussian kernel used to smooth image prior to sampling
%
%      x_pts - 2 x nd vector of sampling locations within an SxS grid, 1st
%      row of x_pts is assumed to be rows, 2nd row assumed to be cols
%
%      y_pts - 2 x nd vector of sampling locations within an SxS grid, 1st
%      row of x_pts is assumed to be rows, 2nd row assumed to be cols
%
%      *optional arguments if x_pts/y_pts not specified*
%
%      S - the size of patch in which to comput x_pts/y_pts
%
%      nd - number of sample pts
%
%      method - {1,2,3} determines how x_pts/y_pts are randomly sampled
%      (see paper for more details)
%
%
% Outputs:
%      brief_codes - N x nd matrix, each row of which is the BRIEF codes
%      for a particular feature point
%
%      x_pts/y_pts - returned if sampling for the first time
%
%
% Example:
%   [brief_codes, x_pts, y_pts] = compute_brief(image_in, f_rows, f_cols, 2, x_pts, y_pts)
%   [brief_codes, x_pts, y_pts] = compute_brief(image_in, f_rows, f_cols, 2, [], [], 32, 256, 1)
%
% Notes: Based on the paper "BRIEF: Binary Robust Independent Elementary
% Features" by Colander, Lepetit, Fua, ECCV 2010
%
% See also:
%
% Created: 30-Sep-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if nargin < 10
    debug_flag = 0;
end

%Generate coordinates of sampling points if we've not alreday been given
%these
if isempty(x_pts) || isempty(y_pts)
    
    switch method %Methods GI-III from paper (IV and V omitted)
        
        case 1 %Uniform sampling from SxS patch
            x_pts = round((rand(2, nd) * S) - 0.5*S);
            y_pts = round((rand(2, nd) * S) - 0.5*S);
        case 2 %Gaussian sampling X,Y~N(0, S^2/25)
            x_pts = round(randn(2, nd) * S/5);
            y_pts = round(randn(2, nd) * S/5);
            x_pts(x_pts > S/2) = S/2;
            y_pts(y_pts > S/2) = S/2;
            x_pts(x_pts < -S/2) = -S/2;
            y_pts(y_pts < -S/2) = -S/2;
        case 3 %Local Gaussian sampling X~N(0, S^2/25), Y~(xi, S^2/100)
            x_pts = round(randn(2, nd) * S/5);
            y_pts = round(x_pts + randn(2, nd) * S/10);
            x_pts(x_pts > S/2) = S/2;
            y_pts(y_pts > S/2) = S/2;
            x_pts(x_pts < -S/2) = -S/2;
            y_pts(y_pts < -S/2) = -S/2;
        otherwise %Default method 2
            x_pts = round(randn(2, nd) * S/5);
            y_pts = round(randn(2, nd) * S/5);
            x_pts(x_pts > S/2) = S/2;
            y_pts(y_pts > S/2) = S/2;
            x_pts(x_pts < -S/2) = -S/2;
            y_pts(y_pts < -S/2) = -S/2;
    end
    
    if debug_flag
        figure; plot([x_pts(1,:); y_pts(1,:)], [x_pts(2,:); y_pts(2,:)]); 
        axis equal; axis(0.5*[-S S -S S]);
    end
end

%Make sure f_rows and f_cols are column vectors
f_rows = f_rows(:);
f_cols = f_cols(:);

%Compute row/col subscripts for all the x/y sampling points
f_rows_x = bsxfun(@plus, f_rows, x_pts(1,:));
f_cols_x = bsxfun(@plus, f_cols, x_pts(2,:));
f_rows_y = bsxfun(@plus, f_rows, y_pts(1,:));
f_cols_y = bsxfun(@plus, f_cols, y_pts(2,:));

%Convert subscriupts to indices to allow easy sampling from image -  - note we've not checked
%whether these pts are all sampled from within the image, could force
%this now other sub2ind might error...
f_idx_x = sub2ind(size(image_in), f_rows_x, f_cols_x);
f_idx_y = sub2ind(size(image_in), f_rows_y, f_cols_y);

%Smooth image (note method calls for smoothing patches but easier to do it
%this way)
image_in = imfilter(image_in, fspecial('gaussian', 9, sigma), 'replicate');

%Sample pixel intensities at x and y pts
p_x = image_in(f_idx_x);
p_y = image_in(f_idx_y);

%Compute binary BRIEF codes
brief_codes = p_x < p_y;

if debug_flag
    figure; imagesc(image_in); axis image; colormap(gray(256)); hold on;
    plot(f_cols, f_rows, 'rx');
    plot([f_cols_x(:) f_cols_y(:)]', [f_rows_x(:) f_rows_y(:)]');
end
    


            
            
function [image_out, label, label_centre] = create_ellipse_spicule(halfwidth, contrast, orientation, row, col)
%
%CREATE_ELLIPSE_SPICULE create an image containing a spicule
%
% USAGE:
%   [image_out, label, label_centre] = create_ellipse_spicule(halfwidth,
%   contrast, orientation, row, col)
%
% Inputs:
%      halfwidth - halfwidth of Gaussian profile at half its maximum height
%
%      contrast - maximum height of (scaled) Gaussian profile
%
%      orientation - orientation of bar in image in degrees
%
%      row - number of rows in image
%
%      col - number of columns in image
%
%
% Outputs:
%      image_out - image containing spicule
%
%      label - label of spicule (1) vs background (0) 
%
%      label_centre - the label of centre line (1) vs background (0)
%
% Example:
%
% Notes:
%
% Created: 04-February-2010
% Author: Zezhi Chen
% Email : zezhi.chen@manchester.ac.uk
% Phone : +44 (0)161 275 7669
% Copyright: (C) University of Manchester
% 

% create a square image
max_rl=max(row, col);
rw=ceil(sqrt(2*max_rl*max_rl)/2)*2;
cl=rw;

temp_im=zeros(rw, cl);
temp_label=temp_im;
tlabel=temp_im;
temp_label_centre=temp_label;
for k=1:cl
    curr_radi = round(halfwidth*(cl-k)/cl);
    curr_contrast = contrast*(cl-k)/cl;
    dx = -curr_radi+1:curr_radi-1;
    ex = (curr_contrast/curr_radi)*sqrt(curr_radi*curr_radi-dx.*dx);
    if length(dx)>3
        temp_im(round(rw/2)+dx, k)=ex;
        temp_label(round(rw/2)+dx, k)=1;
        temp_label_centre(round(rw/2), k)=1;
    end
end

H = fspecial('average', 7); % 7x7 window
temp_im = imfilter(temp_im,H);

%rotate the image
temp_label= imrotate(imfilter(temp_label,H), orientation, 'bilinear','crop');
temp_im = imrotate(temp_im, orientation, 'bilinear','crop');
temp_label_centre = imrotate(temp_label_centre, orientation, 'nearest','crop');

tlabel(temp_label>0.5)=1;

% Crop the image as the right size
pad_row=ceil((rw-row)/2);
pad_col=ceil((cl-col)/2);
label=tlabel(pad_row+1:rw-pad_row, pad_col+1:cl-pad_col);
image_out=temp_im(pad_row+1:rw-pad_row, pad_col+1:cl-pad_col);
label_centre = temp_label_centre(pad_row+1:rw-pad_row, pad_col+1:cl-pad_col);

return
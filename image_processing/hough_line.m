function [line_scores rhos thetas] = hough_line(bw, rhos, thetas)
%HOUGH_LINE given a binary edge map compute hough transform to determine
%votes for lines in the image
%
% Inputs:
%      bw - a binary edge map of the image
%
%      rhos - a 1D vector of possible rho values
%
%      thetas - a 1D vector of possible rho values
%
% Outputs:
%
%      line_scores - an array of votes where each row corresponds to a
%      possible centre and each column to a radius given that centre
%
%
%
% Example:
%
% Notes:
%
% See also: HOUGH
%
% Created: 11-Oct-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Get image size
[rows cols] = size(bw);

%Set defaults for missing arguments - can't be used to a full
%name/parameter search
if nargin < 2 
    %set centres to all possible (x,y) locations in image
    rhos = -ceil(norm([rows cols])):ceil(norm([rows cols]));
end
if nargin < 3
    %set default minimum radius to 5
    thetas = linspace(-pi/2, pi/2, 180);
end
    
%Find all edge pts in BW
[ey ex] = find(bw);

%Create table to store votes in
line_scores = zeros(length(rhos), length(thetas)); 
min_rho = min(rhos);

%For each edge pt
for ii = 1:length(ex)

    for jj = 1:length(thetas)
        th = thetas(jj);
        rho = -ex(ii)*sin(th) + ey(ii)*cos(th);
        
        line_scores(ceil(rho)-min_rho+1, jj) = line_scores(ceil(rho)-min_rho+1, jj) + 1;
        line_scores(floor(rho)-min_rho+1, jj) = line_scores(floor(rho)-min_rho+1, jj) + 1;
        %line_scores(round(rho)-min_rho+1, jj) = line_scores(round(rho)-min_rho+1, jj) + 1;
    end
end

        
        